# AI-assisted (Claude), exploration code for CS581 project
"""Alignment error measures against the true alignment, including tree-oriented ones.

    python3 measures.py [--lanes 4]   -> data/measures.jsonl (one row per (key, method) in data/aln_paths.jsonl,
                                         plus 'true'; restartable)

Every residue r (taxon, position in the ungapped sequence) has a true column T(r) and an estimated column E(r)
(-1 if a masked alignment dropped it). n[c, e] = residues in true column c and estimated column e.
  SPFN, SPFP          1 - shared pairs / true pairs, 1 - shared pairs / estimated pairs (= FastSP)
  TC                  true columns (>= 2 residues) reproduced exactly (column recall)
  col_prec            estimated columns (>= 2 residues) equal to a true column (column precision)
  oversplit_cols      true columns (>= 2 res) whose residues occupy > 1 estimated column
  overmerge_cols      estimated columns (>= 2 res) holding residues of > 1 true column
  split_res           residues outside the plurality estimated column of their true column  (residue level)
  misplaced_res       residues outside the plurality true column of their estimated column (residue level)
  *_pi                the same restricted to parsimony-informative columns (true PI columns for FN-type
                      measures, estimated PI columns for FP-type ones); SPFN_pi, SPFP_pi; pi_cols_ratio
  gap_open_ratio      internal gap openings (per row, est / true): gap-placement volume
  bnd_err / int_err   misplaced-or-split rate of residues next to a true indel (boundary) vs interior ones
  len_ratio           estimated / true length
  fitch_excess        (Fitch parsimony length of the alignment on the TRUE tree - that of the true alignment)
                      / that of the true alignment: homoplasy that alignment error adds on the true tree
  fitch_excess_pi     the same over PI columns only
  fp_clade            among false pairs, the share whose two taxa are each other's nearest neighbours in the
                      true tree (cheap, clade-local) vs spread across the tree: here as the mean true-tree
                      path length (edges) between the taxa of false pairs, sampled; and of missed pairs
  tax_err_p90         90th percentile over taxa of the per-taxon misplaced+split residue rate
"""
import argparse
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
W = "/opt/work/treecrit/reps"
GAP = set("-.")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def tree_file(rep):
    for f in ("true_tree.nwk", "tree.nwk", "true.tree"):
        if os.path.exists(os.path.join(rep, f)):
            return os.path.join(rep, f)
    raise FileNotFoundError(rep)


def parse_newick(s):
    """Minimal Newick parser -> (children list, leaf name per node, root). Branch lengths ignored."""
    s = s.strip().rstrip(";")
    children, names, parent = [[]], [None], [-1]
    stack, cur, i = [], 0, 0
    tok = ""
    while i < len(s):
        ch = s[i]
        if ch == "(":
            children.append([]); names.append(None); parent.append(cur)
            children[cur].append(len(children) - 1)
            stack.append(cur); cur = len(children) - 1
            i += 1
        elif ch in ",)":
            if tok:
                nm = tok.split(":")[0].strip().strip("'")
                if nm and names[cur] is None and not children[cur]:
                    names[cur] = nm
            tok = ""
            if ch == ",":
                p = parent[cur]
                children.append([]); names.append(None); parent.append(p)
                children[p].append(len(children) - 1)
                cur = len(children) - 1
            else:
                cur = parent[cur]
            i += 1
        else:
            tok += ch
            i += 1
    if tok:
        pass
    return children, names, parent


class Tree:
    def __init__(self, path, taxa):
        txt = open(path).read()
        if txt.lstrip().startswith("[&"):
            txt = txt.split("]", 1)[1]
        ch, nm, par = parse_newick(txt)
        # The first node is a placeholder root; if it only has one child, use the child.
        self.children, self.names, self.parent = ch, nm, par
        order, st = [], [0]
        while st:
            v = st.pop(); order.append(v); st.extend(ch[v])
        self.post = order[::-1]
        self.leaf = {nm[v]: v for v in range(len(ch)) if not ch[v]}
        missing = [t for t in taxa if t not in self.leaf]
        assert not missing, missing[:3]
        depth = np.zeros(len(ch), dtype=np.int64)
        for v in order:
            if par[v] >= 0:
                depth[v] = depth[par[v]] + 1
        self.depth = depth

    def pathlen(self, a, b):
        """edges between nodes a and b (vectorised over arrays via climbing)."""
        a, b = np.array(a), np.array(b)
        par = np.array(self.parent)
        da, db = self.depth[a].copy(), self.depth[b].copy()
        d = np.zeros(len(a), dtype=np.int64)
        while True:
            m = da > db
            if m.any():
                a[m] = par[a[m]]; da[m] -= 1; d[m] += 1; continue
            m = db > da
            if m.any():
                b[m] = par[b[m]]; db[m] -= 1; d[m] += 1; continue
            m = a != b
            if not m.any():
                return d
            a[m] = par[a[m]]; b[m] = par[b[m]]; da[m] -= 1; db[m] -= 1; d[m] += 2

    def fitch(self, cols, taxa):
        """cols: (ntaxa, L) uint32 state bitmasks (0 = gap -> all states). Returns per-column Fitch length."""
        L = cols.shape[1]
        allm = np.uint32(0xFFFFFFFF)
        state = {}
        idx = {t: i for i, t in enumerate(taxa)}
        cost = np.zeros(L, dtype=np.int64)
        for v in self.post:
            if not self.children[v]:
                x = cols[idx[self.names[v]]]
                state[v] = np.where(x == 0, allm, x)
            else:
                kids = self.children[v]
                s = state.pop(kids[0])
                for k in kids[1:]:
                    t = state.pop(k)
                    inter = s & t
                    empty = inter == 0
                    cost += empty
                    s = np.where(empty, s | t, inter)
                state[v] = s
        return cost


def residue_cols(aln, taxa):
    """per taxon: array of column index for each residue; and the residue characters."""
    out = {}
    for t in taxa:
        s = aln.get(t)
        if s is None:
            out[t] = None
            continue
        a = np.frombuffer(s.encode(), dtype=np.uint8)
        g = (a == ord("-")) | (a == ord("."))
        out[t] = np.nonzero(~g)[0]
    return out


def state_matrix(aln, taxa, alphabet):
    L = len(next(iter(aln.values())))
    lut = np.zeros(256, dtype=np.uint32)
    for i, c in enumerate(alphabet):
        lut[ord(c)] = np.uint32(1) << np.uint32(i)
        lut[ord(c.lower())] = np.uint32(1) << np.uint32(i)
    M = np.zeros((len(taxa), L), dtype=np.uint32)
    for i, t in enumerate(taxa):
        M[i] = lut[np.frombuffer(aln[t].upper().encode(), dtype=np.uint8)]
    return M


def pi_columns(M):
    """parsimony-informative: >= 2 states each occurring in >= 2 taxa (single-state cells only)."""
    L = M.shape[1]
    cnt = np.zeros((32, L), dtype=np.int32)
    for b in range(32):
        cnt[b] = ((M >> np.uint32(b)) & np.uint32(1)).sum(axis=0)
    return (cnt >= 2).sum(axis=0) >= 2


def c2(x):
    x = x.astype(np.float64)
    return (x * (x - 1) / 2).sum()


def measures(key, method, path):
    rep = os.path.join(W, key)
    true = read_fasta(os.path.join(rep, "true.fasta"))
    est = read_fasta(path)
    taxa = sorted(true)
    dna = key.startswith("1000M") or key.startswith("RNA") or key.startswith("16S")
    alphabet = "ACGT" if dna else "ACDEFGHIKLMNPQRSTVWY"
    rt, re_ = residue_cols(true, taxa), residue_cols(est, taxa)
    T, E, TX, BND = [], [], [], []
    for ti, t in enumerate(taxa):
        a = rt[t]
        b = re_[t]
        n = len(a)
        if b is None or len(b) != n:
            # masked alignment: residues of dropped columns are missing; recover by sequence comparison
            sa = true[t].replace("-", "").replace(".", "")
            sb = est.get(t, "")
            pos = np.full(n, -1, dtype=np.int64)
            if sb:
                arr = np.frombuffer(sb.encode(), dtype=np.uint8)
                keep = np.nonzero((arr != ord("-")) & (arr != ord(".")))[0]
                # the masked file is a column subset of the unmasked one; map via the unmasked out.fasta
                full = read_fasta(path.replace(".masked.fasta", ".fasta"))[t]
                fa = np.frombuffer(full.encode(), dtype=np.uint8)
                fres = np.nonzero((fa != ord("-")) & (fa != ord(".")))[0]
                kept_cols = KEPT.get(path)
                if kept_cols is None:
                    raise RuntimeError("masked map missing")
                newidx = np.full(len(fa), -1, dtype=np.int64)
                newidx[kept_cols] = np.arange(len(kept_cols))
                pos = newidx[fres]
            b = pos
        T.append(a); E.append(b); TX.append(np.full(n, ti))
        d = np.diff(a)
        bnd = np.zeros(n, dtype=bool)
        if n > 1:
            jump = d > 1
            bnd[:-1] |= jump
            bnd[1:] |= jump
        BND.append(bnd)
    T, E, TX, BND = map(np.concatenate, (T, E, TX, BND))
    N = len(T)
    present = E >= 0
    Ltrue = len(next(iter(true.values())))
    Lest = len(next(iter(est.values())))
    # contingency
    code = T[present].astype(np.int64) * (Lest + 1) + E[present]
    u, n_ce = np.unique(code, return_counts=True)
    uc, ue = u // (Lest + 1), u % (Lest + 1)
    tcount = np.bincount(T, minlength=Ltrue)
    ecount = np.bincount(E[present], minlength=Lest)
    shared = c2(n_ce)
    ptrue, pest = c2(tcount), c2(ecount)
    r = {"key": key, "method": method, "SPFN": 1 - shared / ptrue, "SPFP": 1 - shared / pest if pest else 0.0}
    r["avgErr"] = (r["SPFN"] + r["SPFP"]) / 2
    # per true column: number of est columns, max cell; per est column: number of true columns, max cell
    nest_per_t = np.bincount(uc, minlength=Ltrue)
    ntrue_per_e = np.bincount(ue, minlength=Lest)
    maxc_t = np.zeros(Ltrue, dtype=np.int64); np.maximum.at(maxc_t, uc, n_ce)
    maxc_e = np.zeros(Lest, dtype=np.int64); np.maximum.at(maxc_e, ue, n_ce)
    tpres = np.bincount(T[present], minlength=Ltrue)
    multi_t = tcount >= 2
    multi_e = ecount >= 2
    exact_t = (nest_per_t == 1) & (tpres == tcount)
    # a true column is reproduced exactly if all its residues sit in one est column holding nothing else
    ecol_of_t = np.full(Ltrue, -1); ecol_of_t[uc] = ue  # valid where nest_per_t == 1
    exact_t &= ecount[np.maximum(ecol_of_t, 0)] == tcount
    r["TC"] = float(exact_t[multi_t].mean())
    tcol_of_e = np.full(Lest, -1); tcol_of_e[ue] = uc
    exact_e = (ntrue_per_e == 1) & (tcount[np.maximum(tcol_of_e, 0)] == ecount)
    r["col_prec"] = float(exact_e[multi_e].mean())
    r["oversplit_cols"] = float(((nest_per_t > 1) | (tpres < tcount))[multi_t].mean())
    r["overmerge_cols"] = float((ntrue_per_e > 1)[multi_e].mean())
    r["split_res"] = float((tcount - maxc_t).sum() / N)
    r["misplaced_res"] = float((ecount - maxc_e).sum() / N)
    # PI restriction
    Mt = state_matrix(true, taxa, alphabet)
    Me = state_matrix(est, taxa, alphabet)
    pit, pie = pi_columns(Mt), pi_columns(Me)
    sel = pit[uc]
    r["SPFN_pi"] = 1 - c2(n_ce[sel]) / c2(tcount[pit])
    sel = pie[ue]
    r["SPFP_pi"] = 1 - c2(n_ce[sel]) / c2(ecount[pie]) if pie.any() else 0.0
    r["split_res_pi"] = float((tcount - maxc_t)[pit].sum() / tcount[pit].sum())
    r["misplaced_res_pi"] = float((ecount - maxc_e)[pie].sum() / ecount[pie].sum())
    r["oversplit_cols_pi"] = float(((nest_per_t > 1) | (tpres < tcount))[pit].mean())
    r["overmerge_cols_pi"] = float((ntrue_per_e > 1)[pie].mean())
    r["pi_cols_true"], r["pi_cols_est"] = int(pit.sum()), int(pie.sum())
    r["pi_cols_ratio"] = r["pi_cols_est"] / r["pi_cols_true"]
    r["len_ratio"] = Lest / Ltrue
    r["LenEst"], r["LenRef"] = Lest, Ltrue
    # residue-level error flag: residue in the plurality cell of both its true and its est column = ok
    cell_ok_t = np.zeros(len(u), dtype=bool)
    # plurality est column of each true column (first max)
    best_e = np.full(Ltrue, -1)
    order = np.lexsort((-n_ce, uc))
    first = np.ones(len(order), dtype=bool); first[1:] = uc[order][1:] != uc[order][:-1]
    best_e[uc[order][first]] = ue[order][first]
    best_t = np.full(Lest, -1)
    order = np.lexsort((-n_ce, ue))
    first = np.ones(len(order), dtype=bool); first[1:] = ue[order][1:] != ue[order][:-1]
    best_t[ue[order][first]] = uc[order][first]
    err = np.ones(N, dtype=bool)
    ok = present.copy()
    ok[present] = (best_e[T[present]] == E[present]) & (best_t[E[present]] == T[present])
    err = ~ok
    r["res_err"] = float(err.mean())
    r["bnd_frac"] = float(BND.mean())
    r["bnd_err"] = float(err[BND].mean()) if BND.any() else 0.0
    r["int_err"] = float(err[~BND].mean())
    r["bnd_share_of_err"] = float(BND[err].mean()) if err.any() else 0.0
    per_tax = np.bincount(TX[err], minlength=len(taxa)) / np.maximum(np.bincount(TX, minlength=len(taxa)), 1)
    r["tax_err_p90"] = float(np.percentile(per_tax, 90))
    r["tax_err_gini"] = float(gini(per_tax))
    # gap openings (internal)
    def gap_opens(aln):
        tot = 0
        for t in taxa:
            a = np.frombuffer(aln[t].encode(), dtype=np.uint8)
            g = (a == ord("-")) | (a == ord("."))
            res = np.nonzero(~g)[0]
            if len(res) < 2:
                continue
            gg = g[res[0]:res[-1] + 1]
            tot += int(((~gg[:-1]) & gg[1:]).sum())
        return tot
    go_t, go_e = gap_opens(true), gap_opens(est)
    r["gap_open_ratio"] = go_e / go_t if go_t else float("nan")
    # Fitch on the true tree
    tree = Tree(tree_file(rep), taxa)
    ft = tree.fitch(Mt, taxa)
    fe = tree.fitch(Me, taxa)
    r["fitch_true"], r["fitch_est"] = int(ft.sum()), int(fe.sum())
    r["fitch_excess"] = (fe.sum() - ft.sum()) / ft.sum()
    r["fitch_excess_pi"] = (fe[pie].sum() - ft[pit].sum()) / ft[pit].sum()
    # tree distance between the taxa of false pairs and of missed pairs (sampled residue pairs)
    rng = np.random.default_rng(1)
    leaf = np.array([tree.leaf[t] for t in taxa])
    def sample_pairs(colid, mask, k=20000):
        idx = np.nonzero(mask)[0]
        if len(idx) < 2:
            return None
        o = idx[np.argsort(colid[idx], kind="stable")]
        cid = colid[o]
        starts = np.r_[0, np.nonzero(np.diff(cid))[0] + 1]
        sizes = np.diff(np.r_[starts, len(o)])
        w = sizes * (sizes - 1) / 2.0
        if w.sum() == 0:
            return None
        g = rng.choice(len(sizes), size=k, p=w / w.sum())
        i = rng.integers(0, sizes[g]); j = rng.integers(0, sizes[g] - 1); j = j + (j >= i)
        return o[starts[g] + i], o[starts[g] + j]
    # false pairs: same est column, different true column
    fp = fn = tp = None
    pr = sample_pairs(np.where(present, E, -1), present, 60000)
    if pr is not None:
        a, b = pr
        f = T[a] != T[b]
        if f.any():
            fp = tree.pathlen(leaf[TX[a][f]], leaf[TX[b][f]]).mean()
        tp = tree.pathlen(leaf[TX[a][~f]], leaf[TX[b][~f]]).mean()
    pr = sample_pairs(T, np.ones(N, dtype=bool), 60000)
    if pr is not None:
        a, b = pr
        f = ~((E[a] == E[b]) & (E[a] >= 0))
        if f.any():
            fn = tree.pathlen(leaf[TX[a][f]], leaf[TX[b][f]]).mean()
    r["fp_pathlen"], r["fn_pathlen"], r["tp_pathlen"] = (float(x) if x is not None else None for x in (fp, fn, tp))
    return r


def gini(x):
    x = np.sort(np.asarray(x, dtype=float))
    if x.sum() == 0:
        return 0.0
    n = len(x)
    return float((2 * np.arange(1, n + 1) - n - 1).dot(x) / (n * x.sum()))


KEPT = {}


def kept_columns(masked_path):
    """columns of out.fasta kept in out.masked.fasta (mask_columns removes whole columns)."""
    full = read_fasta(masked_path.replace(".masked.fasta", ".fasta"))
    msk = read_fasta(masked_path)
    taxa = sorted(full)
    F = np.array([np.frombuffer(full[t].encode(), dtype=np.uint8) for t in taxa])
    M = np.array([np.frombuffer(msk[t].encode(), dtype=np.uint8) for t in taxa])
    keep, j = [], 0
    for i in range(F.shape[1]):
        if j < M.shape[1] and (F[:, i] == M[:, j]).all():
            keep.append(i); j += 1
    assert j == M.shape[1], (masked_path, j, M.shape)
    return np.array(keep)


def run(arg):
    key, method, path = arg
    try:
        if path.endswith(".masked.fasta"):
            KEPT[path] = kept_columns(path)
        return measures(key, method, path)
    except Exception as e:  # noqa: BLE001
        import traceback
        return {"key": key, "method": method, "error": repr(e), "tb": traceback.format_exc()[-1500:]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lanes", type=int, default=4)
    ap.add_argument("--paths", default=os.path.join(DATA, "aln_paths.jsonl"))
    ap.add_argument("--out", default=os.path.join(DATA, "measures.jsonl"))
    ap.add_argument("--keys")
    a = ap.parse_args()
    jobs = [json.loads(l) for l in open(a.paths)]
    jobs = [(j["key"], j["method"], j["path"]) for j in jobs]
    keys = sorted({k for k, _, _ in jobs})
    jobs += [(k, "true", os.path.join(W, k, "true.fasta")) for k in keys]
    if a.keys:
        jobs = [j for j in jobs if j[0] in a.keys.split(",")]
    done = set()
    if os.path.exists(a.out):
        for l in open(a.out):
            r = json.loads(l)
            if "error" not in r:
                done.add((r["key"], r["method"]))
    jobs = [j for j in jobs if (j[0], j[1]) not in done]
    print(len(jobs), "jobs", flush=True)
    with Pool(a.lanes) as p, open(a.out, "a") as f:
        for r in p.imap_unordered(run, jobs):
            if "error" in r:
                print("ERR", r["key"], r["method"], r["error"], r["tb"], flush=True)
                continue
            f.write(json.dumps(r) + "\n"); f.flush()
            print(r["key"], r["method"], round(r["SPFN"], 4), round(r["SPFP"], 4), flush=True)


if __name__ == "__main__":
    main()
