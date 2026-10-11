# AI-assisted (Claude), exploration code for CS581 project
"""Merge-only runs for gcmwhy: vote.py variants plus reference-ORACLE decompositions, scored and dissected.

    python3 merge.py REP VARIANT [VARIANT ...] [--bb DIR]   -> REP/why/<VARIANT>/..., REP/why_results.jsonl

VARIANT: anything vote.py accepts (magus, es4, hard-bb, ...) or
  or-<F>-false   apply filter F (es4 | hbb) but delete ONLY the edges it removes that are false (oracle)
  or-<F>-true    apply filter F but delete ONLY the edges it removes that are true (oracle)
  or-<F>-falsenear / -falsefar   only the false edges F deletes whose nearest true partner column is <= 2 / > 2 away
  or-allfalse    delete every false cross-subset edge (oracle upper bound of graph cleaning)
  <V>@f<X>       variant V with MCL inflation factor X instead of MAGUS's default 4 (H5)
  pool<K>[-d<D>] exploratory near-miss pooling: keep edge (a, b) if the number of backbones that align a with
                 some b' in [b-D, b+D] or some a' in [a-D, a+D] with b is >= K  (default D = 1)
  poolbb[-d<D>]  near-miss pooling then the beta-binomial vote (hard-bb) on the pooled support
An edge is true if > 50 % of its backbone residue-pair units are reference homologies (calib.py / gcmgen).

For every run it records FastSP (SPFN/SPFP), cross-subset error decomposition (missed true residue pairs
split by what the filtered graph had between the two nodes) and MCL cluster statistics (H4).
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "gcmvote", "code"))
sys.path.insert(0, os.path.join(ROOT, "bbevidence", "code"))
sys.path.insert(0, os.path.join(ROOT, "code"))
sys.path.insert(0, HERE)
import vote  # noqa: E402

ORACLE = ("or-", "pool")


def bbe_node_map(rep, order):
    """vote node id (MAGUS's subset order) -> bbe node id (numeric subset order)."""
    import bbe
    files = [f for f, _ in bbe.subsets(rep)]
    lens = {}
    for f, a in bbe.subsets(rep):
        lens[f] = len(next(iter(a.values())))
    off, o = {}, 0
    for f in files:
        off[f] = o
        o += lens[f]
    return np.concatenate([off[os.path.basename(p)] + np.arange(lens[os.path.basename(p)]) for p in order])


def truth(rep, bbdir):
    """bbe-space edge keys with units and true units (cached)."""
    cache = os.path.join(rep, "why_truth_{}.npz".format(os.path.basename(bbdir.rstrip("/"))))
    if os.path.exists(cache):
        z = np.load(cache)
        return z["keys"], z["w"], z["wt"], int(z["N"])
    import bbe
    from gcmx import fasta
    R = bbe.Rep(rep)
    files = [(f, fasta.upper(fasta.read(os.path.join(bbdir, f)))) for f in sorted(os.listdir(bbdir))]
    keys, w, wt, _ = bbe.graph_edges(R, files)
    np.savez(cache, keys=keys, w=w, wt=wt, N=R.nnodes)
    return keys, w, wt, R.nnodes


def false_offsets(rep, a, b):
    """bbe node pairs -> distance to the nearest truly paired column (min over both directions), as why.py."""
    import bbe
    import why
    R = bbe.Rep(rep)
    N = R.nnodes
    R.node_sub = np.empty(N, dtype=np.int64)
    lo, hi, o = [], [], 0
    for i, s_ in enumerate(R.subsets):
        L = len(next(iter(s_.values())))
        R.node_sub[o:o + L] = i
        lo.append(o); hi.append(o + L)
        o += L
    R.sub_lo, R.sub_hi = np.array(lo), np.array(hi)
    tk = why.true_pair_keys(R)
    return np.minimum(why.nearest_offset(R, tk, a, b), why.nearest_offset(R, tk, b, a))


def pooled(total_per_bb, r, c, sub, d):
    """Backbones aligning a with b' in [b-d, b+d] (same subset) or a' in [a-d, a+d] with b."""
    n = len(sub)
    kp = np.zeros(len(r), dtype=np.int64)
    for M in total_per_bb:
        hit = np.zeros(len(r), dtype=bool)
        for da, db in [(0, x) for x in range(-d, d + 1)] + [(x, 0) for x in range(-d, d + 1) if x]:
            aa, bb = r + da, c + db
            ok = (aa >= 0) & (bb >= 0) & (aa < n) & (bb < n)
            aa, bb = np.clip(aa, 0, n - 1), np.clip(bb, 0, n - 1)
            ok &= (sub[aa] == sub[r]) & (sub[bb] == sub[c])
            v = np.asarray(M[aa, bb]).ravel() > 0
            hit |= ok & v
        kp += hit
    return kp


def install(variant, vd, rep, bbdir):
    from magus.align.merge.graph_build import graph_builder as gb
    orig_transform = vote.transform
    per_bb = []
    if variant.startswith("pool"):
        orig_build = vote.build

        def build(context):  # also keep each backbone's own matrix (for pooling)
            from gcmx.fastgraph import _column_matrix
            from gcmx.weighting import _backbone_alignmap
            from magus.tasks import task
            files = [t.outputFile for t in task.asCompleted(context.backboneTasks)]
            for path in context.backbonePaths:
                if path not in files:
                    files.append(path)
            for f in sorted(files):
                A = _column_matrix(_backbone_alignmap(context, f), context.graph.matrixSize)
                per_bb.append((A.T @ A).tocsr())
            return orig_build(context)
        vote.build = build

    def transform(context, total, support, cover, nbb, var, outdir):
        order = [os.path.abspath(p) for p in context.subalignmentPaths]
        json.dump(order, open(os.path.join(vd, "subsets.json"), "w"))
        if not var.startswith(ORACLE):
            return orig_transform(context, total, support, cover, nbb, var, outdir)
        r, c, w, k, nexp, sub = vote.edge_arrays(context, total, support, cover)
        K, Nn, C, inv = vote.table(k, nexp)
        fit, rc = vote.fit_mixture(K, Nn, C, "bb")
        post = rc[inv]
        info = {"variant": var, "edges": int(len(k))}
        if var.startswith("or-"):
            vmap = bbe_node_map(rep, order)
            keys, uw, uwt, NB = truth(rep, bbdir)
            a, b = vmap[r], vmap[c]
            q = np.minimum(a, b).astype(np.int64) * NB + np.maximum(a, b)
            pos = np.clip(np.searchsorted(keys, q), 0, len(keys) - 1)
            ok = keys[pos] == q
            assert ok.mean() > 0.999, ok.mean()
            istrue = ok & (uwt[pos] / np.maximum(uw[pos], 1) > 0.5)
            if var == "or-allfalse":
                drop = ~istrue
            else:
                _, f, which = var.split("-")
                fdrop = (k < 4) if f == "es4" else (post <= 0.5)
                drop = fdrop & (istrue if which == "true" else ~istrue)
                if which in ("falsenear", "falsefar"):  # split the false deletions by offset to the true partner
                    off = np.full(len(drop), np.inf)
                    sel = np.flatnonzero(drop)
                    off[sel] = false_offsets(rep, a[sel], b[sel])
                    near = off <= 2
                    drop &= near if which == "falsenear" else ~near
                    info["offset_le2_share_of_false_drops"] = float(near[sel].mean())
        else:  # pooling
            name = var
            d = 1
            if "-d" in name:
                name, d = name.split("-d")
                d = int(d)
            kp = pooled(per_bb, r, c, sub, d)
            kp = np.minimum(kp, nexp)  # (pooled votes cannot exceed exposure for the vote model)
            if name == "poolbb":
                K2, N2, C2, inv2 = vote.table(np.maximum(kp, 1), nexp)
                fit2, rc2 = vote.fit_mixture(K2, N2, C2, "bb")
                drop = rc2[inv2] <= 0.5
            else:
                drop = kp < int(name[4:])
            info["mean_kp_minus_k"] = float((kp - np.minimum(k, nexp)).mean())
        info["dropped_edges"] = int(drop.sum())
        info["dropped_weight_frac"] = float(w[drop].sum() / w.sum())
        json.dump(info, open(os.path.join(outdir, "model.json"), "w"), indent=1)
        neww = np.where(drop, 0, w)
        coo = total.tocoo()
        cross = sub[coo.row] != sub[coo.col]
        nz = neww > 0
        rows = np.concatenate([coo.row[~cross], r[nz], c[nz]])
        cols = np.concatenate([coo.col[~cross], c[nz], r[nz]])
        vals = np.concatenate([coo.data[~cross], neww[nz], neww[nz]])
        return sp.csr_matrix((vals, (rows, cols)), shape=total.shape)
    vote.transform = transform
    base = variant if not variant.startswith(ORACLE) else "magus"
    vote.install(base, vd)
    # vote.install's buildMatrix calls vote.transform with ITS variant name; pass ours instead
    from gcmx import fastgraph
    import numpy as _np

    def buildMatrix(context):
        total, support, cover, nbb = vote.build(context)
        total = vote.transform(context, total, support, cover, nbb, variant, vd).tocsr()
        total.indices = total.indices.astype(_np.int32)
        context.graph.matrix = fastgraph.CSRGraph(total)
    gb.buildMatrix = buildMatrix


# ---------------------------------------------------------------- analysis of a finished run

def dissect(rep, vd, out, bbdir):
    """Cross-subset error decomposition at node-pair level + MCL cluster stats (in bbe node space)."""
    import bbe
    from gcmx import fasta
    R = bbe.Rep(rep)
    N = R.nnodes
    order = json.load(open(os.path.join(vd, "subsets.json")))
    vmap = bbe_node_map(rep, order)
    node_sub = np.empty(N, dtype=np.int64)
    o = 0
    for i, s in enumerate(R.subsets):
        L = len(next(iter(s.values())))
        node_sub[o:o + L] = i
        o += L
    # node -> final column
    fc = bbe.final_node_cols(R, out)
    # true pair counts between nodes: T = M M^T with residue counts
    rows, cols, vals = [], [], []
    for node, cnt in R.node_ref.items():
        for rcol, v in cnt.items():
            rows.append(node); cols.append(rcol); vals.append(v)
    M = sp.csr_matrix((np.array(vals, dtype=np.int64), (rows, cols)), shape=(N, max(cols) + 1))
    T = sp.triu(M @ M.T, k=1).tocoo()
    keep = node_sub[T.row] != node_sub[T.col]
    ta, tb, tv = T.row[keep], T.col[keep], T.data[keep]
    together = (fc[ta] == fc[tb]) & (fc[ta] >= 0)
    # graph state between the two nodes: k in the full graph, and kept in this variant's graph?
    keys, uw, uwt, _ = truth(rep, bbdir)
    q = ta.astype(np.int64) * N + tb
    pos = np.clip(np.searchsorted(keys, q), 0, len(keys) - 1)
    inG = keys[pos] == q
    kept = np.zeros(len(q), dtype=bool)
    gfile = os.path.join(vd, "kept_keys.npy")
    if os.path.exists(gfile):
        kk = np.load(gfile)
        p2 = np.clip(np.searchsorted(kk, q), 0, len(kk) - 1)
        kept = kk[p2] == q
    tot = tv.sum()
    miss = ~together
    res = {"x_true_pairs": int(tot), "xSPFN": float(tv[miss].sum() / tot),
           "miss_noedge": float(tv[miss & ~inG].sum() / tot),
           "miss_edge_dropped": float(tv[miss & inG & ~kept].sum() / tot),
           "miss_edge_kept": float(tv[miss & inG & kept].sum() / tot)}
    # false pairs placed together, by graph state
    sizes = R.node_size
    colnodes = {}
    for nd in np.flatnonzero(fc >= 0):
        colnodes.setdefault(int(fc[nd]), []).append(nd)
    fa, fb = [], []
    for nds in colnodes.values():
        if len(nds) < 2:
            continue
        nds = np.array(nds)
        i, j = np.triu_indices(len(nds), 1)
        a, b = nds[i], nds[j]
        s = node_sub[a] != node_sub[b]
        fa.append(np.minimum(a[s], b[s])); fb.append(np.maximum(a[s], b[s]))
    fa, fb = np.concatenate(fa), np.concatenate(fb)
    allp = sizes[fa] * sizes[fb]
    Tcsr = sp.csr_matrix((tv, (ta, tb)), shape=(N, N))
    tp = np.asarray(Tcsr[fa, fb]).ravel()
    fp = allp - tp
    q2 = fa.astype(np.int64) * N + fb
    p3 = np.clip(np.searchsorted(keys, q2), 0, len(keys) - 1)
    inG2 = keys[p3] == q2
    res.update(x_est_pairs=int(allp.sum()), xSPFP=float(fp.sum() / allp.sum()),
               fp_noedge=float(fp[~inG2].sum() / allp.sum()), fp_edge=float(fp[inG2].sum() / allp.sum()))
    # MCL clusters (pre-trace) and trace clusters
    for nm in ("clusters", "trace"):
        f = os.path.join(vd, nm + ".txt")
        if not os.path.exists(f):
            continue
        cl = [vmap[np.array([int(x) for x in line.split()], dtype=np.int64)]
              for line in open(f) if line.strip()]
        cl = [c for c in cl if len(c)]
        multi = [c for c in cl if len(c) > 1]
        viol = [c for c in multi if len(set(node_sub[c].tolist())) < len(c)]
        covered = np.zeros(N, dtype=bool)
        for c in cl:
            covered[c] = True
        # over-splitting: true pairs whose nodes are in different (clustered) clusters
        lab = np.full(N, -1, dtype=np.int64)
        for i, c in enumerate(cl):
            lab[c] = i
        both = (lab[ta] >= 0) & (lab[tb] >= 0)
        res[nm] = {"clusters": len(multi), "mean_size": float(np.mean([len(c) for c in multi])) if multi else 0,
                   "violating_frac": len(viol) / max(len(multi), 1),
                   "nodes_in_violating_frac": float(sum(len(c) for c in viol) / max(sum(len(c) for c in multi), 1)),
                   "nodes_unclustered_frac": float(1 - covered.mean()),
                   "true_pairs_split_between_clusters": float(tv[both & (lab[ta] != lab[tb])].sum() / tot),
                   "true_pairs_unclustered": float(tv[~both].sum() / tot)}
    return res


def flags(v):
    """MAGUS merge flags; a variant suffix @fX sets the MCL inflation factor (H5)."""
    fl = list(vote.MERGE_FLAGS)
    if "@f" in v:
        fl[fl.index("-f") + 1] = v.split("@f")[1]
    return fl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rep")
    ap.add_argument("variants", nargs="+")
    ap.add_argument("--bb")
    ap.add_argument("--name")
    a = ap.parse_args()
    rep = os.path.abspath(a.rep)
    bbdir = os.path.abspath(a.bb or os.path.join(rep, "inputs", "backbones"))
    res = os.path.join(rep, "why_results.jsonl")
    done = {json.loads(l)["variant"] for l in open(res)} if os.path.exists(res) else set()
    for v in a.variants:
        if v in done:
            continue
        vd = os.path.join(rep, "why", v)
        shutil.rmtree(vd, ignore_errors=True)
        os.makedirs(vd)
        out = os.path.join(vd, "out.fasta")
        t0 = time.time()
        with open(os.path.join(vd, "magus.log"), "w") as log:
            r = subprocess.run([sys.executable, os.path.abspath(__file__), "--_inner", v, vd, rep, bbdir,
                                "-np", "4", "-d", os.path.join(vd, "work"), "-s", os.path.join(rep, "inputs", "subalignments"),
                                "-b", bbdir, "-o", out] + flags(v), cwd=os.path.join(ROOT, "code"),
                               stdout=log, stderr=subprocess.STDOUT)
        if r.returncode:
            print("FAILED", rep, v, flush=True)
            continue
        wall = round(time.time() - t0, 1)
        for nm in ("clusters", "trace"):
            shutil.copy(os.path.join(vd, "work", "graph", nm + ".txt"), os.path.join(vd, nm + ".txt"))
        # kept cross-subset edges of this variant's graph (bbe keys), from graph.txt
        order = json.load(open(os.path.join(vd, "subsets.json")))
        vmap = bbe_node_map(rep, order)
        import pandas as pd
        g = pd.read_csv(os.path.join(vd, "work", "graph", "graph.txt"), sep=r"\s+", header=None, usecols=[0, 1],
                        dtype=np.int64).values
        ga, gb_ = vmap[g[:, 0]], vmap[g[:, 1]]
        _, _, _, NB = truth(rep, bbdir)
        np.save(os.path.join(vd, "kept_keys.npy"),
                np.unique(np.minimum(ga, gb_).astype(np.int64) * NB + np.maximum(ga, gb_)))
        shutil.rmtree(os.path.join(vd, "work"), ignore_errors=True)
        from gcmx.bbtool_bench import acc_ref
        s = acc_ref(os.path.join(rep, "true.fasta"), out)
        row = {"rep": a.name or os.path.basename(rep), "variant": v, "merge_wall": wall, **s}
        row.update(dissect(rep, vd, out, bbdir))
        m = os.path.join(vd, "model.json")
        if os.path.exists(m):
            mj = json.load(open(m))
            row["model"] = {kk: mj[kk] for kk in ("dropped_edges", "dropped_weight_frac", "kept_edges", "kept_weight_frac",
                                                   "mean_kp_minus_k", "offset_le2_share_of_false_drops") if kk in mj}
        with open(res, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps({kk: row[kk] for kk in ("rep", "variant", "SPFN", "SPFP") if kk in row}), flush=True)


def inner():
    variant, vd, rep, bbdir = sys.argv[2:6]
    variant = variant.split("@f")[0]
    sys.argv = [sys.argv[0]] + sys.argv[6:]
    sys.path.insert(0, os.path.join(ROOT, "code"))
    install(variant, vd, rep, bbdir)
    from magus.main import main as magus_main
    magus_main()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--_inner":
        inner()
    else:
        main()
