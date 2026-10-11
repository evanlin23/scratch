# AI-assisted (Claude), exploration code for CS581 project
"""Why does deleting low-support GCM edges help proteins but not DNA/RNA?  Edge-level diagnostics.

    python3 why.py REPDIR [--bb DIR] [--name NAME] [--type protein|dna]   -> results/edges/<NAME>.json

REPDIR in the rep-bank layout (inputs/subalignments, inputs/backbones, true.fasta). For every
cross-subset edge (a, b) of MAGUS's graph built from the B backbones:
  w  = backbone residue pairs (MAGUS's weight), wt = of those, true homologies in the reference,
  k  = backbones contributing >= 1 pair, n = backbones covering both nodes (exposure).
An edge is "true" if wt / w > 0.5 (the gcmgen/calib edge-level definition); units = residue pairs.

H1 headroom : precision by k (edges and units); share of weight / false weight / true weight in k <= 3;
              what es4 and hard-bb (beta-binomial vote, refit here) remove: false units, true units.
H2 near miss: for false edges, distance (in subset-B columns) from b to the nearest column of B that
              truly pairs with a, and symmetrically a; offset = min of both. Also the pooled support
              kp(a, b) = backbones that align a with b-1, b or b+1 (or a-1, a, a+1 with b).
H3 features : mean p-distance, gap fraction (reference and subsets), backbone set Jaccard, low-k share.
"""
import argparse
import itertools
import json
import os
import random
import sys

import numpy as np
import scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "bbevidence", "code"))
sys.path.insert(0, os.path.join(ROOT, "gcmvote", "code"))
import bbe  # noqa: E402
from bbe import fasta, residue_columns, evidence_mask  # noqa: E402
import vote  # noqa: E402

KBINS = [(1, 1), (2, 2), (3, 3), (4, 5), (6, 8), (9, 99)]
OBINS = [(0, 0), (1, 1), (2, 2), (3, 5), (6, 10**9)]


def per_backbone(R, files):
    """Per backbone: unique cross-subset edge keys, their units and true units; node coverage."""
    N = R.nnodes
    out = []
    cover = np.zeros(N, dtype=np.int64)
    for bi, (_, a) in enumerate(files):
        rc = residue_columns(a)
        cols, nodes, refs, subs = [], [], [], []
        for t in a:
            m = evidence_mask(a[t])
            cols.append(rc[t][m]); nodes.append(R.subcol[t][m]); refs.append(R.refcol[t][m])
            subs.append(np.full(int(m.sum()), R.subset_of[t], dtype=np.int32))
        cols, nodes, refs, subs = map(np.concatenate, (cols, nodes, refs, subs))
        cover[nodes] |= (1 << bi)
        order = np.argsort(cols, kind="stable")
        cols, nodes, refs, subs = cols[order], nodes[order], refs[order], subs[order]
        bounds = np.flatnonzero(np.diff(cols)) + 1
        ks, ts = [], []
        for lo, hi in zip(np.r_[0, bounds], np.r_[bounds, len(cols)]):
            if hi - lo < 2:
                continue
            i, j = np.triu_indices(hi - lo, 1)
            i, j = i + lo, j + lo
            keep = subs[i] != subs[j]
            i, j = i[keep], j[keep]
            na, nb = np.minimum(nodes[i], nodes[j]), np.maximum(nodes[i], nodes[j])
            ks.append(na.astype(np.int64) * N + nb)
            ts.append(refs[i] == refs[j])
        k = np.concatenate(ks); t = np.concatenate(ts)
        uk, inv = np.unique(k, return_inverse=True)
        out.append((uk, np.bincount(inv), np.bincount(inv, weights=t)))
    return out, cover


def popcount(x):
    x = x.astype(np.uint64).copy()
    c = np.zeros(len(x), dtype=np.int64)
    while x.any():
        c += (x & np.uint64(1)).astype(np.int64)
        x >>= np.uint64(1)
    return c


def true_pair_keys(R):
    """Sorted keys a * N + b (both directions) of node pairs from different subsets that share >= 1
    reference column (i.e. hold at least one true homologous residue pair)."""
    N = R.nnodes
    rows, cols = [], []
    for node, cnt in R.node_ref.items():
        for rc in cnt:
            rows.append(node); cols.append(rc)
    ncol = max(cols) + 1
    M = sp.csr_matrix((np.ones(len(rows), dtype=np.int32), (rows, cols)), shape=(N, ncol))
    T = (M @ M.T).tocoo()
    sub = R.node_sub
    keep = sub[T.row] != sub[T.col]
    keys = T.row[keep].astype(np.int64) * N + T.col[keep]
    return np.unique(keys)


def nearest_offset(R, tkeys, a, b):
    """Distance from b to the nearest node of b's subset that truly pairs with a (inf if none)."""
    N = R.nnodes
    q = a.astype(np.int64) * N + b
    pos = np.searchsorted(tkeys, q)
    best = np.full(len(q), np.inf)
    lo_sub, hi_sub = R.sub_lo[R.node_sub[b]], R.sub_hi[R.node_sub[b]]
    for p in (pos - 1, pos):
        ok = (p >= 0) & (p < len(tkeys))
        kk = tkeys[np.clip(p, 0, len(tkeys) - 1)]
        ra, cb = kk // N, kk % N
        ok &= (ra == a) & (cb >= lo_sub) & (cb < hi_sub)
        best = np.where(ok, np.minimum(best, np.abs(cb - b)), best)
    return best


def pooled_support(R, bbs, a, b, d=1):
    """Backbones that align a with some b' in [b-d, b+d] or some a' in [a-d, a+d] with b (same subset)."""
    N = R.nnodes
    kp = np.zeros(len(a), dtype=np.int64)
    sa, sb = R.node_sub[a], R.node_sub[b]
    for uk, _, _ in bbs:
        hit = np.zeros(len(a), dtype=bool)
        for da, db in [(0, x) for x in range(-d, d + 1)] + [(x, 0) for x in range(-d, d + 1) if x]:
            aa, bb = a + da, b + db
            ok = (aa >= 0) & (bb >= 0) & (aa < N) & (bb < N)
            aa, bb = np.clip(aa, 0, N - 1), np.clip(bb, 0, N - 1)
            ok &= (R.node_sub[aa] == sa) & (R.node_sub[bb] == sb)
            lo, hi = np.minimum(aa, bb), np.maximum(aa, bb)
            q = lo.astype(np.int64) * N + hi
            p = np.clip(np.searchsorted(uk, q), 0, len(uk) - 1)
            hit |= ok & (uk[p] == q)
        kp += hit
    return kp


def pdist_stats(ref, nsamp=400, seed=1):
    names = list(ref)
    rng = random.Random(seed)
    arr = np.frombuffer("".join(ref[t].upper() for t in names).encode(), dtype=np.uint8).reshape(len(names), -1)
    gap = (arr == ord("-")) | (arr == ord("."))
    pd = []
    for _ in range(nsamp):
        i, j = rng.sample(range(len(names)), 2)
        both = ~gap[i] & ~gap[j]
        if both.sum():
            pd.append(float((arr[i][both] != arr[j][both]).mean()))
    return float(np.mean(pd)), float(gap.mean())


def summarize(name, dtype, R, bbs, cover, B, ref, rep, bbdir):
    N = R.nnodes
    allk = np.concatenate([uk for uk, _, _ in bbs])
    keys, inv = np.unique(allk, return_inverse=True)
    w = np.bincount(inv, weights=np.concatenate([x for _, x, _ in bbs]))
    wt = np.bincount(inv, weights=np.concatenate([x for _, _, x in bbs]))
    k = np.bincount(inv)
    a, b = keys // N, keys % N
    n = popcount(cover[a] & cover[b])
    true = wt / w > 0.5
    wf = w - wt
    res = {"rep": name, "type": dtype, "B": B, "nodes": int(N), "edges": int(len(keys)),
           "units": float(w.sum()), "unit_prec": float(wt.sum() / w.sum()), "edge_prec": float(true.mean())}
    # ---- H1: precision by k
    rows = []
    for lo, hi in KBINS:
        s = (k >= lo) & (k <= hi)
        if not s.any():
            continue
        rows.append({"k": "{}-{}".format(lo, hi) if lo != hi else str(lo), "edges": int(s.sum()),
                     "edge_share": float(s.mean()), "unit_share": float(w[s].sum() / w.sum()),
                     "edge_prec": float(true[s].mean()), "unit_prec": float(wt[s].sum() / w[s].sum()),
                     "units_per_edge": float(w[s].mean())})
    res["by_k"] = rows
    low = k <= 3
    res["low_share_edges"] = float(low.mean())
    res["low_share_units"] = float(w[low].sum() / w.sum())
    res["false_units_share_total"] = float(wf.sum() / w.sum())
    res["low_false_of_false_units"] = float(wf[low].sum() / wf.sum())
    res["low_true_of_true_units"] = float(wt[low].sum() / wt.sum())
    res["false_edges_low_share"] = float(low[~true].mean())
    # ---- what the filters remove (graph level): es4 and hard-bb refit
    K, Nn, C, cinv = vote.table(k, n)
    fit, rcell = vote.fit_mixture(K, Nn, C, "bb")
    post = rcell[cinv]
    res["bb_fit"] = fit
    for fname, drop in (("es4", k < 4), ("hard-bb", post <= 0.5)):
        res["cut_" + fname] = {
            "edges_removed": float(drop.mean()), "units_removed": float(w[drop].sum() / w.sum()),
            "false_units_removed_of_false": float(wf[drop].sum() / wf.sum()),
            "true_units_removed_of_true": float(wt[drop].sum() / wt.sum()),
            "removed_prec": float(wt[drop].sum() / max(w[drop].sum(), 1)),
            "true_edges_removed_of_true": float((drop & true).sum() / true.sum()),
            "false_edges_removed_of_false": float((drop & ~true).sum() / (~true).sum()),
            # precision of the graph before/after
            "prec_after_units": float(wt[~drop].sum() / w[~drop].sum())}
    # ---- H6 (added): which TRUE edges the filters delete, by node size, and orphaned nodes
    size = R.node_size
    ms = np.minimum(size[a], size[b])
    res["h6_sizebins"] = []
    for lo, hi in ((1, 1), (2, 3), (4, 10), (11, 10**9)):
        s_ = true & (ms >= lo) & (ms <= hi)
        row = {"min_node_size": "{}-{}".format(lo, hi if hi < 10**9 else ""), "true_edges": int(s_.sum()),
               "true_units_share": float(wt[s_].sum() / wt.sum())}
        for fname, drop in (("es4", k < 4), ("hard-bb", post <= 0.5)):
            row[fname + "_true_edges_removed"] = float(drop[s_].mean()) if s_.any() else None
            row[fname + "_true_units_removed"] = float(wt[s_ & drop].sum() / max(wt[s_].sum(), 1))
        row["mean_k"] = float(k[s_].mean()) if s_.any() else None
        row["mean_n"] = float(n[s_].mean()) if s_.any() else None
        res["h6_sizebins"].append(row)
    res["h6_residues_in_small_nodes"] = float(size[size <= 3].sum() / size.sum())
    res["h6_nodes_small"] = float((size <= 3).mean())
    for fname, drop in (("es4", k < 4), ("hard-bb", post <= 0.5)):
        has = np.zeros(N, dtype=bool); has[a[true]] = True; has[b[true]] = True
        kept = np.zeros(N, dtype=bool); kt = true & ~drop
        kept[a[kt]] = True; kept[b[kt]] = True
        orphan = has & ~kept
        # true partner subsets lost: per node, number of subsets it has a true edge to, before/after
        res["h6_orphan_" + fname] = {"nodes": float(orphan.sum() / has.sum()),
                                     "residues": float(size[orphan].sum() / size[has].sum())}
    # ---- H2: offsets of false edges to the nearest true partner
    R.node_sub = np.empty(N, dtype=np.int64)
    R.sub_lo, R.sub_hi = [], []
    o = 0
    for i, s in enumerate(R.subsets):
        L = len(next(iter(s.values())))
        R.node_sub[o:o + L] = i
        R.sub_lo.append(o); R.sub_hi.append(o + L)
        o += L
    R.sub_lo, R.sub_hi = np.array(R.sub_lo), np.array(R.sub_hi)
    tkeys = true_pair_keys(R)
    fa = ~true
    off = np.minimum(nearest_offset(R, tkeys, a[fa], b[fa]), nearest_offset(R, tkeys, b[fa], a[fa]))
    kf, wff = k[fa], wf[fa]
    tab = []
    for lo, hi in KBINS + [(1, 3), (4, 99)]:
        s = (kf >= lo) & (kf <= hi)
        if not s.any():
            continue
        row = {"k": "{}-{}".format(lo, hi) if lo != hi else str(lo), "false_edges": int(s.sum()),
               "false_units": float(wff[s].sum())}
        for olo, ohi in OBINS:
            t = s & (off >= olo) & (off <= ohi)
            row["off{}".format(olo if olo == ohi else "{}-{}".format(olo, ohi if ohi < 10**9 else ""))] = \
                float(wff[t].sum() / wff[s].sum())
        row["off_none"] = float(wff[s & ~np.isfinite(off)].sum() / wff[s].sum())
        row["off_le2_edges"] = float(((off[s] <= 2)).mean())
        tab.append(row)
    res["h2_offsets_units"] = tab
    # pooled support (near-miss pooling, d = 1)
    kp = pooled_support(R, bbs, a, b, 1)
    res["h2_pool"] = {}
    for thr in (4,):
        for rule, keep in (("k>=4", k >= thr), ("kpool>=4", kp >= thr),
                           ("hard-bb", post > 0.5)):
            res["h2_pool"][rule] = {"true_units_kept": float(wt[keep].sum() / wt.sum()),
                                    "false_units_kept": float(wf[keep].sum() / wf.sum()),
                                    "prec_kept": float(wt[keep].sum() / w[keep].sum())}
    # vote splitting: edges hard-bb deletes, true vs false: own support k vs pooled support kp (d = 1)
    dh = post <= 0.5
    res["h2_split"] = {}
    for lab, sel in (("true_deleted", dh & true), ("false_deleted", dh & ~true), ("true_kept", ~dh & true),
                     ("false_kept", ~dh & ~true)):
        res["h2_split"][lab] = {"edges": int(sel.sum()), "mean_k": float(k[sel].mean()) if sel.any() else None,
                                "mean_kp": float(kp[sel].mean()) if sel.any() else None,
                                "mean_n": float(n[sel].mean()) if sel.any() else None,
                                "share_kp_ge_n_minus_1": float((kp[sel] >= n[sel] - 1).mean()) if sel.any() else None}
    # rescued by pooling: low-k true edges
    lt = low & true
    res["h2_pool"]["lowk_true_edges"] = int(lt.sum())
    res["h2_pool"]["lowk_true_rescued_frac_edges"] = float((kp[lt] >= 4).mean()) if lt.any() else None
    lf = low & ~true
    res["h2_pool"]["lowk_false_rescued_frac_edges"] = float((kp[lf] >= 4).mean()) if lf.any() else None
    res["h2_pool"]["mean_kp_minus_k_true_low"] = float((kp[lt] - k[lt]).mean()) if lt.any() else None
    res["h2_pool"]["mean_kp_minus_k_false_low"] = float((kp[lf] - k[lf]).mean()) if lf.any() else None
    # ---- H3 features
    pd, rg = pdist_stats(ref)
    sg = np.mean([np.mean([s.count("-") / len(s) for s in sa.values()]) for sa in R.subsets])
    res["feat"] = {"pdist": pd, "ref_gapfrac": rg, "sub_gapfrac": float(sg), "ntaxa": len(ref)}
    setdir = os.path.join(rep, "sets")
    bbsets = []
    for f in sorted(os.listdir(bbdir)):
        bbsets.append(set(fasta.read(os.path.join(bbdir, f))))
    jac = [len(x & y) / len(x | y) for x, y in itertools.combinations(bbsets, 2)]
    res["feat"]["bb_jaccard"] = float(np.mean(jac))
    res["feat"]["bb_size"] = float(np.mean([len(x) for x in bbsets]))
    res["feat"]["mean_exposure"] = float(np.average(n, weights=w))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rep")
    ap.add_argument("--bb")
    ap.add_argument("--name")
    ap.add_argument("--type", required=True)
    ap.add_argument("--out", default=os.path.join(ROOT, "gcmwhy", "results", "edges"))
    x = ap.parse_args()
    rep = os.path.abspath(x.rep)
    name = x.name or os.path.basename(rep)
    bbdir = os.path.abspath(x.bb or os.path.join(rep, "inputs", "backbones"))
    R = bbe.Rep(rep)
    files = [(f, fasta.upper(fasta.read(os.path.join(bbdir, f)))) for f in sorted(os.listdir(bbdir))]
    bbs, cover = per_backbone(R, files)
    ref = fasta.read(os.path.join(rep, "true.fasta"))
    res = summarize(name, x.type, R, bbs, cover, len(files), ref, rep, bbdir)
    mj = os.path.join(rep, "magus.json")
    if os.path.exists(mj):
        res["magus_json"] = json.load(open(mj)).get("magus")
    os.makedirs(x.out, exist_ok=True)
    json.dump(res, open(os.path.join(x.out, name + ".json"), "w"), indent=1, default=float)
    print(json.dumps({kk: res[kk] for kk in ("rep", "edges", "unit_prec", "low_share_units",
                                               "low_false_of_false_units", "low_true_of_true_units")}))


if __name__ == "__main__":
    main()
