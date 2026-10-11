# AI-assisted (Claude), exploration code for CS581 project
"""Calibration of the vote posterior against the reference (diagnostic only).

    python3 calib.py REP [TAG]    -> REP/calib.json  (TAG = vote/<TAG>/edges.npz, default magus)

Truth per edge = gcmgen's definition: of the edge's backbone residue-pair units, the fraction that are true
homologies in the reference (bbe.graph_edges). Reports, per posterior bin (binomial and beta-binomial):
edges, mean posterior, unit precision, edge-level precision (share of edges with > 50 % true units);
and per support k: precision (as in gcmgen) for comparison.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gcmgen", "code"))
import gg  # noqa: E402
bbe = gg.bbe


def main(rep, tag="magus"):
    vd = os.path.join(rep, "vote", tag)
    E = np.load(os.path.join(vd, "edges.npz"))
    order = json.load(open(os.path.join(vd, "subsets.json")))
    R = bbe.Rep(rep)
    files = [f for f, _ in bbe.subsets(rep)]
    # MAGUS's node ids (vote) -> bbe node ids
    lens = [len(next(iter(a.values()))) for a in R.subsets]
    off_bbe = np.concatenate([[0], np.cumsum(lens)])[:-1]
    idx = {f: i for i, f in enumerate(files)}
    vmap = np.concatenate([off_bbe[idx[os.path.basename(p)]] + np.arange(lens[idx[os.path.basename(p)]])
                           for p in order])
    bbfiles = sorted(os.listdir(os.path.join(rep, "inputs", "backbones")))
    if "@B" in tag:
        B = int(tag.split("@B")[1])
        d = os.path.join(rep, "bb{}".format(B))
        bbfiles = [os.path.join(d, f) for f in sorted(os.listdir(d))]
    else:
        bbfiles = [os.path.join(rep, "inputs", "backbones", f) for f in bbfiles]
    from gcmx import fasta
    keys, w, wt, nbb = bbe.graph_edges(R, [(os.path.basename(f), fasta.upper(fasta.read(f))) for f in bbfiles])
    a, b = vmap[E["a"]], vmap[E["b"]]
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    vk = lo.astype(np.int64) * R.nnodes + hi
    pos = np.searchsorted(keys, vk)
    ok = (pos < len(keys)) & (keys[np.minimum(pos, len(keys) - 1)] == vk)
    assert ok.mean() > 0.999, ok.mean()
    pos = pos[ok]
    assert (nbb[pos] == E["k"][ok]).mean() > 0.999
    uw, ut = w[pos].astype(float), wt[pos].astype(float)
    frac = ut / np.maximum(uw, 1)
    out = {"rep": os.path.basename(rep), "tag": tag, "edges": int(ok.sum())}
    bins = np.linspace(0, 1, 11)
    for kind in ("post_binom", "post_bb"):
        p = E[kind][ok]
        rows = []
        for lo_, hi_ in zip(bins[:-1], bins[1:]):
            s = (p >= lo_) & ((p < hi_) if hi_ < 1 else (p <= hi_))
            if s.sum() == 0:
                continue
            rows.append({"bin": [round(lo_, 1), round(hi_, 1)], "edges": int(s.sum()), "mean_post": float(p[s].mean()),
                         "unit_prec": float(ut[s].sum() / uw[s].sum()), "edge_prec": float((frac[s] > 0.5).mean())})
        out[kind] = rows
        # expected calibration error (edge-level), Brier
        out[kind + "_ece"] = float(sum(r["edges"] * abs(r["mean_post"] - r["edge_prec"]) for r in rows) / ok.sum())
        out[kind + "_kept_true_frac"] = float(ut[p > 0.5].sum() / ut.sum())  # recall of true units kept by hard
        out[kind + "_kept_prec"] = float(ut[p > 0.5].sum() / max(uw[p > 0.5].sum(), 1))
    k = E["k"][ok]
    out["by_k"] = [{"k": int(kk), "edges": int((k == kk).sum()), "unit_prec": float(ut[k == kk].sum() / uw[k == kk].sum())}
                   for kk in np.unique(k)]
    n = E["n"][ok]
    out["by_kn"] = [{"k": int(kk), "n": int(nn), "edges": int(s.sum()), "edge_prec": float((frac[s] > 0.5).mean()),
                     "post_binom": float(E["post_binom"][ok][s].mean()), "post_bb": float(E["post_bb"][ok][s].mean())}
                    for nn in np.unique(n) for kk in range(1, nn + 1) for s in [(k == kk) & (n == nn)] if s.sum() > 0]
    out["overall_edge_prec"] = float((frac > 0.5).mean())
    json.dump(out, open(os.path.join(rep, "calib{}.json".format("" if tag == "magus" else "_" + tag)), "w"), indent=1)
    return out


if __name__ == "__main__":
    o = main(os.path.abspath(sys.argv[1]), *(sys.argv[2:3]))
    for kind in ("post_binom", "post_bb"):
        print(kind, "ECE", round(o[kind + "_ece"], 3))
        for r in o[kind]:
            print("  {bin} n={edges} post={mean_post:.2f} unit={unit_prec:.2f} edge={edge_prec:.2f}".format(**r))
