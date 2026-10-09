#!/usr/bin/env python3
"""Aggregate result.json files into the report tables.

Replicates 01-05 are the development set (looked at while building methods);
06-20 are held out and are the only ones used for the paired tests.
Error = normalized RF = (FN+FP) / (2 * #internal edges of the true tree).
Paired test: two-sided Wilcoxon signed-rank on per-replicate error differences (ties dropped).
W/T/L: candidate error lower / equal (tie band: exactly equal, errors are discrete) / higher.
"""
import glob, json, os, sys
from collections import defaultdict
import numpy as np
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
DEV = {f"{i:02d}" for i in range(1, 6)}

RATE = {"0": "0", "0.000000002": "0.08", "0.000000005": "0.2", "0.00000002": "0.8", "0.0000002": "8", "0.0000005": "20"}
PAIRS = [("ls[astral4]", "astral4"), ("ls[treeqmc]", "treeqmc"), ("ls[wqfm]", "wqfm"), ("ls[astrid]", "astrid"),
         ("harvest", "astral4"), ("astral3+e", "astral3"), ("astral3", "astral4"), ("astral3+e", "astral4"), ("capminor", "astral4"), ("vote", "astral4"), ("reweight", "astral4"),
         ("treeqmc", "astral4"), ("wqfm", "astral4"), ("astrid", "astral4"), ("wastral", "astral4")]
METHODS = ["astral4", "astral3", "astral3+e", "wastral", "treeqmc", "wqfm", "astrid", "ls[astral4]", "ls[treeqmc]", "ls[wqfm]", "ls[astrid]",
           "harvest", "capminor", "vote", "reweight"]


def load():
    rows = []  # (dataset, condition, ngenes, rep, results)
    for f in glob.glob("/opt/runs/hgt/*/*/*/result.json"):
        parts = f.split("/")
        model, rep, ng = parts[-4], parts[-3], parts[-2]
        cond = "HGT " + RATE[model.split("0.000001.", 1)[1]] + "/gene"
        rows.append(("hgt", cond, ng, rep, json.load(open(f))))
    for f in glob.glob("/opt/runs/miss/*/*/*/result.json"):
        parts = f.split("/")
        c, rep, ng = parts[-4], parts[-3], parts[-2]
        rows.append(("ils", "ILS " + c, ng + "-raxml", rep, json.load(open(f))))
    return rows


def err(r):
    return (r["fn"] + r["fp"]) / (2 * r["nint"])


def fmt_p(p):
    return "n/a" if p is None or np.isnan(p) else (f"{p:.2g}" if p >= 1e-3 else f"{p:.0e}")


def paired(groups, cand, base):
    d, ds = [], []
    for res in groups:
        if cand in res and base in res:
            d.append(err(res[cand]) - err(res[base]))
            ds.append(res[cand]["nscore"] - res[base]["nscore"])
    if not d:
        return None
    d = np.array(d)
    w, t, l = int((d < 0).sum()), int((d == 0).sum()), int((d > 0).sum())
    p = wilcoxon(d[d != 0]).pvalue if (d != 0).sum() >= 1 else float("nan")
    ds = np.array(ds)
    return dict(n=len(d), mean=d.mean(), w=w, t=t, l=l, p=p, dscore=ds.mean(),
                sup=int((ds > 1e-12).sum()), sdown=int((ds < -1e-12).sum()))


def main():
    rows = load()
    out = []
    by = defaultdict(list)
    for ds, cond, ng, rep, res in rows:
        by[(ds, ng, cond)].append((rep, res))
    keys = sorted(by, key=lambda k: (k[0], k[1], k[2]))
    out.append("## Mean error (normalized RF, %) per method; held-out replicates 06-20\n")
    hdr = "| data | genes | condition | n | " + " | ".join(METHODS) + " |"
    out.append(hdr)
    out.append("|" + "---|" * (4 + len(METHODS)))
    for k in keys:
        test = [r for rep, r in by[k] if rep not in DEV]
        if not test:
            continue
        cells = []
        for m in METHODS:
            v = [err(r[m]) for r in test if m in r]
            cells.append(f"{100*np.mean(v):.1f}" if v else "-")
        out.append(f"| {k[0]} | {k[1]} | {k[2]} | {len(test)} | " + " | ".join(cells) + " |")

    out.append("\n## Mean normalized quartet score (x100) per method; held-out replicates\n")
    out.append(hdr)
    out.append("|" + "---|" * (4 + len(METHODS)))
    for k in keys:
        test = [r for rep, r in by[k] if rep not in DEV]
        if not test:
            continue
        cells = []
        for m in METHODS:
            v = [r[m]["nscore"] for r in test if m in r]
            cells.append(f"{100*np.mean(v):.3f}" if v else "-")
        out.append(f"| {k[0]} | {k[1]} | {k[2]} | {len(test)} | " + " | ".join(cells) + " |")

    out.append("\n## Paired comparisons, pooled over conditions within each dataset/gene count (held-out reps)\n")
    out.append("Δerr = mean(candidate − baseline) normalized RF in percentage points (negative = candidate better). "
               "W/T/L = candidate better / tie (identical RF) / worse. p = two-sided Wilcoxon signed-rank (ties dropped). "
               "Δscore = mean normalized quartet score difference (x100); ↑/↓ = #reps where candidate score is higher/lower.\n")
    out.append("| data | genes | candidate | baseline | n | Δerr (pp) | W/T/L | p | Δscore x100 | ↑/↓ |")
    out.append("|---|---|---|---|---|---|---|---|---|---|")
    pooled = defaultdict(list)
    for (ds, ng, cond), lst in by.items():
        pooled[(ds, ng)] += [r for rep, r in lst if rep not in DEV]
    for (ds, ng) in sorted(pooled):
        for c, b in PAIRS:
            s = paired(pooled[(ds, ng)], c, b)
            if s:
                out.append(f"| {ds} | {ng} | {c} | {b} | {s['n']} | {100*s['mean']:+.2f} | {s['w']}/{s['t']}/{s['l']} | "
                           f"{fmt_p(s['p'])} | {100*s['dscore']:+.4f} | {s['sup']}/{s['sdown']} |")

    out.append("\n## Per-condition paired comparisons for the HGT-aware variants (held-out reps)\n")
    out.append("| data | genes | condition | candidate | n | Δerr (pp) | W/T/L | p |")
    out.append("|---|---|---|---|---|---|---|---|")
    for k in keys:
        test = [r for rep, r in by[k] if rep not in DEV]
        for c in ["capminor", "vote", "reweight", "harvest", "ls[treeqmc]"]:
            s = paired(test, c, "astral4" if c != "ls[treeqmc]" else "treeqmc")
            if s:
                out.append(f"| {k[0]} | {k[1]} | {k[2]} | {c} | {s['n']} | {100*s['mean']:+.2f} | "
                           f"{s['w']}/{s['t']}/{s['l']} | {fmt_p(s['p'])} |")

    out.append("\n## Headroom: does ASTRAL-IV already reach the best quartet score found by any method? (held-out reps)\n")
    out.append("best = max normalized quartet score over all methods and local searches on that replicate. "
               "'err(best)' = mean error of the best-scoring tree (ties broken by the ASTRAL-IV tree).\n")
    out.append("| data | genes | n | ASTRAL-IV = best | mean gap x100 | max gap x100 | err ASTRAL-IV % | err best-score tree % | err true-tree-closest method % |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    for (ds, ng) in sorted(pooled):
        L = [r for r in pooled[(ds, ng)] if "astral4" in r]
        if not L:
            continue
        gaps, eb, ea, eo = [], [], [], []
        for r in L:
            ms = [m for m in METHODS if m in r]
            best = max(ms, key=lambda m: (r[m]["nscore"], m == "astral4"))
            gaps.append(r[best]["nscore"] - r["astral4"]["nscore"])
            eb.append(err(r[best])); ea.append(err(r["astral4"])); eo.append(min(err(r[m]) for m in ms))
        g = np.array(gaps)
        out.append(f"| {ds} | {ng} | {len(L)} | {int((g <= 1e-12).sum())}/{len(L)} | {100*g.mean():.4f} | {100*g.max():.4f} | "
                   f"{100*np.mean(ea):.2f} | {100*np.mean(eb):.2f} | {100*np.mean(eo):.2f} |")

    out.append("\n## Is the quartet score aligned with accuracy? Score of the TRUE species tree vs ASTRAL-IV (held-out reps)\n")
    out.append("| data | genes | n | true < ASTRAL-IV | true = | true > | mean (true − ASTRAL-IV) x100 |")
    out.append("|---|---|---|---|---|---|---|")
    for (ds, ng) in sorted(pooled):
        L = [r for r in pooled[(ds, ng)] if "astral4" in r and "truetree" in r]
        if not L:
            continue
        d = np.array([r["truetree"]["nscore"] - r["astral4"]["nscore"] for r in L])
        out.append(f"| {ds} | {ng} | {len(L)} | {int((d < -1e-12).sum())} | {int((abs(d) <= 1e-12).sum())} | "
                   f"{int((d > 1e-12).sum())} | {100*d.mean():+.4f} |")

    out.append("\n## Development replicates 01-05 (pooled; not used for claims)\n")
    out.append("| data | genes | candidate | baseline | n | Δerr (pp) | W/T/L |")
    out.append("|---|---|---|---|---|---|---|")
    dev = defaultdict(list)
    for (ds, ng, cond), lst in by.items():
        dev[(ds, ng)] += [r for rep, r in lst if rep in DEV]
    for (ds, ng) in sorted(dev):
        for c, b in PAIRS[:8]:
            s = paired(dev[(ds, ng)], c, b)
            if s:
                out.append(f"| {ds} | {ng} | {c} | {b} | {s['n']} | {100*s['mean']:+.2f} | {s['w']}/{s['t']}/{s['l']} |")

    out.append("\n## Mean runtime (s, one thread) per method, all replicates\n")
    out.append("| data | genes | " + " | ".join(METHODS) + " |")
    out.append("|" + "---|" * (2 + len(METHODS)))
    allg = defaultdict(list)
    for (ds, ng, cond), lst in by.items():
        allg[(ds, ng)] += [r for rep, r in lst]
    for (ds, ng) in sorted(allg):
        cells = []
        for m in METHODS:
            v = [r[m]["sec"] for r in allg[(ds, ng)] if m in r]
            cells.append(f"{np.mean(v):.1f}" if v else "-")
        out.append(f"| {ds} | {ng} | " + " | ".join(cells) + " |")
    out.append(f"\nReplicates loaded: {len(rows)}")
    txt = "\n".join(out) + "\n"
    open(os.path.join(RES, "tables.md"), "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
