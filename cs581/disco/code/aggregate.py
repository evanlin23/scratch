"""Aggregate run_rep.py results into markdown tables + a CSV.

python aggregate.py RUNROOT OUTPREFIX [--holdout 06,07,08,09,10]
"""
import argparse
import glob
import json
import os

import numpy as np
from scipy.stats import wilcoxon

PAIRS = [("ASTRID-DISCOR", "ASTRID-DISCO"), ("ASTRAL-DISCOR", "ASTRAL-DISCO"),
         ("ASTRID-DISCOR", "ASTRAL-Pro2"), ("ASTRAL-DISCOR", "ASTRAL-Pro2"),
         ("ASTRID-DISCO", "ASTRAL-Pro2"),
         ("ASTRID-DISCOR-it2", "ASTRID-DISCO"), ("ASTRAL-DISCOR-it2", "ASTRAL-DISCO"),
         ("ASTRID-DISCOR-oracle", "ASTRID-DISCO"), ("ASTRAL-DISCOR-oracle", "ASTRAL-DISCO"),
         ("ASTRID-DISCOR-lca", "ASTRID-DISCO"), ("ASTRAL-DISCOR-lca", "ASTRAL-DISCO")]
TIE = 1e-9


def load(root):
    rows = []
    for f in glob.glob(os.path.join(root, "**", "result.json"), recursive=True):
        d = json.load(open(f))
        parts = os.path.normpath(d["rep"]).split(os.sep)
        d["cond"], d["repid"] = parts[-2], parts[-1]
        d["g"] = d["genes"].replace(".trees", "")
        rows.append(d)
    return rows


def compare(rows, a, b, key="rf"):
    x = [(r["rf"][a][key], r["rf"][b][key]) for r in rows if a in r["rf"] and b in r["rf"]]
    if not x:
        return None
    x = np.array(x)
    d = x[:, 0] - x[:, 1]
    w, t, l = int((d < -TIE).sum()), int((abs(d) <= TIE).sum()), int((d > TIE).sum())
    p = wilcoxon(d[abs(d) > TIE]).pvalue if (abs(d) > TIE).sum() >= 1 else 1.0
    return {"n": len(d), "mean_a": x[:, 0].mean(), "mean_b": x[:, 1].mean(), "diff": d.mean(),
            "W": w, "T": t, "L": l, "p": p}


def fmt_cmp(c):
    if c is None:
        return "–"
    return f"{c['diff']:+.4f} | {c['W']}/{c['T']}/{c['L']} | {c['p']:.3g}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("out")
    ap.add_argument("--holdout", default="06,07,08,09,10")
    a = ap.parse_args()
    rows = load(a.root)
    hold = set(a.holdout.split(","))
    lines = []
    # CSV of RF
    with open(a.out + "_rf.csv", "w") as f:
        f.write("cond,rep,genes,method,rf,fn,fp\n")
        for r in sorted(rows, key=lambda r: (r["cond"], r["repid"], r["g"])):
            for m, v in r["rf"].items():
                f.write(f"{r['cond']},{r['repid']},{r['g']},{m},{v['rf']:.5f},{v['fn']:.5f},{v['fp']:.5f}\n")
    gsets = sorted(set(r["g"] for r in rows))
    conds = sorted(set(r["cond"] for r in rows))
    methods = ["ASTRAL-Pro2", "ASTRID-DISCO", "ASTRAL-DISCO", "ASTRID-DISCOR", "ASTRAL-DISCOR",
               "ASTRID-DISCOR-it2", "ASTRAL-DISCOR-it2", "ASTRID-DISCOR-lca", "ASTRAL-DISCOR-lca",
               "ASTRID-DISCOR-oracle", "ASTRAL-DISCOR-oracle"]
    for g in gsets:
        R = [r for r in rows if r["g"] == g]
        lines.append(f"\n### Mean species-tree RF, gene trees = {g}\n")
        lines.append("| condition | n | " + " | ".join(m.replace("ASTRID-", "AD-").replace("ASTRAL-", "AL-") for m in methods) + " |")
        lines.append("|---|---|" + "---|" * len(methods))
        for c in conds:
            rc = [r for r in R if r["cond"] == c]
            if not rc:
                continue
            vals = []
            for m in methods:
                v = [r["rf"][m]["rf"] for r in rc if m in r["rf"]]
                vals.append(f"{np.mean(v):.4f}" if v else "–")
            lines.append(f"| {c} | {len(rc)} | " + " | ".join(vals) + " |")
        for label, sub in [("all replicates", R), ("held-out replicates " + a.holdout, [r for r in R if r["repid"] in hold])]:
            lines.append(f"\n**Paired comparisons, {g}, {label}** (diff = mean RF(A) − RF(B); W/T/L = A better/tie/worse, tie = identical RF; two-sided Wilcoxon signed-rank on non-tied pairs)\n")
            lines.append("| A vs B | n | mean A | mean B | diff | W/T/L | p |")
            lines.append("|---|---|---|---|---|---|---|")
            for A, B in PAIRS:
                c = compare(sub, A, B)
                if c:
                    lines.append(f"| {A} vs {B} | {c['n']} | {c['mean_a']:.4f} | {c['mean_b']:.4f} | {c['diff']:+.4f} | {c['W']}/{c['T']}/{c['L']} | {c['p']:.3g} |")
        # per-condition key comparisons
        lines.append(f"\n**Per condition, {g}** (diff | W/T/L | p)\n")
        lines.append("| condition | AD-DISCOR vs AD-DISCO | AL-DISCOR vs AL-DISCO | AD-DISCOR vs Pro2 | AD-oracle vs AD-DISCO |")
        lines.append("|---|---|---|---|---|")
        for c in conds:
            rc = [r for r in R if r["cond"] == c]
            if rc:
                lines.append(f"| {c} | {fmt_cmp(compare(rc, 'ASTRID-DISCOR', 'ASTRID-DISCO'))} | {fmt_cmp(compare(rc, 'ASTRAL-DISCOR', 'ASTRAL-DISCO'))} | {fmt_cmp(compare(rc, 'ASTRID-DISCOR', 'ASTRAL-Pro2'))} | {fmt_cmp(compare(rc, 'ASTRID-DISCOR-oracle', 'ASTRID-DISCO'))} |")
    # tagging
    T = [r for r in rows if "tags" in r]
    if T:
        lines.append("\n### Tagging accuracy (pooled over replicates; pair acc = fraction of cross-species gene pairs whose ortholog/paralog call matches the locus tree)\n")
        lines.append("| condition | genes | DISCO pair acc | DISCO-R pair acc | DISCO-R-lca | oracle-S | DISCO orth P/R | DISCO-R orth P/R | DISCO root acc* | DISCO-R root acc* |")
        lines.append("|---|---|---|---|---|---|---|---|---|---|")
        for c in conds:
            for g in gsets:
                rc = [r for r in T if r["cond"] == c and r["g"] == g]
                if not rc:
                    continue

                def m(meth, k):
                    v = [r["tags"][meth].get(k) for r in rc if k in r["tags"][meth]]
                    return f"{np.mean(v):.3f}" if v else "–"
                lines.append(f"| {c} | {g} | {m('DISCO','acc')} | {m('DISCOR','acc')} | {m('DISCOR-lca','acc')} | {m('DISCOR-oracle','acc')} | "
                             f"{m('DISCO','orth_prec')}/{m('DISCO','orth_rec')} | {m('DISCOR','orth_prec')}/{m('DISCOR','orth_rec')} | "
                             f"{m('DISCO','root_acc')} | {m('DISCOR','root_acc')} |")
        lines.append("\n\\* root acc: fraction of multi-copy true gene trees (input rooting randomised) whose chosen root bipartition equals SimPhy's.")
    # runtime
    lines.append("\n### Runtime (mean seconds per replicate, 1 thread)\n")
    keys = ["astral-pro2", "DISCO-decomp", "DISCO-astrid", "DISCO-astral", "s0-root", "DISCOR-decomp", "DISCOR-astrid", "DISCOR-astral"]
    lines.append("| condition | genes | " + " | ".join(keys) + " |")
    lines.append("|---|---|" + "---|" * len(keys))
    for c in conds:
        for g in gsets:
            rc = [r for r in rows if r["cond"] == c and r["g"] == g]
            if rc:
                lines.append(f"| {c} | {g} | " + " | ".join(f"{np.mean([r['time'].get(k, np.nan) for r in rc]):.1f}" for k in keys) + " |")
    # S0 rooting
    s0 = [r.get("s0_root_correct") for r in rows if "s0_root_correct" in r]
    lines.append(f"\nS0 root correct (min-DL rooting of ASTRAL-Pro2 tree): {sum(s0)}/{len(s0)} runs.")
    open(a.out + ".md", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
