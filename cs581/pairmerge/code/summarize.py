"""Tables for REPORT.md from results/merge.jsonl.

    python summarize.py ../results/merge.jsonl > ../results/merge_summary.md

Per (condition, dataset-group, merger): mean SPFN, SPFP, avgErr (=(SPFN+SPFP)/2),
cross-half FN/FP, seconds. Then paired comparisons of every merger against a
reference merger on avgErr: mean difference, W/T/L (tie band 0.001 = 0.1 point)
and a two-sided Wilcoxon signed-rank p.
"""

import collections
import json
import sys

import numpy as np
from scipy.stats import wilcoxon

TIE = 0.001
ORDER = ["opal", "muscle3", "mafft-merge", "gcm", "progdp"]


def group(dataset):
    return dataset.rsplit("_R", 1)[0]


def paired(rows, cond, a, b, key="avgErr"):
    x = {r["dataset"]: r[key] for r in rows if r["condition"] == cond and r["merger"] == a and key in r}
    y = {r["dataset"]: r[key] for r in rows if r["condition"] == cond and r["merger"] == b and key in r}
    common = sorted(set(x) & set(y))
    d = np.array([x[k] - y[k] for k in common])
    if not len(d):
        return None
    w, t, l = int((d < -TIE).sum()), int((abs(d) <= TIE).sum()), int((d > TIE).sum())
    p = wilcoxon(d).pvalue if np.any(d != 0) and len(d) >= 2 else float("nan")
    return len(d), d.mean(), w, t, l, p


def main(path):
    rows = [json.loads(l) for l in open(path)]
    errors = [r for r in rows if "error" in r]
    rows = [r for r in rows if "avgErr" in r]
    conds = sorted({r["condition"] for r in rows}, key=lambda c: c != "oracle")
    mergers = [m for m in ORDER if any(r["merger"] == m for r in rows)]
    for cond in conds:
        print("\n### Condition: {}\n".format(cond))
        groups = sorted({group(r["dataset"]) for r in rows if r["condition"] == cond})
        print("Mean (SPFN+SPFP)/2 in %, full alignment; n = replicates\n")
        print("| merger | " + " | ".join(groups) + " | all |")
        print("|---" * (len(groups) + 2) + "|")
        for m in mergers:
            cells = []
            for g in groups + ["all"]:
                v = [r["avgErr"] for r in rows if r["condition"] == cond and r["merger"] == m
                     and (g == "all" or group(r["dataset"]) == g)]
                cells.append("{:.2f} (n={})".format(100 * np.mean(v), len(v)) if v else "-")
            print("| {} | {} |".format(m, " | ".join(cells)))
        print("\nAll replicates pooled, in %: full-alignment SPFN/SPFP, cross-half pairs only, runtime\n")
        print("| merger | SPFN | SPFP | cross FN | cross FP | constraints kept | mean seconds |")
        print("|---|---|---|---|---|---|---|")
        for m in mergers:
            R = [r for r in rows if r["condition"] == cond and r["merger"] == m]
            if not R:
                continue
            f = lambda k: 100 * np.mean([r[k] for r in R])
            print("| {} | {:.2f} | {:.2f} | {:.2f} | {:.2f} | {}/{} | {:.1f} |".format(
                m, f("SPFN"), f("SPFP"), f("crossFN"), f("crossFP"),
                sum(bool(r.get("constraintsKept")) for r in R), len(R), np.mean([r["seconds"] for r in R])))
        ref = "opal"
        print("\nPaired vs `{}` on avgErr (negative = better than {}; W/T/L = better/tie/worse, tie band {} points)\n".format(
            ref, ref, 100 * TIE))
        print("| merger | n | mean diff (points) | W/T/L | Wilcoxon p |")
        print("|---|---|---|---|---|")
        for m in mergers:
            if m == ref:
                continue
            res = paired(rows, cond, m, ref)
            if res:
                n, md, w, t, l, p = res
                print("| {} | {} | {:+.2f} | {}/{}/{} | {:.2g} |".format(m, n, 100 * md, w, t, l, p))
        res = paired(rows, cond, "progdp", "gcm")
        if res:
            n, md, w, t, l, p = res
            print("\nprogdp vs gcm (same graph): n={}, mean diff {:+.2f} points, W/T/L {}/{}/{}, p={:.2g}".format(
                n, 100 * md, w, t, l, p))
    if errors:
        print("\nFailed runs: " + ", ".join("{}/{}/{}: {}".format(r["dataset"], r["condition"], r["merger"],
                                                                   r["error"][:80]) for r in errors))


if __name__ == "__main__":
    main(sys.argv[1])
