"""Summarize merge-variant results across replicates (Markdown).

    python -m gcmx.summarize RUNS_DIR > SUMMARY.md

RUNS_DIR/<CONDITION>_R<k>/pilot.jsonl rows come from gcmx.pilot: every
variant was run on the *same* subalignments and backbones as the default
MAGUS merge of that replicate, so differences are paired.

Tables:
  1. mean average error (SPFN+SPFP)/2 in % per variant and model condition
  2. paired comparison with `default`: mean difference, wins/ties/losses
     (|diff| <= 0.05 points is a tie) and a two-sided Wilcoxon signed-rank p
  3. oracle decomposition: how much error remains with true backbones / true
     subset alignments
"""

import collections
import glob
import json
import os
import sys

from scipy.stats import wilcoxon

ROSE = ["1000L3", "1000S1", "1000M2", "RNASim", "1000L1", "1000S2", "1000S3", "1000M3", "1000L2", "1000M4"]


def group(dataset):
    cond = dataset.rsplit("_R", 1)[0]
    return "BAliBASE" if cond.startswith("BBA") else cond


def pct(values):
    return "{:.2f}".format(100 * sum(values) / len(values)) if values else "–"


def main():
    rows = []
    for path in sorted(glob.glob(os.path.join(sys.argv[1], "*_R*", "pilot.jsonl"))):
        rows += [r for r in map(json.loads, open(path)) if "avgErr" in r]
    if not rows:
        print("no results yet")
        return
    groups = [g for g in ROSE + ["BAliBASE", "16S.M", "16S.T", "16S.3", "RNASim10K"]
              if any(group(r["dataset"]) == g for r in rows)]
    variants = list(dict.fromkeys(r["variant"] for r in rows))
    by = collections.defaultdict(list)
    for r in rows:
        by[r["variant"], group(r["dataset"])].append(r["avgErr"])
        by[r["variant"], "all"].append(r["avgErr"])

    print("### 1. Average error (%) per model condition (replicates in parentheses)\n")
    print("| variant | " + " | ".join(groups) + " | all |")
    print("|---|" + "---|" * (len(groups) + 1))
    for v in variants:
        cells = ["{} ({})".format(pct(by[v, g]), len(by[v, g])) if by[v, g] else "–" for g in groups + ["all"]]
        print("| {} | {} |".format(v, " | ".join(cells)))

    base = {r["dataset"]: r for r in rows if r["variant"] == "default"}
    print("\n### 2. Paired comparison with default MAGUS merge (same inputs)\n")
    print("Negative difference = lower error than MAGUS. W/T/L counts replicates (tie: |diff| <= 0.05 points).\n")
    print("| variant | n | mean diff (pts) | W/T/L | Wilcoxon p | mean length ratio | mean seconds |")
    print("|---|---|---|---|---|---|---|")
    for v in variants:
        if v == "default" or v.startswith("oracle"):
            continue
        pairs = [(r, base[r["dataset"]]) for r in rows if r["variant"] == v and r["dataset"] in base]
        if not pairs:
            continue
        diffs = [100 * (r["avgErr"] - b["avgErr"]) for r, b in pairs]
        w = sum(d < -0.05 for d in diffs)
        l = sum(d > 0.05 for d in diffs)
        try:
            p = "{:.3g}".format(wilcoxon(diffs).pvalue) if any(diffs) and len(diffs) >= 5 else "–"
        except ValueError:
            p = "–"
        ratios = [r["LenEst"] / r["LenRef"] for r, _ in pairs if r.get("LenRef")]
        print("| {} | {} | {:+.2f} | {}/{}/{} | {} | {:.2f} | {:.0f} |".format(
            v, len(diffs), sum(diffs) / len(diffs), w, len(diffs) - w - l, l, p,
            sum(ratios) / len(ratios) if ratios else float("nan"), sum(r["seconds"] for r, _ in pairs) / len(pairs)))

    print("\n### 3. Oracle decomposition (average error %, same subsets)\n")
    print("| condition | MAGUS | true backbones | true full alignment as backbone | true subset alignments | split m=3 + true full backbone |")
    print("|---|---|---|---|---|---|")
    for g in groups:
        cells = [pct(by[v, g]) for v in ("default", "oracle-estSub-trueBB", "oracle-estSub-trueFull",
                                         "oracle-trueSub-estBB", "oracle-soft-m3-trueFull")]
        print("| {} | {} |".format(g, " | ".join(cells)))


if __name__ == "__main__":
    main()
