"""Markdown tables (mean over replicates) from runtrees.py output.

    python summarize.py RESULTS.jsonl [--metric fn_rate|rf_rate|seconds|lnl_tool]
"""
import argparse
import json
from collections import defaultdict

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--metric", default="fn_rate")
    a = ap.parse_args()
    rows = [json.loads(l) for f in a.files for l in open(f)]
    cell = defaultdict(list)
    datasets, cols = [], []
    for r in rows:
        key = (r["aln"], r["method"])
        if r["dataset"] not in datasets:
            datasets.append(r["dataset"])
        if key not in cols:
            cols.append(key)
        cell[(r["dataset"],) + key].append(r.get(a.metric))
    print("| alignment | method | " + " | ".join(datasets) + " |")
    print("|---|---|" + "---|" * len(datasets))
    for aln, m in sorted(cols):
        vals = []
        for ds in datasets:
            v = cell.get((ds, aln, m))
            if not v or any(x is None for x in v):
                vals.append("-")
            elif a.metric.endswith("rate"):
                vals.append("%.1f%% (n=%d)" % (100 * np.mean(v), len(v)))
            else:
                vals.append("%.0f (n=%d)" % (np.mean(v), len(v)))
        print("| %s | %s | %s |" % (aln, m, " | ".join(vals)))


if __name__ == "__main__":
    main()
