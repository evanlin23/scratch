"""Mean species-tree FN (30 taxa, 27 internal branches) vs number of true gene families.
Usage: python analyze_curve.py curve.jsonl out.md"""
import json
import sys

import numpy as np

rows = {}
for l in open(sys.argv[1]):
    r = json.loads(l)
    rows[(r["setting"], r["rep"], r["nfam"])] = r  # keep last (re-runs overwrite)
ms = ["multi-mean", "pro-mean-truetags", "pro-mean-truetags-countroot", "pro-mean", "pro-min", "astrid-disco", "astral-pro"]
L = ["# Pure-GDL true gene trees, 30-taxon Yule species trees (height 1): mean FN vs #families\n",
     "Settings (lam, mu): critical = (2,2) on all branches; supercrit = (2,1); adversarial = internal (3,0.5), "
     "terminal (0,3). Families with >= 4 species. FN out of 27 internal branches; mean over replicates.\n",
     "| setting | families | reps | mean leaves | " + " | ".join(ms) + " |", "|---|---|---|---|" + "---|" * len(ms)]
for s in ["critical", "supercrit", "adversarial"]:
    for k in sorted({k for (ss, _, k) in rows if ss == s}):
        sub = [r for (ss, _, kk), r in rows.items() if ss == s and kk == k]
        vals = []
        for m in ms:
            v = [r["methods"][m]["FN"] for r in sub if m in r["methods"]]
            vals.append("%.2f" % np.mean(v) if v else "-")
        L.append("| %s | %d | %d | %.0f | " % (s, k, len(sub), np.mean([r["mean_leaves"] for r in sub])) + " | ".join(vals) + " |")
open(sys.argv[2], "w").write("\n".join(L) + "\n")
print("\n".join(L))
