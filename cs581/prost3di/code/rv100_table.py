"""Markdown table of the RV100 runs (avg FastSP error %, lower is better)."""
import json, os, sys
cols = ["full-linsi", "full-linsi3di", "full-fm", "merge-A", "merge-B", "merge-C", "merge-D", "merge-E"]
print("| set (identity) | " + " | ".join(cols) + " | MAGUS published |")
print("|---" * (len(cols) + 2) + "|")
pub = {"BBA0067": 25.6, "BBA0101": 28.5, "BBA0117": 12.1, "BBA0081": 57.0}
ident = {"BBA0067": 0.219, "BBA0101": 0.222, "BBA0117": 0.246, "BBA0081": 0.139}
for s in sys.argv[1:]:
    f = f"/opt/p3d/rv100runs/{s}/results.json"
    if not os.path.exists(f):
        continue
    d = json.load(open(f))
    cells = [f"{d[c]['avgErr']*100:.2f}" if c in d else "" for c in cols]
    print(f"| {s} ({ident[s]:.2f}) | " + " | ".join(cells) + f" | {pub[s]} |")
