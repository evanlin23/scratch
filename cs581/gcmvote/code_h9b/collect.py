# AI-assisted (Claude), exploration code for CS581 project
"""Gather /opt/work/h9b/reps/<DATASET>_B<B>/results.jsonl into cs581/gcmvote/results_h9b/sens.jsonl
(one row per dataset, B, variant; plus the end-to-end 20-backbone MAGUS run as variant magus-e2e)."""
import glob
import json
import os

W = "/opt/work/h9b"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results_h9b", "sens.jsonl")
rows = []
for res in sorted(glob.glob(os.path.join(W, "reps", "*_B*", "results.jsonl"))):
    ds, b = os.path.basename(os.path.dirname(res)).rsplit("_B", 1)
    for line in open(res):
        r = json.loads(line)
        rows.append({"dataset": ds, "B": int(b), "variant": r["variant"], "nbb": r["nbb"], "SPFN": r["SPFN"],
                     "SPFP": r["SPFP"], "avgErr": r["avgErr"], "merge_wall": r["merge_wall"]})
if os.path.exists(os.path.join(W, "magus.jsonl")):
    seen = set()
    for line in open(os.path.join(W, "magus.jsonl")):
        r = json.loads(line)
        if r["dataset"] in seen:  # a resumed bbtool_bench run appends its row again
            continue
        seen.add(r["dataset"])
        m = r["magus"]
        rows.append({"dataset": r["dataset"], "B": 20, "variant": "magus-e2e", "nbb": 20, "SPFN": m["SPFN"],
                     "SPFP": m["SPFP"], "avgErr": m["avgErr"], "wall": m["wall"]})
order = {"1000M2_R0": 0, "1000L1_R0": 1, "16S.M_R0": 2}
rows.sort(key=lambda r: (order.get(r["dataset"], 9), r["B"]))
with open(OUT, "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")
print(len(rows), "rows ->", OUT)
