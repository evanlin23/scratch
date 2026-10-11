# AI-assisted (Claude), exploration code for CS581 project
"""REP/results.jsonl of every /opt/work/h6b replicate -> results_h6b/sens.jsonl
(dataset, B, variant, SPFN, SPFP, avgErr from FastSP via acc_ref, merge time)."""
import json
import os
import re

W = "/opt/work/h6b/reps"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results_h6b", "sens.jsonl")
rows = []
for ds in sorted(os.listdir(W)):
    p = os.path.join(W, ds, "results.jsonl")
    if not os.path.exists(p):
        continue
    for r in map(json.loads, open(p)):
        v = r["variant"]
        m = re.match(r"linsi(?:~(\d+))?(.*)$", v)
        if m and (m.group(2) == "" or m.group(2).startswith("#es")):
            B = int(m.group(1) or 20)
            k = int(m.group(2)[3:]) if m.group(2) else None
            var = "magus" if k is None else "es{} ({:.2f}B)".format(k, k / B)
        else:  # vote variants: B from their own row if recorded, else nbb
            B, var = r.get("B", r.get("nbb")), v
        rows.append({"dataset": ds, "B": B, "variant": var, "gg_variant": v, "SPFN": r["SPFN"], "SPFP": r["SPFP"],
                     "avgErr": r["avgErr"], "merge_wall": r["merge_wall"]})
    p = os.path.join(W, ds, "vote_results.jsonl")
    for r in (map(json.loads, open(p)) if os.path.exists(p) else []):
        rows.append({"dataset": ds, "B": r["B"], "variant": "vote:" + r["variant"], "gg_variant": None,
                     "SPFN": r["SPFN"], "SPFP": r["SPFP"], "avgErr": r["avgErr"], "merge_wall": r["merge_wall"]})
with open(OUT, "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")
print(len(rows), "rows ->", OUT)
