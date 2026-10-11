# AI-assisted (Claude), exploration code for CS581 project
"""One vote.py variant (B=10) on a replicate, scored as run.py / gg.py do (acc_ref vs REP/true.fasta).
Row -> REP/vote_rows/VARIANT.json (separate files so lanes can run in parallel; skipped if present).

    python3 vote_one.py REP VARIANT
"""
import json
import os
import subprocess
import sys

S = "/home/user/scratch/cs581"
sys.path.insert(0, os.path.join(S, "code"))
from gcmx.bbtool_bench import acc_ref  # noqa: E402

rep, v = os.path.abspath(sys.argv[1]), sys.argv[2]
dst = os.path.join(rep, "vote_rows", v + ".json")
if os.path.exists(dst):
    sys.exit(0)
os.makedirs(os.path.dirname(dst), exist_ok=True)
r = subprocess.run([sys.executable, os.path.join(S, "gcmvote", "code", "vote.py"), rep, v,
                    "--bb", os.path.join(rep, "inputs", "backbones"), "--tag", v], capture_output=True, text=True)
if r.returncode:
    print("FAILED", rep, v, r.stderr[-2000:], flush=True)
    sys.exit(1)
info = json.loads(r.stdout.strip().splitlines()[-1])
m = json.load(open(os.path.join(rep, "vote", v, "model.json")))
s = acc_ref(os.path.join(rep, "true.fasta"), info["out"])
row = {"rep": os.path.basename(rep), "variant": v, "B": 10, "merge_wall": info["merge_wall"],
       "kept_edges": m["kept_edges"], "edges": m["edges"], "kept_weight_frac": m["kept_weight_frac"],
       "cutoff_by_n": m["cutoff_by_n"], **s}
open(dst, "w").write(json.dumps(row) + "\n")
print(json.dumps(row), flush=True)
