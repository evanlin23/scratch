# AI-assisted (Claude), exploration code for CS581 project
"""Pre-registered vote-model variants (cs581/gcmvote/code/vote.py) on an h9 rep dir whose
inputs/backbones holds exactly the B backbones of that rep; scored with FastSP as gg.py does.
Appends to REP/results.jsonl (same file as gg.py; variant names prefixed 'vote:'). Restartable.

    python3 vote_run.py REP VARIANT [VARIANT ...]
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VOTE = os.path.join(HERE, "..", "code", "vote.py")
sys.path.insert(0, os.path.join(HERE, "..", "..", "code"))
from gcmx.bbtool_bench import acc_ref  # noqa: E402

rep = os.path.abspath(sys.argv[1])
res = os.path.join(rep, "results.jsonl")
done = {json.loads(l)["variant"] for l in open(res)} if os.path.exists(res) else set()
bb = os.path.join(rep, "inputs", "backbones")
nbb = len([f for f in os.listdir(bb) if f.endswith("_mafft.txt")])
for v in sys.argv[2:]:
    name = "vote:" + v
    if name in done:
        continue
    r = subprocess.run([sys.executable, VOTE, rep, v, "--bb", bb, "--tag", v], capture_output=True, text=True)
    if r.returncode:
        print("FAILED", rep, v, r.stderr[-2000:], flush=True)
        continue
    info = json.loads(r.stdout.strip().splitlines()[-1])
    m = json.load(open(os.path.join(rep, "vote", v, "model.json")))
    s = acc_ref(os.path.join(rep, "true.fasta"), info["out"])
    row = {"rep": os.path.basename(rep), "variant": name, "nbb": nbb, "merge_wall": info["merge_wall"],
           "kept_edges": m.get("kept_edges"), "edges": m.get("edges"), **s}
    with open(res, "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)
