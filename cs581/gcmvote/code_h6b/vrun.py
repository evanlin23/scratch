# AI-assisted (Claude), exploration code for CS581 project
"""Vote-model variants (../code/vote.py, copied from claude/cs581-gcmvote) on an h6 replicate with B backbones.

    python3 vrun.py REP B VARIANT [...]   -> REP/vote_results.jsonl (restartable on (variant, B))

REP/inputs/backbones holds MAGUS's 20 L-INS-i backbones; B uses backbones 1..B (REP/bbB), the same
backbones gg.py's linsi~B variants use."""
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VOTE = os.path.join(HERE, "..", "code", "vote.py")
sys.path.insert(0, os.path.join(HERE, "..", "..", "code"))
from gcmx.bbtool_bench import acc_ref  # noqa: E402


def bbdir(rep, B):
    src = os.path.join(rep, "inputs", "backbones")
    d = os.path.join(rep, "bb{}".format(B))
    if not os.path.isdir(d):
        os.makedirs(d + ".tmp", exist_ok=True)
        for f in sorted(os.listdir(src), key=lambda f: int(f.split("_")[1]))[:B]:
            shutil.copy(os.path.join(src, f), d + ".tmp")
        os.rename(d + ".tmp", d)
    return d


rep, B = os.path.abspath(sys.argv[1]), int(sys.argv[2])
res = os.path.join(rep, "vote_results.jsonl")
done = {(r["variant"], r["B"]) for r in map(json.loads, open(res))} if os.path.exists(res) else set()
for v in sys.argv[3:]:
    if (v, B) in done:
        continue
    tag = "{}@B{}".format(v, B)
    r = subprocess.run([sys.executable, VOTE, rep, v, "--bb", bbdir(rep, B), "--tag", tag], capture_output=True, text=True)
    if r.returncode:
        print("FAILED", rep, tag, r.stderr[-2000:], flush=True)
        continue
    info = json.loads(r.stdout.strip().splitlines()[-1])
    m = json.load(open(os.path.join(rep, "vote", tag, "model.json")))
    s = acc_ref(os.path.join(rep, "true.fasta"), info["out"])
    row = {"rep": os.path.basename(rep), "variant": v, "B": B, "merge_wall": info["merge_wall"],
           "kept_edges": m.get("kept_edges"), "edges": m.get("edges"), **s}
    with open(res, "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)
