# AI-assisted (Claude), exploration code for CS581 project
"""Run vote2.py variants on a replicate, score with FastSP vs the reference, append to REP/results2.jsonl.

    python3 run2.py REP VARIANT [VARIANT ...]

Backbones: REP/bb10 when REP/inputs/backbones holds 20 (h6b bank reps), else REP/inputs/backbones.
Restartable: variants already in results2.jsonl are skipped.
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, "..", "..", "code"))
sys.path.insert(0, CODE)
from gcmx.bbtool_bench import acc_ref  # noqa: E402


def bbdir(rep):
    src = os.path.join(rep, "inputs", "backbones")
    if len(os.listdir(src)) != 10 and os.path.isdir(os.path.join(rep, "bb10")):
        return os.path.join(rep, "bb10")
    return src


def main():
    rep = os.path.abspath(sys.argv[1])
    res = os.path.join(rep, "results2.jsonl")
    done = {json.loads(l)["variant"] for l in open(res)} if os.path.exists(res) else set()
    for v in sys.argv[2:]:
        if v in done:
            continue
        r = subprocess.run([sys.executable, os.path.join(HERE, "vote2.py"), rep, v, "--bb", bbdir(rep)],
                           capture_output=True, text=True)
        if r.returncode:
            print("FAILED", rep, v, r.stderr[-2000:], flush=True)
            continue
        info = json.loads(r.stdout.strip().splitlines()[-1])
        m = json.load(open(os.path.join(rep, "vote2", v, "model.json")))
        s = acc_ref(os.path.join(rep, "true.fasta"), info["out"])
        row = {"rep": os.path.basename(rep), "variant": v, "merge_wall": info["merge_wall"],
               "kept_edges": m["kept_edges_final"], "edges": m["edges"], "kept_weight_frac": m["kept_weight_frac"],
               "guard": m.get("guard"), "median_neff": m.get("median_neff"), **s}
        with open(res, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
