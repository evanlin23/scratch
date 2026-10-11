# AI-assisted (Claude), exploration code for CS581 project
"""Parallel version of gcmvote/code/run.py (B=10 only): runs vote.py variants on a replicate in parallel and
scores each (FastSP vs REP/true.fasta, as gg.py / run.py do). Rows (run.py's format) go to REP/vote_results.jsonl,
kept apart from gg.py's REP/results.jsonl. Restartable: variants already in vote_results.jsonl are skipped.

    python3 par.py REP NPROC VARIANT [VARIANT ...]
"""
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from threading import Lock

HERE = os.path.dirname(os.path.abspath(__file__))
VOTE = os.path.join(HERE, "..", "code", "vote.py")
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "code")))
from gcmx.bbtool_bench import acc_ref  # noqa: E402

lock = Lock()


def one(rep, v, res):
    r = subprocess.run([sys.executable, VOTE, rep, v, "--bb", os.path.join(rep, "inputs", "backbones"), "--tag", v, "--threads", "1"],
                       capture_output=True, text=True)
    if r.returncode:
        print("FAILED", rep, v, r.stderr[-2000:], flush=True)
        return
    info = json.loads(r.stdout.strip().splitlines()[-1])
    m = json.load(open(os.path.join(rep, "vote", v, "model.json")))
    s = acc_ref(os.path.join(rep, "true.fasta"), info["out"])
    row = {"rep": os.path.basename(rep), "variant": v, "B": 10, "merge_wall": info["merge_wall"],
           "kept_edges": m["kept_edges"], "edges": m["edges"], "kept_weight_frac": m["kept_weight_frac"],
           "cutoff_by_n": m["cutoff_by_n"], **s}
    with lock, open(res, "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)


def main():
    rep, nproc, variants = os.path.abspath(sys.argv[1]), int(sys.argv[2]), sys.argv[3:]
    res = os.path.join(rep, "vote_results.jsonl")
    done = {json.loads(l)["variant"] for l in open(res)} if os.path.exists(res) else set()
    with ThreadPoolExecutor(nproc) as ex:
        for v in variants:
            if v not in done:
                ex.submit(one, rep, v, res)


if __name__ == "__main__":
    main()
