"""Iterative soft refinement: derive evidence from the current alignment, re-merge, repeat.

    python -m gcmx.iterate REP_DIR START_ALIGNMENT SPLIT_DIR OUT_DIR [--rounds 3] [--jobs 2]

Round r: pseudo-backbones = current alignment restricted to 8 random sequences
per subset (new random draw each round), HMM-extended to all sequences; then a
soft-constraint GCM merge of SPLIT_DIR's groups with that evidence. The output
becomes the next round's current alignment. Every round is scored against
REP_DIR/true.fasta (for analysis only; the procedure never looks at it).
"""

import argparse
import json
import os
import subprocess
import sys
import time

from . import score
from .pilot import self_backbones

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("rep")
    parser.add_argument("start")
    parser.add_argument("split")
    parser.add_argument("out")
    parser.add_argument("--rounds", type=int, default=3)
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    true = os.path.join(args.rep, "true.fasta")
    current = args.start
    history = [{"round": 0, "avgErr": score.fastsp(true, current)["avgErr"]}]
    print(json.dumps(history[-1]), flush=True)
    for r in range(1, args.rounds + 1):
        start = time.time()
        bb, ext = os.path.join(args.out, "bb_r{}".format(r)), os.path.join(args.out, "ext_r{}".format(r))
        self_backbones(current, os.path.join(args.rep, "inputs", "subalignments"), bb, seed=100 + r)
        subprocess.run([sys.executable, "-m", "gcmx.extend", bb, os.path.join(args.rep, "unaligned.fasta"), ext,
                        "--jobs", str(args.jobs)], check=True, cwd=CODE, stdout=subprocess.DEVNULL)
        out = os.path.join(args.out, "round_{}.fasta".format(r))
        with open(os.path.join(args.out, "round_{}.log".format(r)), "w") as log:
            subprocess.run([sys.executable, "-m", "gcmx.run_magus", "-np", str(args.jobs), "-d",
                            os.path.join(args.out, "work_r{}".format(r)), "-s", args.split, "-b", ext, "-o", out],
                           check=True, cwd=CODE, stdout=log, stderr=subprocess.STDOUT)
        s = score.fastsp(true, out)
        history.append({"round": r, "avgErr": s["avgErr"], "SPFN": s["SPFN"], "SPFP": s["SPFP"],
                        "LenEst": s["LenEst"], "seconds": round(time.time() - start, 1)})
        print(json.dumps(history[-1]), flush=True)
        current = out
        subprocess.run(["rm", "-rf", os.path.join(args.out, "work_r{}".format(r))])
    json.dump(history, open(os.path.join(args.out, "history.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
