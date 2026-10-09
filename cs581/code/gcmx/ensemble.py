"""Multi-decomposition MAGUS ensemble + soft-constraint consensus.

    python -m gcmx.ensemble REP_DIR BASE_ALIGNMENT OUT_DIR [--subsets 10,25,50] [--threads 4]

REP_DIR holds inputs/backbones, unaligned.fasta and true.fasta (gcmx.prep /
fan-out layout). BASE_ALIGNMENT is MAGUS's own output for that replicate.

1. Re-estimate a tree from BASE_ALIGNMENT with FastTree (-nt -gtr).
2. For each K in --subsets: MAGUS with that tree, K subsets (paper's other
   settings), and the replicate's EXISTING backbones (-b), so only the cheap
   subset alignments are recomputed. K=25 alone = one extra MAGUS iteration.
3. GCM soft-constraint consensus (gcmx.consensus) of BASE + all K runs.
Every alignment is scored against true.fasta; results -> OUT_DIR/results.json.
"""

import argparse
import json
import os
import subprocess
import sys
import time

from . import consensus, score

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = ["--maxsubsetsize", "0", "--decompstrategy", "pastastyle", "--decompskeletonsize", "300",
         "--graphbuildmethod", "mafft", "--graphbuildhmmextend", "false", "--graphclustermethod", "mcl",
         "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false", "-f", "4"]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("rep")
    parser.add_argument("base")
    parser.add_argument("out")
    parser.add_argument("--subsets", default="10,25,50")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    true = os.path.join(args.rep, "true.fasta")
    results = {"base": score.fastsp(true, args.base)["avgErr"]}

    tree = os.path.join(args.out, "iter_tree.nwk")
    if not os.path.exists(tree):
        start = time.time()
        with open(tree, "w") as f:
            subprocess.run(["FastTree", "-nt", "-gtr", "-nosupport", "-quiet", args.base], stdout=f,
                           stderr=subprocess.DEVNULL, check=True)
        results["tree_seconds"] = round(time.time() - start, 1)

    members = [args.base]
    for k in [int(x) for x in args.subsets.split(",")]:
        out = os.path.join(args.out, "magus_k{}.fasta".format(k))
        if not os.path.exists(out):
            start = time.time()
            with open(os.path.join(args.out, "magus_k{}.log".format(k)), "w") as log:
                subprocess.run([sys.executable, "-m", "gcmx.run_magus", "-np", str(args.threads),
                                "-d", os.path.join(args.out, "work_k{}".format(k)),
                                "-i", os.path.join(args.rep, "unaligned.fasta"), "-t", tree,
                                "--maxnumsubsets", str(k), "-b", os.path.join(args.rep, "inputs", "backbones"),
                                "-o", out] + PAPER, cwd=CODE, stdout=log, stderr=subprocess.STDOUT, check=True)
            results["k{}_seconds".format(k)] = round(time.time() - start, 1)
        results["k{}".format(k)] = score.fastsp(true, out)["avgErr"]
        members.append(out)
        print(json.dumps(results), flush=True)

    for name, inputs in (("consensus_all", members), ("consensus_base_k25", [args.base, members[2]] if len(members) > 2 else None)):
        if not inputs:
            continue
        out = os.path.join(args.out, name + ".fasta")
        info = consensus.build(out, os.path.join(args.out, name), inputs)
        results[name] = score.fastsp(true, out)["avgErr"]
        results[name + "_primary"] = info["primary"]
    json.dump(results, open(os.path.join(args.out, "results.json"), "w"), indent=1)
    print(json.dumps(results), flush=True)


if __name__ == "__main__":
    main()
