"""Generate additional MAGUS-style backbones for a replicate.

    python -m gcmx.more_backbones SUBALIGNMENT_DIR OUT_DIR N [--per-subset 8] [--seed 1] [--jobs 2]

Same scheme as MAGUS (graph_builder.assignBackboneTaxa + external_tools.runMafft):
each backbone takes `--per-subset` random sequences from every subset
(200 / 25 = 8 by default) and aligns them with MAFFT-L-INS-i, MAGUS's command.
"""

import argparse
import concurrent.futures
import os
import random

from . import fasta
from .local_backbones import align


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("subalignments")
    parser.add_argument("outdir")
    parser.add_argument("n", type=int)
    parser.add_argument("--per-subset", type=int, default=8)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    subs = [fasta.read(os.path.join(args.subalignments, f)) for f in sorted(os.listdir(args.subalignments))]
    os.makedirs(args.outdir, exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        jobs = []
        for b in range(args.n):
            seqs = {}
            for sub in subs:
                for t in rng.sample(list(sub), min(args.per_subset, len(sub))):
                    seqs[t] = sub[t].replace("-", "").replace(".", "").upper()
            jobs.append(pool.submit(align, seqs, os.path.join(args.outdir, "extra_backbone_{}.txt".format(b + 1))))
        for job in jobs:
            print("wrote", job.result(), flush=True)


if __name__ == "__main__":
    main()
