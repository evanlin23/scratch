"""MAGUS-lite: no MAFFT backbones; evidence = alignments MAGUS already has, merged with soft constraints.

    python -m gcmx.lite TRUE_ALIGNMENT OUT_DIR NUM_SUBSETS [--threads 4] [--groups 3] [--rounds 1]

MAGUS spends ~55-60% of its runtime on 10 MAFFT-linsi backbone alignments
that only serve as evidence for the alignment graph. MAGUS-lite drops them:

  1. MAGUS with --graphbuildmethod initial (MAGUS's own option): the guide-tree
     alignment (MAFFT-linsi skeleton of 300 + HMM-added rest), which MAGUS
     computes anyway, is the only evidence. Output A0.
  2. Self-derived evidence: A0 restricted to 8 random sequences per subset
     (10 pseudo-backbones), each HMM-extended to all sequences.
  3. Soft constraints: every subset alignment split into GROUPS similarity
     groups; GCM merges all groups at once with evidence = initial alignment +
     self-derived evidence. Output A1. With --rounds R > 1, steps 2-3 repeat
     from the latest output.

Same decomposition, subset aligner (MAFFT-linsi) and GCM settings as the
paper's MAGUS. Records wall-clock per stage and SPFN/SPFP of every output.
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time

from . import fasta, score
from .pilot import self_backbones

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def flags(k):
    return ["--maxsubsetsize", "0", "--maxnumsubsets", str(k), "--decompstrategy", "pastastyle",
            "--decompskeletonsize", "300", "--graphclustermethod", "mcl", "--graphtracemethod", "minclusters",
            "--graphtraceoptimize", "false", "-r", "10", "-m", "200", "-f", "4"]


def run(cmd, log):
    start = time.time()
    with open(log, "w") as f:
        subprocess.run([sys.executable, "-m"] + cmd, cwd=CODE, stdout=f, stderr=subprocess.STDOUT, check=True)
    return round(time.time() - start, 1)


def stage_seconds(magus_log):
    text = open(magus_log).read()
    get = lambda pattern: float(re.findall(pattern, text)[-1]) if re.findall(pattern, text) else None
    return {"initial_tree": get(r"Built initial tree on .* in ([0-9.]+) sec"),
            "merge": get(r"Merged \d+ subalignments into .* in ([0-9.]+) sec")}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("true")
    parser.add_argument("out")
    parser.add_argument("k")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--groups", type=int, default=3)
    parser.add_argument("--rounds", type=int, default=1)
    parser.add_argument("--mclthreads", type=int, default=0)
    args = parser.parse_args()
    T = str(args.threads)
    os.makedirs(args.out, exist_ok=True)
    true = os.path.join(args.out, "true.fasta")
    unaligned = os.path.join(args.out, "unaligned.fasta")
    ref = fasta.upper(fasta.read(args.true))
    fasta.write(ref, true)
    fasta.write(fasta.ungap(ref), unaligned)
    acc = lambda path: {k: score.fastsp(true, path)[k] for k in ("SPFN", "SPFP", "avgErr")}
    mcl = ["--gcmx-mclthreads", str(args.mclthreads)]
    row = {"true": args.true, "threads": args.threads, "groups": args.groups, "mclthreads": args.mclthreads}

    # 1. MAGUS, initial alignment as the only evidence
    work = os.path.join(args.out, "magus_initial")
    shutil.rmtree(work, ignore_errors=True)
    a0 = os.path.join(args.out, "A0.fasta")
    wall = run(["gcmx.run_magus", "-np", T, "-d", work, "-i", unaligned, "-o", a0, "--graphbuildmethod", "initial"]
               + mcl + flags(args.k), os.path.join(args.out, "A0.log"))
    row["A0"] = {"wall": wall, **stage_seconds(os.path.join(work, "log.txt")), **acc(a0)}
    print(json.dumps({"A0": row["A0"]}), flush=True)
    subs = os.path.join(work, "subalignments")
    initial = os.path.join(work, "decomposition", "initial_tree", "initial_insert_align.txt")

    # 2-3. soft-constraint merge with initial + self-derived evidence
    split = os.path.join(args.out, "split")
    row["split_wall"] = run(["gcmx.split", subs, split, str(args.groups)], os.path.join(args.out, "split.log"))
    current = a0
    for r in range(1, args.rounds + 1):
        start = time.time()
        bb, ev = os.path.join(args.out, "self_r{}".format(r)), os.path.join(args.out, "ev_r{}".format(r))
        self_backbones(current, subs, bb, seed=100 + r)
        subprocess.run([sys.executable, "-m", "gcmx.extend", bb, unaligned, ev, "--jobs", T], check=True, cwd=CODE,
                       stdout=subprocess.DEVNULL)
        shutil.copy(initial, os.path.join(ev, "initial_insert_align.txt"))
        ev_wall = round(time.time() - start, 1)
        out = os.path.join(args.out, "A{}.fasta".format(r))
        mwork = os.path.join(args.out, "merge_r{}".format(r))
        shutil.rmtree(mwork, ignore_errors=True)
        m_wall = run(["gcmx.run_magus", "-np", T, "-d", mwork, "-s", split, "-b", ev, "-o", out] + mcl,
                     os.path.join(args.out, "A{}.log".format(r)))
        row["A{}".format(r)] = {"evidence_wall": ev_wall, "merge_wall": m_wall, **acc(out)}
        print(json.dumps({"A{}".format(r): row["A{}".format(r)]}), flush=True)
        current = out
        shutil.rmtree(mwork, ignore_errors=True)  # graph files are large
    row["total_wall"] = round(row["A0"]["wall"] + row["split_wall"] + sum(
        row["A{}".format(r)]["evidence_wall"] + row["A{}".format(r)]["merge_wall"] for r in range(1, args.rounds + 1)), 1)
    with open(os.path.join(args.out, "lite.json"), "w") as f:
        json.dump(row, f, indent=1)
    print(json.dumps({"total_wall": row["total_wall"]}), flush=True)


if __name__ == "__main__":
    main()
