"""Controlled runtime benchmark: MAGUS(Fast), MAGUS(Slow) and Slow+soft, end to end.

    python -m gcmx.timing_bench JOBFILE OUT.jsonl WORKDIR [--threads 4] [--groups 3]

JOBFILE lines: NAME TRUE_ALIGNMENT NUM_SUBSETS (as in fanout/jobs_*.txt).
Run on an otherwise idle machine. For every job, sequentially and with the
same thread count:
  fast       full MAGUS with the paper's MAGUS(Fast) flags
  slow       full MAGUS with the paper's MAGUS(Slow) flags (--graphbuildhmmextend true)
  slow-soft  MAGUS(Fast)'s decomposition, subset alignments and backbones, then
             HMM extension + similarity split + soft-constraint merge (this project)
Wall-clock and CPU seconds (children's user+sys) are recorded per stage, plus
accuracy (SPFN/SPFP via FastSP). slow-soft total = fast total - fast merge
time + extension + split + soft merge. fast/slow run MAGUS's original code
(--gcmx-fastgraph false); slow-soft needs gcmx.fastgraph's compact graph
(MAGUS's dict graph runs out of memory for soft constraints), which also makes
graph construction faster -- state this when comparing runtimes.
"""

import argparse
import glob
import json
import os
import re
import resource
import shutil
import subprocess
import sys
import time

from . import fasta, score

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def flags(k):
    return ["--maxsubsetsize", "0", "--maxnumsubsets", str(k), "--decompstrategy", "pastastyle",
            "--decompskeletonsize", "300", "--graphbuildmethod", "mafft", "--graphclustermethod", "mcl",
            "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false", "-r", "10", "-m", "200", "-f", "4"]


def timed(cmd, log):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.time()
    with open(log, "w") as f:
        subprocess.run(cmd, cwd=CODE, stdout=f, stderr=subprocess.STDOUT, check=True)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    return round(time.time() - start, 1), round(cpu, 1)


def merge_seconds(magus_log):
    text = open(magus_log).read()
    found = re.findall(r"Merged \d+ subalignments into .* in ([0-9.]+) sec", text)
    return float(found[-1]) if found else 0.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("jobs")
    parser.add_argument("out")
    parser.add_argument("work")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--groups", type=int, default=3)
    args = parser.parse_args()
    done = set()
    if os.path.exists(args.out):
        done = {json.loads(l)["dataset"] for l in open(args.out)}
    py = [sys.executable, "-m"]
    T = str(args.threads)

    for line in open(args.jobs):
        if not line.strip():
            continue
        name, src, k = line.split()[:3]
        if name in done:
            continue
        repo = os.path.dirname(os.path.dirname(CODE))
        src = src if os.path.isabs(src) else os.path.join(repo, src)
        w = os.path.join(args.work, name)
        shutil.rmtree(w, ignore_errors=True)
        os.makedirs(w)
        true = fasta.upper(fasta.read(src))
        fasta.write(true, os.path.join(w, "true.fasta"))
        fasta.write(fasta.ungap(true), os.path.join(w, "unaligned.fasta"))
        row = {"dataset": name, "threads": args.threads, "nproc": os.cpu_count()}

        for mode, extra in (("fast", ["--graphbuildhmmextend", "false"]), ("slow", ["--graphbuildhmmextend", "true"])):
            out = os.path.join(w, mode + ".fasta")
            # MAGUS's own code path (original graph builder), i.e. the tool the paper timed
            wall, cpu = timed(py + ["gcmx.run_magus", "--gcmx-fastgraph", "false", "-np", T, "-d", os.path.join(w, mode), "-i",
                                    os.path.join(w, "unaligned.fasta"), "-o", out] + flags(k) + extra,
                              os.path.join(w, mode + ".log"))
            row[mode] = {"wall": wall, "cpu": cpu, "merge_wall": merge_seconds(os.path.join(w, mode, "log.txt")),
                         **{key: score.fastsp(os.path.join(w, "true.fasta"), out)[key] for key in ("SPFN", "SPFP", "avgErr")}}
            print(json.dumps({name: row[mode]}), flush=True)

        fw = os.path.join(w, "fast")
        subs, bbs = os.path.join(w, "soft_inputs", "subalignments"), os.path.join(w, "soft_inputs", "backbones")
        os.makedirs(subs)
        os.makedirs(bbs)
        for p in glob.glob(os.path.join(fw, "subalignments", "subalignment_subset_*.txt")):
            shutil.copy(p, subs)
        for p in glob.glob(os.path.join(fw, "graph", "backbone_*_mafft.txt")):
            shutil.copy(p, bbs)
        ext, split = os.path.join(w, "ext"), os.path.join(w, "split")
        e_wall, e_cpu = timed(py + ["gcmx.extend", bbs, os.path.join(w, "unaligned.fasta"), ext, "--jobs", T],
                              os.path.join(w, "extend.log"))
        s_wall, s_cpu = timed(py + ["gcmx.split", subs, split, str(args.groups)], os.path.join(w, "split.log"))
        out = os.path.join(w, "slow-soft.fasta")
        m_wall, m_cpu = timed(py + ["gcmx.run_magus", "-np", T, "-d", os.path.join(w, "softmerge"),
                                    "-s", split, "-b", ext, "-o", out], os.path.join(w, "softmerge.log"))
        shared_wall = row["fast"]["wall"] - row["fast"]["merge_wall"]
        row["slow-soft"] = {"wall": round(shared_wall + e_wall + s_wall + m_wall, 1),
                            "extend_wall": e_wall, "split_wall": s_wall, "merge_wall": m_wall,
                            "extra_cpu": round(e_cpu + s_cpu + m_cpu, 1),
                            **{key: score.fastsp(os.path.join(w, "true.fasta"), out)[key] for key in ("SPFN", "SPFP", "avgErr")}}
        print(json.dumps({name: row["slow-soft"]}), flush=True)
        with open(args.out, "a") as f:
            f.write(json.dumps(row) + "\n")
        for d in ("fast", "slow", "softmerge", "ext", "split", "soft_inputs"):
            shutil.rmtree(os.path.join(w, d), ignore_errors=True)  # keep only outputs and logs


if __name__ == "__main__":
    main()
