"""Measured end-to-end runtime + accuracy: PASTA vs MAGUS vs soft MAGUS, same machine.

    python -m gcmx.e2e_bench JOBFILE OUT.jsonl WORKDIR [--threads 4]

JOBFILE lines: NAME TRUE_ALIGNMENT NUM_SUBSETS (as in fanout/jobs_*.txt).
Run on an otherwise idle machine. For every job, sequentially, all with the
same thread count, every pipeline run for real from unaligned sequences:

  pasta       PASTA default (3 iterations, its own starting tree), as in the MAGUS paper
  magus       MAGUS(Fast), the paper's flags, MAGUS's own code
  magus-slow  MAGUS(Slow), the paper's flags
  self-soft   magus, then: self-derived evidence (magus output restricted to 8 seqs
              per subset, 10 pseudo-backbones, HMM-extended) + 3 similarity groups
              per subset + soft-constraint GCM merge (MCL threaded)
  slow-soft   magus, then: HMM-extension of magus's 10 MAFFT backbones + 3 groups
              per subset + soft-constraint GCM merge (MCL threaded)

The soft pipelines' totals are magus's measured wall-clock plus the measured
wall-clock of their extra steps (they literally run after magus, on its files);
nothing is estimated. Accuracy is FastSP SPFN/SPFP against the reference.
"""

import argparse
import collections
import glob
import json
import os
import resource
import shutil
import subprocess
import sys
import time

from . import fasta, score
from .pilot import self_backbones

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = os.path.dirname(os.path.dirname(CODE))
PASTA_PY = "/opt/mm/root/envs/pasta183/bin/python"  # PASTA 1.8.3 (the MAGUS paper's version), see setup.sh
PASTA = "/opt/src/pasta/run_pasta.py"


def magus_flags(k):
    return ["--maxsubsetsize", "0", "--maxnumsubsets", str(k), "--decompstrategy", "pastastyle",
            "--decompskeletonsize", "300", "--graphbuildmethod", "mafft", "--graphclustermethod", "mcl",
            "--graphtracemethod", "minclusters", "--graphtraceoptimize", "false", "-r", "10", "-m", "200", "-f", "4"]


def timed(cmd, log, env=None):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.time()
    with open(log, "w") as f:
        proc = subprocess.run(cmd, cwd=CODE, stdout=f, stderr=subprocess.STDOUT, env=env)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = (after.ru_utime - before.ru_utime) + (after.ru_stime - before.ru_stime)
    if proc.returncode != 0:
        raise RuntimeError("failed ({}): {}".format(proc.returncode, " ".join(cmd)))
    return round(time.time() - start, 1), round(cpu, 1)


def seq_type(seqs):
    """dna / rna / protein by composition: nucleotide data may carry a few IUPAC ambiguity codes
    (16S.M has N, Y, R, W, ...), so a residue set check misclassifies it as protein."""
    counts = collections.Counter("".join(seqs.values()).upper())
    total = sum(v for k, v in counts.items() if k not in "-.") or 1
    if sum(counts[c] for c in "ACGTUN") / total < 0.9:
        return "protein"
    return "rna" if counts["U"] > counts["T"] else "dna"


def is_protein(seqs):
    return seq_type(seqs) == "protein"


def acc(true, path):
    s = score.fastsp(true, path)
    return {k: s[k] for k in ("SPFN", "SPFP", "avgErr", "LenEst", "LenRef")}


def run_pasta(w, unaligned, datatype, threads):
    out = os.path.join(w, "pasta")
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    cmd = [PASTA_PY, PASTA,
           "-i", unaligned, "-o", out, "-d", datatype, "--num-cpus", str(threads),
           "--iter-limit", "3", "--temporaries", os.path.join(w, "pasta_tmp"), "-j", "pastajob"]
    wall, cpu = timed(cmd, os.path.join(w, "pasta.log"))
    alns = [p for p in glob.glob(os.path.join(out, "pastajob.marker001.*.aln")) if "masked" not in p]
    if len(alns) != 1:
        raise RuntimeError("PASTA output not found: {}".format(alns))
    shutil.copy(alns[0], os.path.join(w, "pasta.fasta"))
    return wall, cpu


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("jobs")
    parser.add_argument("out")
    parser.add_argument("work")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--groups", type=int, default=3)
    parser.add_argument("--only", default="pasta,magus,magus-slow,self-soft,slow-soft")
    args = parser.parse_args()
    only = set(args.only.split(","))
    done = set()
    if os.path.exists(args.out):
        done = {json.loads(l)["dataset"] for l in open(args.out)}
    T = str(args.threads)
    py = [sys.executable, "-m"]

    for line in open(args.jobs):
        if not line.strip():
            continue
        name, src, k = line.split()[:3]
        if name in done:
            continue
        src = src if os.path.isabs(src) else os.path.join(REPO, src)
        w = os.path.join(args.work, name)
        os.makedirs(w, exist_ok=True)
        # restartable per method (cloud sessions stop background commands after 2 h)
        state_path = os.path.join(w, "state.json")
        row = json.load(open(state_path)) if os.path.exists(state_path) else {}
        ref = fasta.upper(fasta.read(src))
        true, unaligned = os.path.join(w, "true.fasta"), os.path.join(w, "unaligned.fasta")
        fasta.write(ref, true)
        fasta.write(fasta.ungap(ref), unaligned)
        datatype = seq_type(ref)
        row.update({"dataset": name, "threads": args.threads, "nproc": os.cpu_count(), "nseq": len(ref),
                    "protein": datatype == "protein", "datatype": datatype})

        def log_row(method, data):
            row[method] = data
            with open(state_path + ".tmp", "w") as f:
                json.dump(row, f)
            os.replace(state_path + ".tmp", state_path)  # atomic: a kill never leaves half a file
            print(json.dumps({name: {method: data}}), flush=True)

        def fresh(*paths):
            for path in paths:
                shutil.rmtree(os.path.join(w, path), ignore_errors=True)

        if "pasta" in only and "pasta" not in row:
            fresh("pasta", "pasta_tmp")
            try:
                wall, cpu = run_pasta(w, unaligned, datatype, args.threads)
                log_row("pasta", {"wall": wall, "cpu": cpu, **acc(true, os.path.join(w, "pasta.fasta"))})
            except RuntimeError as e:
                log_row("pasta", {"error": str(e)})

        for mode, extend in (("magus", "false"), ("magus-slow", "true")):
            if mode not in only and not (mode == "magus" and only & {"self-soft", "slow-soft"}):
                continue
            if mode in row and (mode != "magus" or os.path.isdir(os.path.join(w, "magus"))):
                continue
            fresh(mode)
            out = os.path.join(w, mode + ".fasta")
            wall, cpu = timed(py + ["gcmx.run_magus", "--gcmx-fastgraph", "false", "-np", T, "-d", os.path.join(w, mode),
                                    "-i", unaligned, "-o", out, "--graphbuildhmmextend", extend] + magus_flags(k),
                              os.path.join(w, mode + ".log"))
            log_row(mode, {"wall": wall, "cpu": cpu, **acc(true, out)})

        mw = os.path.join(w, "magus")
        subs = os.path.join(mw, "subalignments")
        split = os.path.join(w, "split")
        if only & {"self-soft", "slow-soft"} and ("split" not in row or not os.path.isdir(split)):
            fresh("split")
            s_wall, s_cpu = timed(py + ["gcmx.split", subs, split, str(args.groups)], os.path.join(w, "split.log"))
            log_row("split", {"wall": s_wall, "cpu": s_cpu})
        for mode in ("self-soft", "slow-soft"):
            if mode not in only or mode in row:
                continue
            bb, ev = os.path.join(w, mode + "_bb"), os.path.join(w, mode + "_ev")
            fresh(mode + "_bb", mode + "_ev", mode + "_merge")
            start = time.time()
            if mode == "self-soft":
                self_backbones(os.path.join(w, "magus.fasta"), subs, bb)
            else:
                os.makedirs(bb)
                for p in glob.glob(os.path.join(mw, "graph", "backbone_*_mafft.txt")):
                    shutil.copy(p, bb)
            prep = time.time() - start
            e_wall, e_cpu = timed(py + ["gcmx.extend", bb, unaligned, ev, "--jobs", T], os.path.join(w, mode + "_extend.log"))
            e_wall = round(e_wall + prep, 1)
            out = os.path.join(w, mode + ".fasta")
            m_wall, m_cpu = timed(py + ["gcmx.run_magus", "-np", T, "--gcmx-mclthreads", T, "-d", os.path.join(w, mode + "_merge"),
                                        "-s", split, "-b", ev, "-o", out], os.path.join(w, mode + "_merge.log"))
            s_wall, s_cpu = row["split"]["wall"], row["split"]["cpu"]
            log_row(mode, {"wall": round(row["magus"]["wall"] + e_wall + s_wall + m_wall, 1),
                           "cpu": round(row["magus"]["cpu"] + e_cpu + s_cpu + m_cpu, 1),
                           "evidence_wall": e_wall, "split_wall": s_wall, "merge_wall": m_wall, **acc(true, out)})

        with open(args.out, "a") as f:
            f.write(json.dumps(row) + "\n")
        for d in glob.glob(os.path.join(w, "*")):
            if os.path.isdir(d):
                shutil.rmtree(d, ignore_errors=True)  # keep only alignments, logs and state


if __name__ == "__main__":
    main()
