"""Consensus-alignment pilot across replicates.

    python -m gcmx.consensus_study RUNS_DIR WORK_DIR OUT.jsonl [--jobs 2] [--groups 3]

For every replicate in RUNS_DIR (fan-out results: inputs.tar.xz + prep.json),
build a GCM consensus of three independent MAGUS alignments of the same data:
the paper's published MAGUS(Fast) and MAGUS(Slow) alignments (Illinois Data
Bank Results.zip) and our own MAGUS(Fast) rerun. The consensus uses our
subsets split into `--groups` similarity groups (soft constraints) and, as a
control, our subsets as hard constraints. Every input and output is scored.
"""

import argparse
import io
import json
import os
import subprocess
import sys
import tarfile
import zipfile

from . import fasta, score
from .validate_published import URL, HttpFile

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROSE = "/opt/data/Datasets/ROSE"


def true_path(name):
    cond, rep = name.rsplit("_", 1)
    if cond.startswith("BBA"):
        return os.path.join(CODE, "..", "data", "balibase_clean", "RV100_{}.fasta".format(cond)), ("balibase", "RV100_" + cond)
    if cond == "RNASim":
        return "/opt/data/Datasets/RNASim/1000/{}/true_align.txt".format(rep), ("RNASim", rep)
    if cond.startswith("16S"):
        return "/opt/data/Datasets/Gutell/{}/{}/true_align_clean.txt".format(cond, rep), (cond, rep)
    return os.path.join(ROSE, cond, rep, "rose.aln.true.fasta"), (cond, rep)


def magus(args, workdir, output, extra):
    cmd = [sys.executable, "-m", "gcmx.run_magus", "-np", "1", "-d", workdir, "-o", output] + extra
    with open(workdir + ".log", "w") as log:
        subprocess.run(cmd, cwd=CODE, stdout=log, stderr=subprocess.STDOUT)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runs")
    parser.add_argument("work")
    parser.add_argument("out")
    parser.add_argument("--groups", type=int, default=3)
    args = parser.parse_args()

    done = set()
    if os.path.exists(args.out):
        done = {(r["dataset"], r["method"]) for r in map(json.loads, open(args.out)) if "avgErr" in r}
    z = zipfile.ZipFile(io.BufferedReader(HttpFile(URL), buffer_size=1 << 22))
    names = set(z.namelist())

    for name in sorted(os.listdir(args.runs)):
        tar = os.path.join(args.runs, name, "inputs.tar.xz")
        if not os.path.exists(tar) or (name, "consensus-soft") in done:
            continue
        src, (pub_ds, pub_rep) = true_path(name)
        w = os.path.join(args.work, name)
        os.makedirs(os.path.join(w, "cons"), exist_ok=True)
        with tarfile.open(tar) as t:
            t.extractall(w, filter="data")
        true = fasta.upper(fasta.read(src))
        fasta.write(true, os.path.join(w, "true.fasta"))
        subs, bbs = os.path.join(w, "inputs", "subalignments"), os.path.join(w, "inputs", "backbones")

        rows = []
        ours = os.path.join(w, "cons", "ours_magus_fast.fasta")
        if not os.path.exists(ours):
            magus(args, os.path.join(w, "magus_default"), ours, ["-s", subs, "-b", bbs])
        for f, label in (("gcm.txt", "pub_magus_fast"), ("gcm_slow.txt", "pub_magus_slow")):
            member = "Outputs/{}/{}/{}".format(pub_ds, pub_rep, f)
            if member in names:
                fasta.write(fasta.upper(_parse(z.read(member).decode())), os.path.join(w, "cons", label + ".fasta"))
        split_dir = os.path.join(w, "split_m{}".format(args.groups))
        if not os.path.exists(split_dir):
            subprocess.run([sys.executable, "-m", "gcmx.split", subs, split_dir, str(args.groups)], cwd=CODE, check=True,
                           stdout=subprocess.DEVNULL)
        cons_dir = os.path.join(w, "cons")
        outputs = {
            "consensus-soft": (["-s", split_dir, "-b", cons_dir], os.path.join(w, "consensus_soft.fasta")),
            "consensus-hard": (["-s", subs, "-b", cons_dir], os.path.join(w, "consensus_hard.fasta")),
        }
        for method, (extra, out) in outputs.items():
            magus(args, os.path.join(w, method), out, extra)
        for label in sorted(os.listdir(cons_dir)):
            rows.append((label[:-6], os.path.join(cons_dir, label)))
        rows += [(m, out) for m, (_, out) in outputs.items()]
        for method, path in rows:
            row = {"dataset": name, "method": method, "inputs": sorted(os.listdir(cons_dir))}
            if os.path.exists(path):
                row.update(score.fastsp(os.path.join(w, "true.fasta"), path))
            print(json.dumps(row), flush=True)
            with open(args.out, "a") as f:
                f.write(json.dumps(row) + "\n")


def _parse(text):
    seqs, name = {}, None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif line and name is not None:
            seqs[name].append(line)
    return {n: "".join(s) for n, s in seqs.items()}


if __name__ == "__main__":
    main()
