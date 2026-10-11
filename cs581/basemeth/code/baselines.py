"""Whole-dataset baselines (no MAGUS): FAMSA2, MUSCLE5, regressive T-Coffee, and PASTA 1.8.3 aligner variants.

    PYTHONPATH=cs581/code python3 cs581/basemeth/code/baselines.py JOBS WORK OUT.jsonl --methods famsa,muscle5,regressive
    ... --methods pasta-mafft,pasta-probcons,pasta-prank,pasta-muscle [--iter 3]

Uses WORK/NAME/{true.fasta, unaligned input} from basemeth.py prep. 4 threads, --cap seconds wall limit
(a run that hits it is recorded as a timeout). Restartable: finished (dataset, method) rows are skipped.
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from basemeth import ALN2, BIO, accuracy, jobs, timed  # noqa: E402
from gcmx import fasta  # noqa: E402
from gcmx.e2e_bench import PASTA, PASTA_PY, seq_type  # noqa: E402

ENV = dict(os.environ, PATH=ALN2 + ":" + BIO + ":" + os.environ["PATH"])


def command(method, inp, out, nseq, datatype, d, iters):
    if method == "famsa":
        return [BIO + "/famsa", "-t", "4", inp, out]
    if method == "muscle5":  # -align is recommended up to ~1,000 sequences, Super5 above
        return [ALN2 + "/muscle", "-align" if nseq <= 1000 else "-super5", inp, "-output", out, "-threads", "4"]
    if method == "regressive":  # mBed guide trees hang in this T-Coffee 12 build; NJ tree, Clustal Omega children
        return [ALN2 + "/t_coffee", "-reg", "-seq", inp, "-reg_nseq", "100", "-reg_tree", "nj",
                "-reg_method", "clustalo_msa", "-outfile", out, "-thread", "4"]
    if method.startswith("pasta-"):
        return [PASTA_PY, PASTA, "-i", inp, "-o", os.path.join(d, "pasta_out"), "-d", datatype, "--num-cpus", "4",
                "--iter-limit", str(iters), "--aligner", method[6:], "--temporaries", os.path.join(d, "tmp"),
                "-j", "job"]
    raise ValueError(method)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("jobs")
    p.add_argument("work")
    p.add_argument("out")
    p.add_argument("--methods", default="famsa,muscle5,regressive")
    p.add_argument("--iter", type=int, default=3, help="PASTA iterations")
    p.add_argument("--cap", type=int, default=1800)
    p.add_argument("--only", default="", help="comma list of dataset names")
    args = p.parse_args()
    done = set()
    if os.path.exists(args.out):
        done = {(r["dataset"], r["method"]) for r in map(json.loads, open(args.out))}
    only = set(filter(None, args.only.split(",")))
    for j in jobs(args.jobs):
        if only and j["name"] not in only:
            continue
        w = os.path.join(args.work, j["name"])
        if not os.path.exists(os.path.join(w, "true.fasta")):
            continue
        ref = fasta.upper(fasta.read(os.path.join(w, "true.fasta")))
        unaligned = j["unaligned"] or os.path.join(w, "unaligned.fasta")
        if not os.path.exists(unaligned):
            fasta.write(fasta.ungap(ref), unaligned)
        nseq = len(fasta.read(unaligned))
        t = seq_type(ref)
        datatype = {"protein": "Protein", "dna": "DNA", "rna": "RNA"}[t]
        for m in args.methods.split(","):
            key = m + ("" if not m.startswith("pasta-") else "-it{}".format(args.iter))
            if (j["name"], key) in done:
                continue
            d = os.path.join(w, "base_" + key)
            shutil.rmtree(d, ignore_errors=True)
            os.makedirs(d)
            out = os.path.join(d, "out.fa")
            cmd = command(m, unaligned, out, nseq, datatype, d, args.iter)
            if args.cap:
                cmd = ["timeout", str(args.cap)] + cmd
            row = {"dataset": j["name"], "method": key, "nseq": nseq}
            try:
                old = os.environ.copy()
                os.environ.update(ENV)
                try:
                    row["wall"], row["cpu"] = timed(cmd, os.path.join(d, "log.txt"), cwd=d)
                finally:
                    os.environ.clear()
                    os.environ.update(old)
                if m.startswith("pasta-"):
                    alns = [x for x in glob.glob(os.path.join(d, "pasta_out", "job.marker001.*.aln"))
                            if "masked" not in x]
                    out = alns[0]
                row.update(accuracy(ref, out))
            except Exception as e:
                row["error"] = repr(e)[:200]
            with open(args.out, "a") as f:
                f.write(json.dumps(row) + "\n")
            print(json.dumps(row), flush=True)
            shutil.rmtree(os.path.join(d, "tmp"), ignore_errors=True)


if __name__ == "__main__":
    main()
