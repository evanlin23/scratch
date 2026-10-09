"""Subset-aligner bake-off: how accurately do different aligners align MAGUS's subsets?

    python -m gcmx.subset_bakeoff TRUE.fa SUBALIGNMENT_DIR OUTDIR [--methods a,b,...] [--jobs 4]

For every subset (taken from MAGUS's subset alignments, gaps removed) and every
method, aligns the subset, scores it against the true alignment restricted to
the subset, and writes OUTDIR/<method>/<subset file> (usable directly as MAGUS
`-s` input) plus OUTDIR/results.jsonl. Errors are pooled over subsets by
homology counts (so large subsets weigh more), as SP error would be.
"""

import argparse
import concurrent.futures
import json
import os
import shutil
import subprocess
import tempfile
import time

from . import fasta, score

MAGUS_MAFFT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "MAGUS", "magus", "tools", "mafft", "mafft")
BIO = "/opt/mm/root/envs/bio/bin"
ALN2 = "/opt/mm/root/envs/aln2/bin"

METHODS = {
    # MAGUS's own subset aligner (exact command)
    "linsi": [MAGUS_MAFFT, "--localpair", "--maxiterate", "1000", "--ep", "0.123", "--quiet", "--anysymbol", "{in}"],
    "ginsi": [MAGUS_MAFFT, "--globalpair", "--maxiterate", "1000", "--quiet", "--anysymbol", "{in}"],
    "einsi": [MAGUS_MAFFT, "--genafpair", "--ep", "0", "--maxiterate", "1000", "--quiet", "--anysymbol", "{in}"],
    "fftnsi": [MAGUS_MAFFT, "--retree", "2", "--maxiterate", "1000", "--quiet", "--anysymbol", "{in}"],
    "muscle5": [BIO + "/muscle", "-align", "{in}", "-output", "{out}"],
    "prank": [ALN2 + "/prank", "-d={in}", "-o={outbase}", "-f=fasta", "-quiet"],
    "prank+F": [ALN2 + "/prank", "-d={in}", "-o={outbase}", "-f=fasta", "-quiet", "+F"],
    "kalign": [ALN2 + "/kalign", "-i", "{in}", "-o", "{out}"],
    "clustalo": ["clustalo", "-i", "{in}", "-o", "{out}", "--force"],
    "famsa": [BIO + "/famsa", "{in}", "{out}"],
    "tcoffee": [ALN2 + "/t_coffee", "{in}", "-output", "fasta_aln", "-outfile", "{out}", "-quiet", "-n_core", "1"],
}


def run_method(method, unaligned, out_path):
    with tempfile.TemporaryDirectory() as tmp:
        inp = os.path.join(tmp, "in.fa")
        shutil.copy(unaligned, inp)
        outbase = os.path.join(tmp, "out")
        out = outbase + ".fa"
        cmd = [a.format(**{"in": inp, "out": out, "outbase": outbase}) for a in METHODS[method]]
        start = time.time()
        if "{out}" in " ".join(METHODS[method]) or "{outbase}" in " ".join(METHODS[method]):
            subprocess.run(cmd, cwd=tmp, check=True, capture_output=True)
        else:
            with open(out, "w") as f:
                subprocess.run(cmd, cwd=tmp, check=True, stdout=f, stderr=subprocess.DEVNULL)
        elapsed = time.time() - start
        if method.startswith("prank"):
            out = outbase + ".best.fas"
        aln = fasta.upper(fasta.read(out))
        fasta.write(aln, out_path)
    return elapsed


def job(args, method, filename, true):
    sub = fasta.read(os.path.join(args.subalignments, filename))
    with tempfile.TemporaryDirectory() as tmp:
        unaligned = os.path.join(tmp, "unaligned.fa")
        fasta.write(fasta.upper(fasta.ungap(sub)), unaligned)
        ref = os.path.join(tmp, "ref.fa")
        fasta.write(fasta.restrict(true, list(sub)), ref)
        out_path = os.path.join(args.outdir, method, filename)
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        row = {"method": method, "subset": filename, "nseq": len(sub)}
        try:
            row["seconds"] = round(run_method(method, unaligned, out_path), 2)
            row.update(score.fastsp(ref, out_path))
        except Exception as e:
            row["error"] = repr(e)[:200]
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("true")
    parser.add_argument("subalignments")
    parser.add_argument("outdir")
    parser.add_argument("--methods", default=",".join(METHODS))
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()

    true = fasta.upper(fasta.read(args.true))
    files = sorted(os.listdir(args.subalignments))
    methods = args.methods.split(",")
    os.makedirs(args.outdir, exist_ok=True)
    rows = []
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        futures = [pool.submit(job, args, m, f, true) for m in methods for f in files]
        for fut in concurrent.futures.as_completed(futures):
            row = fut.result()
            rows.append(row)
            with open(os.path.join(args.outdir, "results.jsonl"), "a") as out:
                out.write(json.dumps(row) + "\n")

    print("| method | subsets | SPFN | SPFP | avg | total seconds |")
    print("|---|---|---|---|---|---|")
    for m in methods:
        ok = [r for r in rows if r["method"] == m and "shared" in r]
        if not ok:
            print("| {} | 0 | failed | | | |".format(m))
            continue
        shared, ref, est = (sum(r[k] for r in ok) for k in ("shared", "refHom", "estHom"))
        spfn, spfp = 1 - shared / ref, 1 - shared / est
        print("| {} | {} | {:.4f} | {:.4f} | {:.4f} | {:.0f} |".format(
            m, len(ok), spfn, spfp, (spfn + spfp) / 2, sum(r["seconds"] for r in ok)))


if __name__ == "__main__":
    main()
