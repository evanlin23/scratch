"""Prepare one benchmark replicate for merge experiments.

    python -m gcmx.prep TRUE_ALIGNMENT OUTDIR [--threads 4] [extra magus args]

Writes OUTDIR/{true.fasta, unaligned.fasta}, runs default MAGUS once
(OUTDIR/magus_default.fasta = baseline), and keeps the subalignments and
backbones it computed in OUTDIR/inputs/ so that merge variants can be rerun
on identical inputs. Scores and timing go to OUTDIR/prep.json.
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time

from . import fasta, score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("true_alignment")
    parser.add_argument("outdir")
    parser.add_argument("--threads", type=int, default=4)
    args, magus_args = parser.parse_known_args()

    out = os.path.abspath(args.outdir)
    os.makedirs(out, exist_ok=True)
    true = fasta.upper(fasta.read(args.true_alignment))
    fasta.write(true, os.path.join(out, "true.fasta"))
    fasta.write(fasta.ungap(true), os.path.join(out, "unaligned.fasta"))

    workdir = os.path.join(out, "magus_work")
    result = os.path.join(out, "magus_default.fasta")
    # MAGUS coordinates tasks through files in workdir/tasks; a killed run leaves
    # "running" entries behind that make a restarted run wait forever.
    shutil.rmtree(os.path.join(workdir, "tasks"), ignore_errors=True)
    start = time.time()
    with open(os.path.join(out, "magus_default.log"), "w") as log:
        subprocess.run([sys.executable, "-m", "gcmx.run_magus", "-np", str(args.threads), "-d", workdir,
                        "-i", os.path.join(out, "unaligned.fasta"), "-o", result] + magus_args,
                       stdout=log, stderr=subprocess.STDOUT, check=True,
                       cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    elapsed = time.time() - start

    for name, pattern in (("subalignments", "subalignments/subalignment_subset_*.txt"),
                          ("backbones", "graph/backbone_*_mafft.txt")):
        dest = os.path.join(out, "inputs", name)
        os.makedirs(dest, exist_ok=True)
        for path in glob.glob(os.path.join(workdir, pattern)):
            shutil.copy(path, dest)

    record = {"outdir": out, "seconds": round(elapsed, 1), "nseq": len(true)}
    record.update(score.fastsp(os.path.join(out, "true.fasta"), result))
    with open(os.path.join(out, "prep.json"), "w") as f:
        json.dump(record, f)
    print(json.dumps(record))


if __name__ == "__main__":
    main()
