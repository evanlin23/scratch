#!/usr/bin/env python3
"""Drop-in replacement for PASTA's MUSCLE merger.

PASTA calls `<muscle path> -in1 A -in2 B -out OUT -quiet -profile`. Point the
[muscle] path of a PASTA config at a shell shim that runs this script and use
`--merger muscle`; then every pairwise merge in PASTA goes through MERGER
(env PAIRMERGE_METHOD: progdp | gcm | opal | muscle3 | mafft-merge).

progdp/gcm need backbone evidence. Like MAGUS, we draw PAIRMERGE_R random
backbones of up to PAIRMERGE_M sequences (half from each side), align them
with MAFFT (PAIRMERGE_BBALIGNER: linsi | auto), and weight A-B column pairs by
how many backbones align their letters.
"""

import os
import random
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mergers  # noqa: E402
from gcmx import fasta  # noqa: E402


def backbones(a, b, outdir, r, m, aligner, seed):
    A, B = fasta.ungap(fasta.read(a)), fasta.ungap(fasta.read(b))
    rng = random.Random(seed)
    os.makedirs(outdir)
    for i in range(r):
        half = m // 2
        na, nb = min(half, len(A)), min(half, len(B))
        na, nb = min(len(A), m - nb), min(len(B), m - na)  # use the slack of a small side
        taxa = rng.sample(sorted(A), na) + rng.sample(sorted(B), nb)
        seqs = {t: (A.get(t) or B[t]) for t in taxa}
        unal = os.path.join(outdir, "..", "bb_{}.fa".format(i))
        fasta.write(seqs, unal)
        cmd = ["mafft", "--quiet", "--thread", "1"]
        cmd += ["--localpair", "--maxiterate", "1000"] if aligner == "linsi" else ["--auto"]
        with open(os.path.join(outdir, "backbone_{}_mafft.txt".format(i + 1)), "w") as o:
            subprocess.run(cmd + [unal], stdout=o, stderr=subprocess.DEVNULL, check=True)


def main(argv):
    a, b, out = argv[argv.index("-in1") + 1], argv[argv.index("-in2") + 1], argv[argv.index("-out") + 1]
    method = os.environ.get("PAIRMERGE_METHOD", "progdp")
    work = tempfile.mkdtemp(prefix="pm_", dir=os.environ.get("PAIRMERGE_TMP"))
    try:
        if method in ("progdp", "gcm"):
            bb = os.path.join(work, "bb")
            backbones(a, b, bb, int(os.environ.get("PAIRMERGE_R", "10")), int(os.environ.get("PAIRMERGE_M", "200")),
                      os.environ.get("PAIRMERGE_BBALIGNER", "linsi"), seed=len(a) + len(b))
            mergers.gcm(a, b, out, work, 1, backbones=bb, trace="progdp" if method == "progdp" else "minclusters")
        elif method == "opal":
            mergers.opal(a, b, out, work, 1)
        elif method == "muscle3":
            mergers.muscle3(a, b, out, work, 1)
        elif method == "mafft-merge":
            mergers.mafft_merge(a, b, out, work, 1)
        else:
            raise SystemExit("unknown PAIRMERGE_METHOD " + method)
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main(sys.argv[1:])
