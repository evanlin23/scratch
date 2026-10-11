# AI-assisted (Claude), exploration code for CS581 project
"""Noise floor of FastTree nRF: the same MAGUS alignment with 8 random columns removed (about 0.1% of the
columns; carries essentially the same information). Also re-runs MAGUS unchanged on 3 draws to check that
trees.py reproduces the logged trees.

    python3 noise.py > data/noise_jobs.tsv
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from measures import read_fasta  # noqa: E402

W = "/opt/work/treecrit"
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")
paths = {(r["key"], r["method"]): r["path"] for r in map(json.loads, open(os.path.join(D, "aln_paths.jsonl")))}
os.makedirs(os.path.join(W, "noise"), exist_ok=True)
keys = sorted({k for k, m in paths if m == "magus" and k.startswith("SIMHIGH")})
for i, k in enumerate(keys):
    src = paths[(k, "magus")]
    aln = read_fasta(src)
    L = len(next(iter(aln.values())))
    rng = random.Random(1000 + i)
    drop = set(rng.sample(range(L), 8))
    dst = os.path.join(W, "noise", k + ".drop8.fasta")
    with open(dst, "w") as f:
        for t, s in aln.items():
            f.write(">{}\n{}\n".format(t, "".join(c for j, c in enumerate(s) if j not in drop)))
    print("{}\tmagus_drop8\t{}".format(k, dst))
    if i < 3:
        print("{}\tmagus_rerun\t{}".format(k, src))
