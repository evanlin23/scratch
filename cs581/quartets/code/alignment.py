#!/usr/bin/env python3
"""Objective-alignment diagnostic: under each objective (raw MQSST, capminor, vote), compare the
score of the TRUE species tree with the score of the ASTRAL-IV tree on held-out replicates.
If the truth scores lower, no optimizer for that objective can recover it.
Writes results/alignment.tsv and prints a summary."""
import glob, os, subprocess, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results", "alignment.tsv")
DEV = {f"{i:02d}" for i in range(1, 6)}
from multiprocessing import Pool


def one(d):
    both = f"{d}/true_vs_astral4.tre"
    with open(both, "w") as f:
        f.write(open(f"{d}/true.tre").read().strip() + "\n" + open(f"{d}/astral4.tre").read().strip() + "\n")
    r = [d]
    for mode in ["raw", "capminor", "vote"]:
        out = subprocess.run(["/opt/runs/qtool", "score", "-m", mode, "-g", f"{d}/genes.tre", "-t", both],
                             capture_output=True, text=True).stdout.split()
        r.append(float(out[0]) - float(out[2]))
    return r


dirs = sorted(glob.glob("/opt/runs/hgt/*/*/*-est") + glob.glob("/opt/runs/hgt/*/*/*-true") + glob.glob("/opt/runs/miss/*/*/*"))
dirs = [d for d in dirs if d.rstrip("/").split("/")[-2] not in DEV and os.path.exists(f"{d}/true.tre")]
with Pool(4) as pool:
    rows = pool.map(one, dirs)
with open(OUT, "w") as f:
    f.write("dir\traw\tcapminor\tvote\n")
    for r in rows:
        f.write("\t".join(map(str, r)) + "\n")
agg = defaultdict(list)
for r in rows:
    p = r[0].split("/")
    key = ("hgt " + p[4].split("0.000001.")[1] if p[3] == "hgt" else "ils " + p[4]) + " " + p[-1]
    agg[key].append(r[1:])
print("condition\tn\t" + "\t".join(f"{m}: true<A4 / = / >" for m in ["raw", "capminor", "vote"]))
for k in sorted(agg):
    L = agg[k]
    cells = []
    for i in range(3):
        v = [x[i] for x in L]
        cells.append(f"{sum(x < -1e-9 for x in v)}/{sum(abs(x) <= 1e-9 for x in v)}/{sum(x > 1e-9 for x in v)}")
    print(f"{k}\t{len(L)}\t" + "\t".join(cells))
