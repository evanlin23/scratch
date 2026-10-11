"""Per-set mean pairwise identity of the MAFFT L-INS-i alignment (no reference used) vs. paired SP gains.

Usage: python3 identity.py RUNROOT scores.csv out.tsv
Also evaluates a reference-free switch: use linsi3di when estimated identity < t, else linsi.
"""
import csv, glob, os, sys, collections
from itertools import combinations
from bbscore import read_fasta

root, sc, out = sys.argv[1:4]
S = collections.defaultdict(dict)
for r in csv.DictReader(open(sc)):
    S[r["set"]][r["method"]] = float(r["SP"])
ident = {}
for d in sorted(glob.glob(os.path.join(root, "BB*"))):
    a = list(read_fasta(os.path.join(d, "linsi.fa")).values())
    ids = []
    for x, y in combinations(a, 2):
        p = [(c, e) for c, e in zip(x.upper(), y.upper()) if c != "-" and e != "-"]
        if p:
            ids.append(sum(c == e for c, e in p) / len(p))
    ident[os.path.basename(d)] = sum(ids) / len(ids)
with open(out, "w") as f:
    f.write("set\tidentity\tlinsi\tlinsi3di\tmc2\tmc3di\n")
    for s in sorted(ident):
        f.write(f"{s}\t{ident[s]:.3f}\t" + "\t".join(f"{S[s].get(m, float('nan')):.4f}" for m in ("linsi", "linsi3di", "mc2", "mc3di")) + "\n")
sets = sorted(ident)
n = len(sets)
base = sum(S[s]["linsi"] for s in sets) / n
for t in (0.15, 0.2, 0.25, 0.3, 0.35):
    sw = sum(S[s]["linsi3di"] if ident[s] < t else S[s]["linsi"] for s in sets) / n
    sw2 = sum(S[s]["mc3di"] if ident[s] < t else S[s]["linsi"] for s in sets) / n
    k = sum(ident[s] < t for s in sets)
    print(f"t={t:.2f}: {k} sets switched; mean SP linsi {base:.4f}  switch->linsi3di {sw:.4f}  switch->mc3di {sw2:.4f}")
for lo, hi in ((0, .15), (.15, .2), (.2, .25), (.25, .3), (.3, .4), (.4, 1)):
    ss = [s for s in sets if lo <= ident[s] < hi]
    if ss:
        d = [S[s]["linsi3di"] - S[s]["linsi"] for s in ss]; d2 = [S[s]["mc2"] - S[s]["linsi"] for s in ss]
        print(f"identity [{lo:.2f},{hi:.2f}): n={len(ss):3d}  dSP linsi3di {sum(d)/len(d):+.3f} (W {sum(x>0 for x in d)}/L {sum(x<0 for x in d)})  dSP mc2 {sum(d2)/len(d2):+.3f}")
