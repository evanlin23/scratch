"""CPU-time summary: ProstT5 3Di prediction (4 threads) and each aligner (1 thread), summed per group.

Usage: python3 runtime.py RUNROOT P3DDIR scores.csv
"""
import csv, glob, os, sys, collections
from bbscore import read_fasta
root, p3d, sc = sys.argv[1:4]
cpu = collections.defaultdict(lambda: collections.defaultdict(float))
res = collections.defaultdict(int)
for r in csv.DictReader(open(sc)):
    g = r["group"][:4]
    cpu[g][r["method"]] += float(r["cpu"])
pt = collections.defaultdict(lambda: [0.0, 0.0])
for f in glob.glob(os.path.join(p3d, "BB*.time")):
    sid = os.path.basename(f)[:-5]; g = "RV" + sid[-5:-3]
    w, u, s = map(float, open(f).read().split()[-3:])
    pt[g][0] += w; pt[g][1] += u + s
    res[g] += sum(len(x) for x in read_fasta(os.path.join(p3d, sid + ".aa.fa")).values())
for g in sorted(pt):
    print(f"{g}: residues {res[g]}, ProstT5 wall {pt[g][0]:.0f}s cpu {pt[g][1]:.0f}s ({1000*pt[g][1]/res[g]:.2f} cpu-ms/residue)")
    for m in sorted(cpu[g]):
        print(f"   {m:14s} cpu {cpu[g][m]:8.1f}s")
