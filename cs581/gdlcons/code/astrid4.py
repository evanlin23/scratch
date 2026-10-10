"""Limiting ASTRID-multi four-point comparison on a 4-taxon pool, with block SEs.
ASTRID-multi = per-family mean internode distance over copy pairs, averaged over families
(all families contain all 4 species here).  Correct split AB|CD needs
S_AB = d(A,B) + d(C,D) to be the smallest of the three sums.
usage: python astrid4.py POOL NFAM BLOCK"""
import os, sys, statistics as st
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core, methods as M
pool, n, B = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
fams = [l for l in open(pool) if ';' in l][:n]
spi = {s: i for i, s in enumerate("ABCD")}
per = []
for f in fams:
    rt, _ = core.true_rt(f, spi)
    ls = [spi[core.sp_of(rt.label[v])] if not rt.children[v] else -1 for v in range(len(rt.parent))]
    m, h = M.gene_distances(rt, [False] * len(rt.parent), ls, 4, mode="multi")
    per.append(m)
per = np.array(per)
def sums(D):
    return D[0, 1] + D[2, 3], D[0, 2] + D[1, 3], D[0, 3] + D[1, 2]
blocks = [sums(per[i:i + B].mean(axis=0)) for i in range(0, len(per) - B + 1, B)]
for name, j in (("AC|BD", 1), ("AD|BC", 2)):
    d = [b[j] - b[0] for b in blocks]   # positive = correct split shorter = good
    print("%s - AB|CD four-point sum: %.4f +- %.4f (blocks=%d)" % (name, st.mean(d), st.stdev(d) / len(d) ** 0.5, len(d)))
print("overall sums AB|CD, AC|BD, AD|BC:", [round(x, 4) for x in sums(per.mean(axis=0))])
