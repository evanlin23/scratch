"""High-precision Monte Carlo of the *limiting* ASTRID-multi (and ideal ASTRID-Pro) matrices
under pure GDL. Simulates families in batches. For each batch it writes the per-species-pair
sums of per-gene averages and the gene counts, so batches from many processes/seeds can be
pooled and a batch-means standard error computed for every four-point margin.
No tree building. Restartable: appends one JSON line per batch.

Usage: python hiprec.py TREE RATES_JSON BATCH NBATCH SEED OUT.jsonl [modes]
modes: comma list from multi,pro_true_root (default both)"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import methods as M  # noqa: E402
from gdlsim import simulate  # noqa: E402
from phylo import parse_newick  # noqa: E402

tree, rates_s, batch, nbatch, seed, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), \
    int(sys.argv[5]), sys.argv[6]
modes = (sys.argv[7] if len(sys.argv) > 7 else "multi,pro_true_root").split(",")
st = parse_newick(tree)
rates = json.loads(rates_s)
lam = [0.0] * len(st.parent)
mu = [0.0] * len(st.parent)
for v in range(len(st.parent)):
    if st.label[v] in rates:
        lam[v], mu[v] = rates[st.label[v]]
species = sorted(st.label[v] for v in st.leaves())
idx = {s: i for i, s in enumerate(species)}
n = len(species)
start = 0
if os.path.exists(out):
    start = sum(1 for _ in open(out))
for b in range(start, nbatch):
    fams, tries, over = simulate(st, lam, mu, batch, seed * 1000003 + b, min_species=2)
    acc = {m: [np.zeros((n, n)), np.zeros((n, n))] for m in modes}
    for f in fams:
        g = parse_newick(f)
        rt, tg, ls = M.root_and_tag(g, M.simphy_species, idx, truetags=True)
        for m in modes:
            if m == "multi":
                D, H = M.gene_distances(rt, tg, ls, n, mode="multi")
            else:
                D, H = M.gene_distances(rt, tg, ls, n, mode="pro", count_root=True)
            acc[m][0] += D
            acc[m][1] += H
    rec = {"batch": b, "seed": seed, "nfam": len(fams), "tries": tries, "overflow": over, "species": species,
           "acc": {m: [a[0].tolist(), a[1].tolist()] for m, a in acc.items()}}
    with open(out, "a") as f:
        f.write(json.dumps(rec) + "\n")
