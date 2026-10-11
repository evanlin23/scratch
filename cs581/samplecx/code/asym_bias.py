"""Asymptotic bias of ASTRID/NJst under iid taxon deletion (Rhodes, Nute & Warnow 2020).

For full MSC gene trees, E[observed internode count | i,j kept] = sum_v (1 - p^K_v) exactly.
Averaging this over many simulated gene trees gives (up to MC error shared with the full-data
matrix) the limit matrix of plain ASTRID under deletion. We compare the tree FastME/NJ build
from the limit matrix with the tree from the complete-data matrix of the SAME gene trees, so
MC noise is paired out. A difference = inconsistency of plain ASTRID at that (tree, p).
The corrected (HT) estimator has the complete-data matrix as its limit by construction.
Output: results/asym_bias.csv
"""
import sys, os, csv, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
import msc, stree
from msc import Node

def yule(n, mean_bl, rng):
    nodes = [Node(label=f"t{i}") for i in range(n)]
    while len(nodes) > 1:
        i, j = rng.choice(len(nodes), 2, replace=False)
        a, b = nodes[i], nodes[j]
        for x in (a, b):
            if x.children:
                x.length = rng.exponential(mean_bl)
        nodes = [x for k, x in enumerate(nodes) if k not in (i, j)] + [Node([a, b])]
    nodes[0].length = np.inf
    for nd in [nodes[0]]:
        pass
    return nodes[0]

def job(args):
    kind, n, mbl, seed, K = args
    rng = np.random.default_rng(seed)
    sp = yule(n, mbl, rng) if kind == "yule" else (msc.caterpillar(n, mbl) if kind == "cat" else msc.balanced(n, mbl))
    true = msc.species_newick(sp)
    taxa = [f"t{i}" for i in range(n)]
    g = msc.sim_genes(sp, K, seed)
    M0, _ = stree.distance_matrix(g, taxa)
    t_full = stree.fastme(M0, taxa); nj_full = stree.nj(M0, taxa)
    rows = []
    for p in [0.3, 0.5, 0.7, 0.9]:
        Mp, _ = stree.distance_matrix(g, taxa, mode="expected", p=p)
        t_p = stree.fastme(Mp, taxa); nj_p = stree.nj(Mp, taxa)
        rows.append(dict(kind=kind, n=n, mbl=mbl, seed=seed, p=p,
                         fn_full_astrid=stree.fn_fp(true, t_full)[0], fn_lim_astrid=stree.fn_fp(true, t_p)[0],
                         diff_astrid=stree.fn_fp(t_full, t_p)[0],
                         fn_full_nj=stree.fn_fp(true, nj_full)[0], fn_lim_nj=stree.fn_fp(true, nj_p)[0],
                         diff_nj=stree.fn_fp(nj_full, nj_p)[0]))
    return rows

if __name__ == "__main__":
    K = int(sys.argv[1]); nproc = int(sys.argv[2])
    jobs = [("yule", n, mbl, s, K) for n in [10, 20, 30] for mbl in [0.05, 0.2, 1.0] for s in range(10)]
    jobs += [(k, n, f, 0, K) for k in ["cat", "bal"] for n in [8, 16, 32] for f in [0.05, 0.2, 1.0]]
    out = os.path.join(os.path.dirname(__file__), "..", "results", "asym_bias.csv")
    with Pool(nproc) as pool:
        rows = [r for rs in pool.map(job, jobs) for r in rs]
    with open(out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    import pandas as pd
    d = pd.DataFrame(rows)
    print(d.groupby(["kind", "mbl", "p"])[["fn_full_astrid", "fn_lim_astrid", "diff_astrid", "diff_nj"]].mean().round(3).to_string())
