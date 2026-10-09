"""Consistency probe on a small species tree: simulate many pure-GDL families with
per-branch rates and report, for each distance mode, the limiting quartet sums
(four-point condition) and the FastME tree; plus ASTRAL-Pro on the same families.

Rates are given per species-tree node label (branch above the node), e.g.
  --tree "(((A:1,B:1)x:1,C:2)y:1,D:3)r;" --rates "x=3,0;A=0,3;B=0,3;C=0,3"
(lam,mu); unspecified branches use --default lam,mu.
"""
import argparse
import itertools
import json
import os
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import methods as M  # noqa: E402
from gdlsim import simulate  # noqa: E402
from phylo import parse_newick, rf_error  # noqa: E402

MODES = [("multi", "mean", False), ("pro", "mean", True), ("pro", "mean", False),
         ("pro", "min", False), ("ortho_all", "mean", False), ("spec_all", "mean", False)]


def quartet_sums(D, idx, q):
    a, b, c, d = [idx[x] for x in q]
    return {"%s%s|%s%s" % (q[0], q[1], q[2], q[3]): D[a, b] + D[c, d],
            "%s%s|%s%s" % (q[0], q[2], q[1], q[3]): D[a, c] + D[b, d],
            "%s%s|%s%s" % (q[0], q[3], q[1], q[2]): D[a, d] + D[b, c]}


def run(tree, rates, default, nfam, seed, min_species, astral=True, quartets=None, root_len=0.0, modes=None):
    st = parse_newick(tree)
    lam = [default[0]] * len(st.parent)
    mu = [default[1]] * len(st.parent)
    for v in range(len(st.parent)):
        if callable(rates):
            lam[v], mu[v] = rates(st, v)
        elif st.label[v] in rates:
            lam[v], mu[v] = rates[st.label[v]]
    species = sorted(st.label[v] for v in st.leaves())
    idx = {s: i for i, s in enumerate(species)}
    t0 = time.time()
    fams, tries, over = simulate(st, lam, mu, nfam, seed, min_species=min_species, root_len=root_len)
    simt = time.time() - t0
    G = [parse_newick(f) for f in fams]
    res = {"nfam": nfam, "tries": tries, "overflow": over, "sim_s": round(simt, 1),
           "mean_leaves": float(np.mean([len(g.leaves()) for g in G])), "methods": {}}
    sp_of = M.simphy_species
    tagged_inf = [M.root_and_tag(g, sp_of, idx) for g in G]
    tagged_true = [M.root_and_tag(g, sp_of, idx, truetags=True) for g in G]
    for mode, agg, tt in (modes or MODES):
        tg = tagged_true if tt else tagged_inf
        per = [M.gene_distances(rt, t, ls, len(species), mode=mode, agg=agg) for rt, t, ls in tg]
        D, nmiss = M.average_matrix(per)
        T = M.fastme_tree(D, species)
        fn = rf_error(T, st)[0]
        name = "%s-%s%s" % (mode, agg, "-truetags" if tt else "")
        r = {"FN": fn, "missing": nmiss}
        if quartets:
            r["quartets"] = {"".join(q): {k: round(v, 4) for k, v in quartet_sums(D, idx, q).items()}
                             for q in quartets}
        res["methods"][name] = r
    # DISCO + ASTRID
    dts = []
    for rt, t, ls in tagged_inf:
        dts += M.disco_decompose(rt, t, sp_of)
    if dts:
        DT = [parse_newick(x) for x in dts]
        T, nm = M.species_tree_from_genes(DT, lambda x: x, species, "multi", "mean")
        res["methods"]["astrid-disco"] = {"FN": rf_error(T, st)[0], "ntrees": len(DT)}
    if astral:
        with tempfile.TemporaryDirectory() as td:
            gf, mf, of = [os.path.join(td, x) for x in ("g.trees", "map.txt", "out.tre")]
            with open(gf, "w") as f:
                f.write("\n".join(fams) + "\n")
            labs = set()
            for g in G:
                labs.update(g.label[v] for v in g.leaves())
            with open(mf, "w") as f:
                for l in sorted(labs):
                    f.write("%s %s\n" % (l, sp_of(l)))
            t0 = time.time()
            T = M.astral_pro(gf, mf, of, threads=1)
            res["methods"]["astral-pro"] = {"FN": rf_error(T, st)[0], "secs": round(time.time() - t0, 1)}
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--rates", default="")
    ap.add_argument("--default", default="0,0")
    ap.add_argument("--nfam", type=int, nargs="+", default=[1000])
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--min-species", type=int, default=2)
    ap.add_argument("--quartets", default="")
    ap.add_argument("--root-len", type=float, default=0.0)
    ap.add_argument("--no-astral", action="store_true")
    a = ap.parse_args()
    rates = {}
    for kv in filter(None, a.rates.split(";")):
        k, v = kv.split("=")
        rates[k] = tuple(float(x) for x in v.split(","))
    default = tuple(float(x) for x in a.default.split(","))
    qs = [tuple(q) for q in a.quartets.split(",")] if a.quartets else None
    for n in a.nfam:
        r = run(a.tree, rates, default, n, a.seed, a.min_species, astral=not a.no_astral,
                quartets=qs, root_len=a.root_len)
        print(json.dumps(r))
