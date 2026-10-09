"""Paired comparison on the FastMulRFS data (estimated RAxML gene-family trees).
For one model condition + replicate, for each (sqln, ngen): run ASTRID-multi,
ASTRID-Pro variants, ASTRID-DISCO, ASTRAL-Pro (ours) and score the published
FastMulRFS / ASTRAL-multi / STAG / MulRF / DupTree trees on the same replicate.
Writes one JSON line per (sqln, ngen). Restartable: skips keys already in out file.
Usage: python run_fmrfs.py COND_DIR REP OUT.jsonl"""
import json
import os
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import methods as M  # noqa: E402
from phylo import parse_newick, read_trees, rf_error  # noqa: E402

cond, rep, out = sys.argv[1], sys.argv[2], sys.argv[3]
SQLN = os.environ.get("SQLN", "25,50,100,250").split(",")
NGEN = [int(x) for x in os.environ.get("NGEN", "25,100,500").split(",")]
d = os.path.join(cond, rep)
true = read_trees(os.path.join(d, "s_tree.trees"))[0]
species = sorted(true.label[v] for v in true.leaves())
idx = {s: i for i, s in enumerate(species)}
sp_of = M.simphy_species
done = set()
import glob  # noqa: E402
for fn in set(glob.glob(os.path.join(os.path.dirname(out), "fmrfs_w*.jsonl")) + [out]):
    if os.path.exists(fn):
        for l in open(fn):
            r = json.loads(l)
            done.add((r["cond"], r["rep"], r["sqln"], r["ngen"]))
PUB = {"fastmulrfs": "fastmulrfs-raxml-sqln-{s}-ngen-{n}.tree.single",
       "astral-multi(pub)": "astral-raxml-sqln-{s}-ngen-{n}.tree",
       "mulrf(pub)": "mulrf-raxml-sqln-{s}-ngen-{n}.tree",
       "stag(pub)": "stag-raxml-sqln-{s}-ngen-{n}.tree"}
cname = os.path.basename(cond.rstrip("/"))


def score(T):
    fn, fp, i1, i2 = rf_error(T, true)
    return {"FN": fn, "FP": fp, "nI": i1, "RF": (fn + fp) / (i1 + i2)}


for s in SQLN:
    allg = read_trees(os.path.join(d, "g_trees-raxml-sqlen-%s.trees" % s))
    for n in NGEN:
        if (cname, rep, s, n) in done:
            continue
        G = allg[:n]
        res = {"cond": cname, "rep": rep, "sqln": s, "ngen": n, "methods": {}, "secs": {}}
        t0 = time.time()
        tagged = [M.root_and_tag(g, sp_of, idx) for g in G]
        res["secs"]["root+tag"] = time.time() - t0
        for name, mode, agg in [("astrid-multi", "multi", "mean"), ("astrid-pro", "pro", "mean"),
                                ("astrid-pro-min", "pro", "min"), ("ortho-allnodes", "ortho_all", "mean"),
                                ("spec-allpairs", "spec_all", "mean")]:
            t0 = time.time()
            T, nm = M.species_tree_from_genes(G, sp_of, species, mode, agg, tagged=tagged)
            res["methods"][name] = dict(score(T), missing=nm)
            res["secs"][name] = time.time() - t0
        t0 = time.time()
        dts = []
        for rt, tg, ls in tagged:
            dts += M.disco_decompose(rt, tg, sp_of)
        DT = [parse_newick(x) for x in dts]
        T, nm = M.species_tree_from_genes(DT, lambda x: x, species, "multi", "mean")
        res["methods"]["astrid-disco"] = dict(score(T), missing=nm, ntrees=len(DT))
        res["secs"]["astrid-disco"] = time.time() - t0 + res["secs"]["root+tag"]
        with tempfile.TemporaryDirectory() as td:
            gf, mf, of = [os.path.join(td, x) for x in ("g.trees", "map.txt", "out.tre")]
            labs = set()
            with open(gf, "w") as f:
                for g in G:
                    f.write(g.newick() + "\n")
                    labs.update(g.label[v] for v in g.leaves())
            with open(mf, "w") as f:
                for l in sorted(labs):
                    f.write("%s %s\n" % (l, sp_of(l)))
            t0 = time.time()
            T = M.astral_pro(gf, mf, of, threads=1)
            res["secs"]["astral-pro"] = time.time() - t0
            res["methods"]["astral-pro"] = score(T)
        for name, pat in PUB.items():
            p = os.path.join(d, pat.format(s=s, n=n))
            if os.path.exists(p) and os.path.getsize(p) > 0:
                tr = read_trees(p)
                if tr:
                    res["methods"][name] = score(tr[0])
        with open(out, "a") as f:
            f.write(json.dumps(res) + "\n")
        print(cname, rep, s, n, {k: v["FN"] for k, v in res["methods"].items()}, flush=True)
