"""Paired comparison on the DISCO data (doi:10.13012/B2IDB-4050038_V1): estimated
gene-family trees g_<sqln>.trees (with branch support) or g_true.trees ('true').
Same methods as run_fmrfs.py plus support-weighted variants. Writes one JSON line per
(sqln, ngen). Usage: python run_disco.py COND_DIR REP OUT.jsonl"""
import json
import os
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import methods as M  # noqa: E402
from phylo import parse_newick, read_trees, rf_error  # noqa: E402

cond, rep, out = sys.argv[1], sys.argv[2], sys.argv[3]
SQLN = os.environ.get("SQLN", "100").split(",")
NGEN = [int(x) for x in os.environ.get("NGEN", "100,1000").split(",")]
ASTRAL = os.environ.get("ASTRAL", "1") == "1"
d = os.path.join(cond, rep)
true = read_trees(os.path.join(d, "s_tree.trees"))[0]
species = sorted(true.label[v] for v in true.leaves())
idx = {s: i for i, s in enumerate(species)}
sp_of = M.simphy_species
done = set()
import glob  # noqa: E402
for fn in set(glob.glob(os.path.join(os.path.dirname(out), "disco_w*.jsonl")) + [out]):
    if os.path.exists(fn):
        for l in open(fn):
            r = json.loads(l)
            done.add((r["cond"], r["rep"], r["sqln"], r["ngen"]))
PUB = {}
cname = os.path.basename(cond.rstrip("/"))


def score(T):
    fn, fp, i1, i2 = rf_error(T, true)
    return {"FN": fn, "FP": fp, "nI": i1, "RF": (fn + fp) / (i1 + i2)}


for s in SQLN:
    allg = read_trees(os.path.join(d, "g_%s.trees" % s))
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
                                ("spec-allpairs", "spec_all", "mean"),
                                ("astrid-multi-w", "multi", "mean-w"), ("astrid-pro-w", "pro", "mean-w")]:
            if agg.endswith("-w") and s == "true":
                continue
            t0 = time.time()
            T, nm = M.species_tree_from_genes(G, sp_of, species, mode, agg.replace("-w", ""), tagged=tagged,
                                              weighted=agg.endswith("-w"))
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
          if ASTRAL:
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
            T = M.astral_pro(gf, mf, of, threads=int(os.environ.get("ATHREADS", "1")))
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
