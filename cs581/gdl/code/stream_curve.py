"""Memory-light: stream gene trees from a file, accumulate ASTRID-multi / ASTRID-Pro (mean, min)
matrices, and score the FastME tree after the first n genes for each n in NGEN.
Usage: python stream_curve.py COND_DIR REP OUT.jsonl  (env NGEN=1000,10000)"""
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import methods as M  # noqa: E402
from phylo import parse_newick, read_trees, rf_error  # noqa: E402

cond, rep, out = sys.argv[1], sys.argv[2], sys.argv[3]
NG = sorted(int(x) for x in os.environ.get("NGEN", "1000,10000").split(","))
d = os.path.join(cond, rep)
true = read_trees(os.path.join(d, "s_tree.trees"))[0]
sp = sorted(true.label[v] for v in true.leaves())
idx = {s: i for i, s in enumerate(sp)}
modes = [("astrid-multi", "multi", "mean"), ("astrid-pro", "pro", "mean"), ("astrid-pro-min", "pro", "min")]
acc = {m: None for m, _, _ in modes}
t0 = time.time()
k = 0
with open(os.path.join(d, "g_true.trees")) as f:
    for line in f:
        if ";" not in line:
            continue
        g = parse_newick(line)
        rt, tg, ls = M.root_and_tag(g, M.simphy_species, idx)
        for name, mode, agg in modes:
            D, H = M.gene_distances(rt, tg, ls, len(sp), mode=mode, agg=agg)
            if acc[name] is None:
                acc[name] = [D.copy(), H.copy()]
            else:
                acc[name][0] += D; acc[name][1] += H
        k += 1
        if k in NG:
            res = {"cond": os.path.basename(cond), "rep": rep, "sqln": "true", "ngen": k, "methods": {},
                   "secs": {"total": time.time() - t0}, "stream": True}
            for name in acc:
                Dm, nm = M.average_matrix([acc[name]])
                fn, fp, i1, i2 = rf_error(M.fastme_tree(Dm, sp), true)
                res["methods"][name] = {"FN": fn, "FP": fp, "nI": i1, "RF": (fn + fp) / (i1 + i2)}
            with open(out, "a") as fo:
                fo.write(json.dumps(res) + "\n")
            print(rep, k, {n: v["FN"] for n, v in res["methods"].items()}, flush=True)
        if k >= NG[-1]:
            break
