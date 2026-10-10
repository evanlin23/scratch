"""Does ASTER's reduced search (-r 1 -s 0, used in the sweeps) lose accuracy vs default search?"""
import sys, os, json
sys.path.insert(0, os.path.dirname(__file__))
import msc, stree
rows = []
for shape, n, f, k in [("cat", 32, 0.1, 200), ("cat", 32, 0.1, 500), ("bal", 32, 0.1, 200), ("cat", 64, 0.2, 100), ("bal", 64, 0.1, 300)]:
    sp = msc.caterpillar(n, f) if shape == "cat" else msc.balanced(n, f)
    true = msc.species_newick(sp)
    same = 0; fr = fd = 0
    for rep in range(12):
        g = msc.sim_genes(sp, k, 99000 + rep + n)
        a = stree.astral(g, ("-r", "1", "-s", "0")); b = stree.astral(g)
        fr += stree.fn_fp(true, a)[0]; fd += stree.fn_fp(true, b)[0]; same += stree.fn_fp(a, b)[0] == 0
    rows.append(dict(shape=shape, n=n, f=f, k=k, reps=12, FN_reduced=fr / 12, FN_default=fd / 12, identical=same / 12))
    print(rows[-1], flush=True)
json.dump(rows, open(os.path.join(os.path.dirname(__file__), "..", "results", "astral_search_check.json"), "w"), indent=1)
