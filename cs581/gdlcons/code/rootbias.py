"""Systematic misrooting (root on a leaf edge of one fixed species) on small trees.
usage: python rootbias.py OUT.jsonl NFAM"""
import json, os, random, sys, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import apply_error, encode, quartet_scores, simulate, true_rt
from phylo import parse_newick
SH = {"cat4": "(((A:1,B:1)x:{T},C:{c})y:{T},D:{d});",
      "cat5": "((((A:1,B:1)x:{T},C:{c})y:{T},D:{d})z:{T},E:{e});",
      "bal4": "((A:1,B:1)x:{T},(C:1,D:1)y:{T});"}
out, nfam = sys.argv[1], int(sys.argv[2])
for shape, T, inner in itertools.product(["cat4", "cat5", "bal4"], [0.05, 0.3], [(3.0, 0.5), (1.0, 1.0), (3.0, 3.0)]):
    nwk = SH[shape].format(T=T, c=1 + T, d=1 + 2 * T, e=1 + 3 * T)
    st = parse_newick(nwk)
    lam, mu = [0.0] * len(st.parent), [0.0] * len(st.parent)
    for v in range(len(st.parent)):
        if st.parent[v] >= 0:
            lam[v], mu[v] = (1.0, 1.0) if not st.children[v] else inner
    fams, tries, over = simulate(st, lam, mu, nfam, 5, min_species=4, cap=4000, max_over=nfam, max_tries=200 * nfam)
    species = sorted(st.label[v] for v in st.leaves()); spi = {s: i for i, s in enumerate(species)}
    base = [true_rt(f, spi) for f in fams]
    rng = random.Random(3)
    res = {}
    for sp in species:
        err = [apply_error("rsp-" + sp, 1.0, rt, tg, spi, rng) for rt, tg in base]
        sc = quartet_scores([encode(r, t) for r, t in err], [r for r, _ in err], True, st.newick())
        per = {}
        for node, t, _, s in sc: per.setdefault(node, {})[t] = s
        res[sp] = [[d["t1"], d["t2"], d["t3"]] for _, d in sorted(per.items())]
    r = {"shape": shape, "T": T, "inner": inner, "nfam": nfam, "scores": res}
    open(out, "a").write(json.dumps(r) + "\n")
    print(shape, T, inner, {k: [round((a - max(b, c)) / max(a + b + c, 1e-9), 3) for a, b, c in v] for k, v in res.items()}, flush=True)
