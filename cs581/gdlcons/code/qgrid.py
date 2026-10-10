"""Small-tree study: ASTRAL-Pro quartet-support margins as a function of the tagging
error rate q (D->S relabelling and symmetric flips) and the rerooting rate p.

Species trees: 4-taxon caterpillar / balanced and 5-taxon caterpillar with internal
branch length T; rates: internal (lam_i, mu_i), terminal (lam_t, mu_t).
usage: python qgrid.py OUT.jsonl NFAM
"""
import itertools
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import apply_error, encode, plain, quartet_scores, simulate, true_rt  # noqa: E402
from phylo import parse_newick  # noqa: E402

SHAPES = {
    "cat4": "(((A:1,B:1)x:{T},C:{c})y:{T},D:{d});",
    "bal4": "((A:1,B:1)x:{T},(C:1,D:1)y:{T});",
    "cat5": "((((A:1,B:1)x:{T},C:{c})y:{T},D:{d})z:{T},E:{e});",
}
QS = [0, 0.1, 0.25, 0.5, 0.75, 1.0]
PS = [0.1, 0.3, 1.0]


def run(shape, T, inner, term, nfam, seed):
    nwk = SHAPES[shape].format(T=T, c=1 + T, d=1 + 2 * T, e=1 + 3 * T)
    st = parse_newick(nwk)
    lam, mu = [0.0] * len(st.parent), [0.0] * len(st.parent)
    for v in range(len(st.parent)):
        if st.parent[v] < 0:
            continue
        lam[v], mu[v] = term if not st.children[v] else inner
    fams, tries, over = simulate(st, lam, mu, nfam, seed, min_species=4, cap=4000,
                                 max_over=nfam, max_tries=200 * nfam)
    species = sorted(st.label[v] for v in st.leaves())
    spi = {s: i for i, s in enumerate(species)}
    base = [true_rt(f, spi) for f in fams]
    topo = st.newick()
    rng = random.Random(seed + 1)
    out = {}
    models = [("d2s", q) for q in QS] + [("flip", q) for q in QS[1:4]] + [("rovl", p) for p in PS] + [("ovl", 0)]
    for model, par in models:
        err = [apply_error(model, par, rt, tg, spi, rng) for rt, tg in base]
        sc = quartet_scores([encode(r, t) for r, t in err], [r for r, _ in err], True, topo)
        per = {}
        for node, t, _, s in sc:
            per.setdefault(node, {})[t] = s
        out["%s:%g" % (model, par)] = [[d["t1"], d["t2"], d["t3"]] for _, d in sorted(per.items())]
    sc = quartet_scores([plain(r) for r, _ in base], [r for r, _ in base], False, topo)
    per = {}
    for node, t, _, s in sc:
        per.setdefault(node, {})[t] = s
    out["own:0"] = [[d["t1"], d["t2"], d["t3"]] for _, d in sorted(per.items())]
    return {"shape": shape, "T": T, "inner": inner, "term": term, "nfam": nfam, "seed": seed,
            "tries": tries, "over": over, "scores": out}


if __name__ == "__main__":
    out, nfam = sys.argv[1], int(sys.argv[2])
    done = set()
    if os.path.exists(out):
        done = {(r["shape"], r["T"], tuple(r["inner"]), tuple(r["term"]), r["seed"])
                for r in map(json.loads, open(out))}
    grid = list(itertools.product(["cat4", "cat5", "bal4"], [0.05, 0.2, 0.5],
                                  [(3.0, 0.5), (1.0, 1.0), (3.0, 3.0)], [(1.0, 1.0), (0.0, 1.0)]))
    for shape, T, inner, term in grid:
        seed = 11
        if (shape, T, inner, term, seed) in done:
            continue
        try:
            r = run(shape, T, inner, term, nfam, seed)
        except RuntimeError as e:
            r = {"shape": shape, "T": T, "inner": inner, "term": term, "seed": seed, "error": str(e)}
        with open(out, "a") as f:
            f.write(json.dumps(r) + "\n")
        if "scores" in r:
            m = {k: [round((a - max(b, c)) / max(a + b + c, 1e-9), 3) for a, b, c in v]
                 for k, v in r["scores"].items()}
            print(shape, T, inner, term, m, flush=True)
