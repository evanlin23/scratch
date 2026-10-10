"""ASTRAL-Pro quartet supports on one 4-taxon caterpillar (((A,B)x,C)y,D) with
branch-specific rates (e.g. an association-inequality violator from exact_assoc.py).
The 'stem' of the 3-taxon analysis is the y-branch (root -> y).

usage: python probe4.py OUT.jsonl NFAM SEED 'A=lam,mu,t;B=...;C=...;x=...;y=...;D=...' [models]
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import apply_error, encode, plain, quartet_scores, simulate, true_rt  # noqa: E402
from core import split_of as core_split  # noqa: E402
from phylo import parse_newick  # noqa: E402

DEFAULT_MODELS = "true:0,ovl:0,d2s:0.1,d2s:0.25,d2s:0.5,d2s:0.75,d2s:1,flip:0.1,flip:0.25,flip:0.5,rovl:0.3,rovl:1,own:0"


def main():
    out, nfam, seed, spec = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    models = (sys.argv[5] if len(sys.argv) > 5 else DEFAULT_MODELS).split(",")
    br = {k: tuple(float(x) for x in v.split(",")) for k, v in (kv.split("=") for kv in spec.split(";"))}
    nwk = "(((A:%g,B:%g)x:%g,C:%g)y:%g,D:%g);" % (br["A"][2], br["B"][2], br["x"][2], br["C"][2],
                                                br["y"][2], br["D"][2])
    st = parse_newick(nwk)
    lam, mu = [0.0] * len(st.parent), [0.0] * len(st.parent)
    for v in range(len(st.parent)):
        lab = st.label[v]
        if lab in br:
            lam[v], mu[v] = br[lab][0], br[lab][1]
    fams, tries, over = simulate(st, lam, mu, nfam, seed, min_species=4, cap=6000, max_over=nfam,
                                 max_tries=2000 * nfam)
    species = ["A", "B", "C", "D"]
    spi = {s: i for i, s in enumerate(species)}
    base = [true_rt(f, spi) for f in fams]
    rng = random.Random(seed + 1)
    res = {}
    for m in models:
        model, par = m.split(":")
        if model == "own":
            sc = quartet_scores([plain(r) for r, _ in base], [r for r, _ in base], False, st.newick())
        else:
            err = [apply_error(model, float(par), rt, tg, spi, rng) for rt, tg in base]
            sc = quartet_scores([encode(r, t) for r, t in err], [r for r, _ in err], True, st.newick())
        d = {core_split(lab): s for _, _, lab, s in sc}
        res[m] = [d.get("AB|CD", 0.0), d.get("AD|BC", 0.0), d.get("AC|BD", 0.0)]  # correct, AD|BC, AC|BD
    rec = {"spec": spec, "nfam": nfam, "seed": seed, "tries": tries, "over": over,
           "mean_leaves": sum(len(r.leaves()) for r in [b[0] for b in base]) / len(base), "scores": res}
    with open(out, "a") as f:
        f.write(json.dumps(rec) + "\n")
    for m, (a, b, c) in res.items():
        print("%-10s correct %12.1f  wrong %12.1f %12.1f  margin %+.4f" % (m, a, b, c, (a - max(b, c)) / max(a + b + c, 1e-9)))


if __name__ == "__main__":
    main()
