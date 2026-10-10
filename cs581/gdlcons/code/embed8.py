"""Embed the 4-taxon counterexample cand2 in an 8-taxon tree: does ASTRAL-Pro3 (own
rooting/tagging, true gene trees) converge to a wrong 8-taxon tree?
usage: python embed8.py OUT.jsonl K NREP"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core
from phylo import parse_newick
out, K, nrep = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
which = sys.argv[4] if len(sys.argv) > 4 else "cand2"
if which == "cand2":
    nwk = "(((((A:1,B:1)x:4,C:4)y:0.5,(D1:1,D2:1)u:0.5)r1:0.5,E:1)r2:0.5,(F1:1,F2:1)v:1);"
    rates = {"y": (8, 0.5), "A": (1, 4), "B": (4, 8), "C": (0.5, 1)}  # everything else event-free
else:  # cand4: ASTRAL-Pro's own rooting fails
    nwk = "(((((A:0.5,B:0.5)x:0.5,C:0.2)y:1,(D1:1,D2:1)u:0.5)r1:0.5,E:1)r2:0.5,(F1:1,F2:1)v:1);"
    rates = {"y": (4, 1), "x": (1, 0.5), "A": (0, 8), "B": (0, 8), "C": (8, 8)}
st = parse_newick(nwk)
lam = [rates.get(st.label[v], (0, 0))[0] for v in range(len(st.parent))]
mu = [rates.get(st.label[v], (0, 0))[1] for v in range(len(st.parent))]
species = sorted(st.label[v] for v in st.leaves())
spi = {s: i for i, s in enumerate(species)}
for rep in range(nrep):
    fams, tries, _ = core.simulate(st, lam, mu, K, 900 + rep, min_species=4, cap=6000, max_over=K)
    base = [core.true_rt(f, spi) for f in fams]
    rec = {"K": K, "rep": rep}
    for m in ["own", "true", "ovl"]:
        if m == "own":
            est = core.run_apro([core.plain(r) for r, _ in base], [r for r, _ in base], fixed=False)
        else:
            import random
            err = [core.apply_error(m, 0, rt, tg, spi, random.Random(rep)) for rt, tg in base]
            est = core.run_apro([core.encode(r, t) for r, t in err], [r for r, _ in err], fixed=True)
        rec[m] = {"FN": core.rf_error(est, st)[0], "tree": est.newick()}
    with open(out, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec), flush=True)
