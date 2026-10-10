"""Error-vs-families curve for ASTRAL-Pro (full search, not just scoring) on one
4-taxon caterpillar, under tagging-error models: fraction of disjoint replicate
datasets of K families on which ASTRAL-Pro returns a wrong species tree.

usage: python curve4.py SPEC POOLSIZE OUT.jsonl models KS
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core  # noqa: E402
from phylo import parse_newick  # noqa: E402

spec, pool_n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
models = sys.argv[4].split(",")
ks = [int(x) for x in sys.argv[5].split(",")]
br = {k: tuple(float(x) for x in v.split(",")) for k, v in (kv.split("=") for kv in spec.split(";"))}
nwk = "(((A:%g,B:%g)x:%g,C:%g)y:%g,D:%g);" % (br["A"][2], br["B"][2], br["x"][2], br["C"][2], br["y"][2], br["D"][2])
st = parse_newick(nwk)
lam, mu = [0.0] * len(st.parent), [0.0] * len(st.parent)
for v in range(len(st.parent)):
    if st.label[v] in br:
        lam[v], mu[v] = br[st.label[v]][0], br[st.label[v]][1]
cache = out + ".pool.nwk"
if os.path.exists(cache):
    fams = [l.strip() for l in open(cache) if ";" in l]
else:
    fams, _, _ = core.simulate(st, lam, mu, pool_n, 4242, min_species=4, cap=6000, max_over=pool_n,
                               max_tries=3000 * pool_n)
    open(cache, "w").write("\n".join(fams) + "\n")
spi = {s: i for i, s in enumerate("ABCD")}
base = [core.true_rt(f, spi) for f in fams]
done = set()
if os.path.exists(out):
    done = {(r["model"], r["K"], r["block"]) for r in map(json.loads, open(out))}
for m in models:
    model, par = m.split(":")
    for k in ks:
        nb = len(base) // k
        for b in range(min(nb, 40)):
            if (m, k, b) in done:
                continue
            blk = base[b * k:(b + 1) * k]
            if model == "own":
                est = core.run_apro([core.plain(r) for r, _ in blk], [r for r, _ in blk], fixed=False)
            elif model == "multi":
                est = core.run_apro([core.plain(r) for r, _ in blk], [r for r, _ in blk], fixed=False, binary=core.ASTRAL)
            else:
                rng = random.Random(1000 * k + b)
                err = [core.apply_error(model, float(par), rt, tg, spi, rng) for rt, tg in blk]
                est = core.run_apro([core.encode(r, t) for r, t in err], [r for r, _ in err], fixed=True)
            fn = core.rf_error(est, st)[0]
            with open(out, "a") as f:
                f.write(json.dumps({"model": m, "K": k, "block": b, "FN": fn}) + "\n")
        print(m, k, flush=True)
