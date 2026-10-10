"""Monte-Carlo check of the 'triplet association' inequalities behind the 4-taxon
argument (results/theory.md).  One gene copy starts at the top of a stem branch of
length t0 above the 3-taxon tree ((A,B)x,C)y.  With indicators a, b, c = 'the copy
leaves a descendant in A / B / C':

  (I1)  P(a b ~c) P(~a ~b c) >= P(a ~b c) P(~a b ~c)   (exact survival sets)
  (I2)  P(a b) P(c)          >= P(a c) P(b)            (containment)

(I1) controls hidden-paralogy quartets under species-overlap tagging; (I2) controls
duplication nodes mislabelled S at random (d2s / flip error models).
usage: python assoc.py N_CONFIGS NREP OUT.jsonl
"""
import json
import math
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gdlsim import simulate_family, leaf_species  # noqa: E402
from phylo import parse_newick  # noqa: E402

LAM = [0, 0.5, 1, 2, 3]
MU = [0, 0.5, 1, 2, 3]
TS = [0.05, 0.2, 0.5, 1.0]


def one(seed, nrep):
    rng = random.Random(seed)
    tx, ta, tb, tc, t0 = (rng.choice(TS) for _ in range(5))
    nwk = "((A:%g,B:%g)x:%g,C:%g)y;" % (ta, tb, tx, tc + tx)
    st = parse_newick(nwk)
    rates = {lab: (rng.choice(LAM), rng.choice(MU)) for lab in ["A", "B", "C", "x", "y"]}
    lam = [rates[st.label[v] or "y"][0] for v in range(len(st.parent))]
    mu = [rates[st.label[v] or "y"][1] for v in range(len(st.parent))]
    cnt = {}
    srng = random.Random(seed * 7 + 1)
    over = 0
    for _ in range(nrep):
        g = simulate_family(st, lam, mu, srng, root_len=t0, cap=3000)
        if g == "OVERFLOW":
            over += 1
            continue
        s = frozenset() if g is None else frozenset(leaf_species(g))
        cnt[s] = cnt.get(s, 0) + 1
    n = sum(cnt.values())
    P = lambda *sets: sum(cnt.get(frozenset(x), 0) for x in sets) / n
    pabnc, pnnc, panc, pnbn = P("AB"), P("C"), P("AC"), P("B")
    pab = P("AB", "ABC")
    pac = P("AC", "ABC")
    pa_b = P("B", "AB", "BC", "ABC")
    pc = P("C", "AC", "BC", "ABC")
    i1 = pabnc * pnnc - panc * pnbn
    i2 = pab * pc - pac * pa_b
    # standard errors via delta method are messy; report a crude binomial scale
    se1 = math.sqrt(max(pabnc * pnnc, panc * pnbn, 1e-12) / n) * 2
    se2 = math.sqrt(max(pab * pc, pac * pa_b, 1e-12) / n) * 2
    return {"seed": seed, "tree": nwk, "t0": t0, "rates": rates, "n": n, "over": over,
            "I1": i1, "I1_se": se1, "I2": i2, "I2_se": se2,
            "dist": {"".join(sorted(k)) or "-": v for k, v in cnt.items()}}


if __name__ == "__main__":
    nconf, nrep, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    s0 = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    for seed in range(s0, s0 + nconf):
        r = one(seed, nrep)
        with open(out, "a") as f:
            f.write(json.dumps(r) + "\n")
        flag = "VIOL" if (r["I1"] < -r["I1_se"] or r["I2"] < -r["I2_se"]) else ""
        print(r["seed"], round(r["I1"], 5), round(r["I1_se"], 5), round(r["I2"], 5), round(r["I2_se"], 5), flag, flush=True)
