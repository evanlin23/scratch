"""Smallest y-branch duplication load lamT at which family R can be inconsistent
(correct root + overlap tags), grid over a=b, c with a,b >= 1e-3; exact interval formula."""
import json, math
from fractions import Fraction as Fr
from exact_proof import run

res = {}
for lamT in [0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 6.0]:
    q0 = Fr(math.exp(-lamT)).limit_denominator(10 ** 12)
    best = None
    for la in [-3, -2.5, -2, -1.5, -1.3, -1]:
        for c in [0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6]:
            a = Fr(10 ** la).limit_denominator(10 ** 6)
            if Fr(c) <= a: continue
            r = run(a, a, Fr(c).limit_denominator(1000), q0)
            tot = r["O"] + r["H_AB"] + r["H_AC"] + r["H_BC"]
            m = float(((r["O"] + r["H_AB"] - r["H_AC"]) / tot).mid)
            if best is None or m < best[0]:
                best = (m, float(a), c)
    res[lamT] = best
    print(lamT, "min rel margin %.4f at a=b=%.4g c=%g" % best, flush=True)
json.dump(res, open("../results/famR_minlam.json", "w"), indent=1)
