"""Family R phase boundary (correct root + overlap tags, exact formula):
(1) min over (a,b,c) of the relative margin, for each lamT (smallest lamT allowing failure);
(2) for a=b=s, the critical c above which ASTRAL-Pro is inconsistent (lamT = 3, inf)."""
import math, json
import numpy as np
from scipy.optimize import minimize, brentq
from famR import margin

def rel(x, lamT):
    a, b, c = 1 / (1 + np.exp(-np.asarray(x)))
    m, r = margin(a, b, c, lamT)
    tot = r["O"] + r["AB"] + r["AC"] + r["BC"]
    return m / tot

out = {}
for lamT in [0.5, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 5.0, math.inf]:
    best = (9, None)
    for x0 in [(-4, -4, -1), (-5, -5, -1.5), (-3, -6, -1), (-6, -6, -2), (-2, -2, 0)]:
        res = minimize(rel, x0, args=(lamT,), method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-6, "maxiter": 400})
        if res.fun < best[0]:
            best = (res.fun, (1 / (1 + np.exp(-res.x))).tolist())
    out[str(lamT)] = best
    print("lamT", lamT, "min rel margin %.4f at a,b,c=" % best[0], ["%.4g" % v for v in best[1]], flush=True)
crit = {}
for lamT in [2.0, 3.0, math.inf]:
    row = {}
    for s in [0.005, 0.01, 0.02, 0.05, 0.1]:
        f = lambda c: margin(s, s, c, lamT)[0]
        cs = np.linspace(s * 1.01, 0.99, 60)
        vals = [f(c) for c in cs]
        neg = [c for c, v in zip(cs, vals) if v < 0]
        row[s] = (min(neg), max(neg)) if neg else None
    crit[str(lamT)] = row
    print("lamT", lamT, "a=b=s: failing c range", row, flush=True)
json.dump({"min_margin": out, "crit": {k: {str(s): v for s, v in r.items()} for k, r in crit.items()}},
          open("../results/famR_phase.json", "w"), indent=1)
