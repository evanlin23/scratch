"""Closed-form limiting ASTRAL-Pro scores in 'family R' (see ../theory.md):
caterpillar (((A,B)x,C)y,D); y-branch pure birth (rate lam, length T); no events on x, D;
tips A,B,C pure loss with per-copy survival a,b,c.  Correct root + species-overlap tags.
Per accepted-or-not family, divided by E[#copies at y] = e^{lam T}:
  score(AB|CD)/e^{lT} = abc + 2 int_{q0}^1 p_q(AB) p_q(C) dq
  score(AC|BD)/e^{lT} =       2 int_{q0}^1 p_q(AC) p_q(B) dq,   q0 = e^{-lam T}
p_q(S) = P(union of k ~ Geom(q) copies has species set exactly S)."""
import itertools, math
import numpy as np
from scipy.integrate import quad

SP = "ABC"

def G(z, q):  # E[z^k], k ~ Geometric(q) on {1,2,...}
    return z * q / (1 - (1 - q) * z)

def p_exact(S, surv, q):
    # surv: dict X->survival; P(union exactly S) by inclusion-exclusion over U subset S
    tot = 0.0
    S = list(S)
    for r in range(len(S) + 1):
        for U in itertools.combinations(S, r):
            z = 1.0
            for X in SP:
                if X not in U:
                    z *= 1 - surv[X]
            tot += (-1) ** (len(S) - r) * G(z, q)
    return tot

def scores(a, b, c, lamT=math.inf):
    s = {"A": a, "B": b, "C": c}
    q0 = math.exp(-lamT) if lamT < math.inf else 0.0
    f = lambda S1, S2: 2 * quad(lambda q: p_exact(S1, s, q) * p_exact(S2, s, q), q0, 1, epsabs=1e-14, epsrel=1e-12, limit=200)[0]
    return {"O": a * b * c, "AB": f("AB", "C"), "AC": f("AC", "B"), "BC": f("BC", "A")}

def margin(a, b, c, lamT=math.inf):
    r = scores(a, b, c, lamT)
    cor = r["O"] + r["AB"]
    return cor - max(r["AC"], r["BC"]), r

if __name__ == "__main__":
    rng = np.random.default_rng(1)
    best = []
    for it in range(4000):
        a, b, c = 10 ** rng.uniform(-3, 0, 3)
        m, r = margin(a, b, c)
        tot = r["O"] + r["AB"] + r["AC"] + r["BC"]
        best.append((m / tot, a, b, c, r))
    best.sort(key=lambda x: x[0])
    for row in best[:15]:
        print("%.4f a=%.4g b=%.4g c=%.4g" % row[:4], {k: "%.3g" % v for k, v in row[4].items()})
    print("frac negative", sum(1 for x in best if x[0] < 0) / len(best))
