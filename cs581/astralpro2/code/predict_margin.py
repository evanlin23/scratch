"""Predicted per-family ASTRAL-Pro quartet supports for the 4-taxon caterpillar
(((A,B)x,C)y,D) under pure GDL when D is single-copy and duplications happen only on
the y-branch (root -> y, length T, rates lam, mu) -- see results/theory.md.

  O       = E[# orthologous SQ classes]            = m(T) * P(A&B | x) * s_C
  H_tau   = E[# hidden-paralog nodes with pattern tau]
          = int_0^T 2 lam m(t) p_{T-t}(S1) p_{T-t}(S2) dt
with p_s(S) the probability that one copy at distance s above y leaves exactly species
set S, m(t) = exp((lam - mu) t).  Correct pattern ({A,B},{C}); wrong ({A,C},{B}),
({B,C},{A}).
  true tags:   correct O,                wrong 0
  ovl / d2sv(q): correct O + q H_AB,     wrong q H_AC, q H_BC   (q = 1 for ovl)
usage: python predict_margin.py 'A=..;B=..;C=..;x=..;y=..' [q ...]
"""
import math
import sys

from exact_assoc import probs, phi


def predict(br, n=400):
    lam, mu, T = br["y"]
    m = lambda t: math.exp((lam - mu) * t)
    # O: one copy at y; x-side needs A and B, C-side needs C
    sA = 1 - phi(0.0, *br["A"])
    sB = 1 - phi(0.0, *br["B"])
    sC = 1 - phi(0.0, *br["C"])
    # P(A&B) for a copy entering x: 1 - F(no A) - F(no B) + F(none)
    fx = lambda a, b: phi(phi(a, *br["A"]) * phi(b, *br["B"]), *br["x"])
    pAB_x = 1 - fx(0, 1) - fx(1, 0) + fx(0, 0)
    O = m(T) * pAB_x * sC
    H = {"AB": 0.0, "AC": 0.0, "BC": 0.0}
    dt = T / n
    for i in range(n):
        t = (i + 0.5) * dt
        b2 = dict(br)
        b2["stem"] = (lam, mu, T - t)
        _, _, ex = probs(b2)
        e = lambda s: ex[frozenset(s)]
        w = 2 * lam * m(t) * dt
        H["AB"] += w * e("AB") * e("C")
        H["AC"] += w * e("AC") * e("B")
        H["BC"] += w * e("BC") * e("A")
    return O, H


if __name__ == "__main__":
    br = {k: tuple(float(x) for x in v.split(",")) for k, v in (kv.split("=") for kv in sys.argv[1].split(";"))}
    qs = [float(x) for x in sys.argv[2:]] or [1.0]
    O, H = predict(br)
    print("O = %.5g  H_AB = %.5g  H_AC = %.5g  H_BC = %.5g" % (O, H["AB"], H["AC"], H["BC"]))
    for q in qs:
        c, w1, w2 = O + q * H["AB"], q * H["AC"], q * H["BC"]
        print("q=%.3g  correct %.5g  AC|BD %.5g  AD|BC %.5g  c-AC %.5g  c-BC %.5g" % (q, c, w1, w2, c - w1, c - w2))
    d = max(H["AC"], H["BC"]) - H["AB"]
    print("threshold q* = O / (H_wrong - H_AB) =", O / d if d > 0 else "none (consistent for all q)")
