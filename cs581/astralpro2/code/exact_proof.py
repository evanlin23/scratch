"""Exact evaluation (rational coefficients + logs, interval arithmetic) of the limiting
ASTRAL-Pro scores in family R with correct root + species-overlap tags (theory.md, Thm 1).

Per family, divided by E[#copies at y] = 1/q0 (q0 = e^{-lam T}):
  O    = a b c
  H_S1|S2 = 2 int_{q0}^{1} p_q(S1) p_q(S2) dq
  p_q(S) = sum_{U subset S} (-1)^{|S-U|} G(z_U),  z_U = prod_{X notin U} (1 - s_X),
  G(z) = E[z^k] = q z / (1 - (1-q) z) = q / (q + beta),  beta = (1-z)/z   (G(1) = 1).
Integrals of products of q/(q+beta1) * q/(q+beta2) are done in closed form and
evaluated in mpmath interval arithmetic, so the printed enclosures are rigorous.
usage: python exact_proof.py a b c q0      (rationals, e.g. 1/100 1/100 1/5 1/20)"""
import itertools
import sys
from fractions import Fraction as Fr

from mpmath import iv

iv.dps = 50


def terms(S, surv):
    """p_q(S) as a list of (coef, beta) meaning coef * q/(q + beta); beta None means constant coef.
    G(z) = E[z^k] = q z / (1 - (1-q) z) = q / (q + beta), beta = (1-z)/z  (G(1) = 1)."""
    out = []
    for r in range(len(S) + 1):
        for U in itertools.combinations(S, r):
            z = Fr(1)
            for X in "ABC":
                if X not in U:
                    z *= 1 - surv[X]
            sign = (-1) ** (len(S) - r)
            out.append((sign, None) if z == 1 else (sign, (1 - z) / z))
    return out


def I(x):
    x = Fr(x)
    return iv.mpf(x.numerator) / iv.mpf(x.denominator)


def integ(b1, b2, q0):
    """int_{q0}^1 f1 f2 dq with f = q/(q+beta) = 1 - beta/(q+beta), or f = 1 when beta is None."""
    lo, hi = I(q0), iv.mpf(1)
    dlog = lambda b: iv.log(hi + I(b)) - iv.log(lo + I(b))        # int 1/(q+b)
    if b1 is None and b2 is None:
        return hi - lo
    if b1 is None or b2 is None:
        b = b1 if b2 is None else b2
        return (hi - lo) - I(b) * dlog(b)
    if b1 == b2:
        prod = 1 / (lo + I(b1)) - 1 / (hi + I(b1))                  # int 1/(q+b)^2
    else:   # 1/((q+b1)(q+b2)) = [1/(q+b1) - 1/(q+b2)] / (b2 - b1)
        prod = (dlog(b1) - dlog(b2)) / I(b2 - b1)
    return (hi - lo) - I(b1) * dlog(b1) - I(b2) * dlog(b2) + I(b1) * I(b2) * prod


def H(S1, S2, surv, q0):
    tot = iv.mpf(0)
    for c1, w1 in terms(S1, surv):
        for c2, w2 in terms(S2, surv):
            tot += I(2 * c1 * c2) * integ(w1, w2, q0)
    return tot


def run(a, b, c, q0):
    surv = {"A": Fr(a), "B": Fr(b), "C": Fr(c)}
    return {"O": I(Fr(a) * b * c), "H_AB": H("AB", "C", surv, q0), "H_AC": H("AC", "B", surv, q0),
            "H_BC": H("BC", "A", surv, q0)}


if __name__ == "__main__":
    a, b, c, q0 = [Fr(x) for x in sys.argv[1:5]]
    r = run(a, b, c, q0)
    for k, v in r.items():
        print("%-5s in [%s, %s]" % (k, v.a, v.b))
    for name, w in [("AC|BD", "H_AC"), ("AD|BC", "H_BC")]:
        d = r["O"] + r["H_AB"] - r[w]
        verdict = "NEGATIVE (inconsistent)" if d.b < 0 else ("positive" if d.a > 0 else "undetermined")
        print("score(AB|CD) - score(%s) in [%s, %s]  -> %s" % (name, d.a, d.b, verdict))
