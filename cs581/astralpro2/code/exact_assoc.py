"""Exact (generating-function) check of the triplet association inequalities
(see results/theory.md, Section 3) for a single gene copy entering the stem of
((A,B)x,C)y under a linear birth-death (duplication lam, loss mu) process with
branch-specific rates.

For a copy at the top of branch e above node v: N copies reach v with pgf phi_e, and
the surviving species set is the union of N iid single-copy sets, so
F_e(T) := P(set within T) = phi_e(F_v(T)), and F_v(T) = prod_children F_c(T & L_c).

  I1 = P(=AB) P(=C) - P(=AC) P(=B)        (exact sets)
  I2 = P(>=AB) P(>=C) - P(>=AC) P(>=B)    (containment)
usage: python exact_assoc.py N_CONFIGS OUT.jsonl
"""
import itertools
import json
import math
import random
import sys


def phi(s, lam, mu, t):
    if t <= 0 or (lam == 0 and mu == 0):
        return s
    if abs(lam - mu) < 1e-12:
        a = lam * t
        return (a * (1 - s) + s) / (a * (1 - s) + 1)
    e = math.exp(-(lam - mu) * t)
    return (mu * (1 - s) - (mu - lam * s) * e) / (lam * (1 - s) - (mu - lam * s) * e)


def F_tree(T, br):
    """br: dict of (lam, mu, t) for branches 'A','B','C','x','stem'. T subset of {'A','B','C'}"""
    leaf = lambda X: phi(1.0 if X in T else 0.0, *br[X])
    Fx = phi(leaf("A") * leaf("B"), *br["x"])
    return phi(Fx * leaf("C"), *br["stem"])


def probs(br):
    U = ("A", "B", "C")
    F = {frozenset(T): F_tree(set(T), br) for r in range(4) for T in itertools.combinations(U, r)}
    # exact-set probabilities by Moebius inversion
    exact = {}
    for S in F:
        exact[S] = sum((-1) ** (len(S) - len(R)) * F[R] for R in F if R <= S)
    contains = lambda *xs: sum(p for S, p in exact.items() if set(xs) <= S)
    i1 = exact[frozenset("AB")] * exact[frozenset("C")] - exact[frozenset("AC")] * exact[frozenset("B")]
    i2 = contains("A", "B") * contains("C") - contains("A", "C") * contains("B")
    return i1, i2, exact


LAM = [0, 0.25, 0.5, 1, 2, 4, 8]
MU = [0, 0.25, 0.5, 1, 2, 4, 8]
TS = [0.01, 0.05, 0.2, 0.5, 1.0, 2.0, 4.0]

if __name__ == "__main__":
    n, out = int(sys.argv[1]), sys.argv[2]
    rng = random.Random(1)
    worst1 = worst2 = (1, None)
    with open(out, "w") as f:
        for i in range(n):
            br = {k: (rng.choice(LAM), rng.choice(MU), rng.choice(TS)) for k in ["A", "B", "C", "x", "stem"]}
            i1, i2, ex = probs(br)
            # scale-free versions: divide by the larger product
            d1 = max(ex[frozenset("AB")] * ex[frozenset("C")], ex[frozenset("AC")] * ex[frozenset("B")], 1e-300)
            r1, r2 = i1 / d1, i2
            if r1 < worst1[0]:
                worst1 = (r1, br)
            if r2 < worst2[0]:
                worst2 = (r2, br)
            if i < 2000 or i1 < -1e-12 or i2 < -1e-12:
                f.write(json.dumps({"br": br, "I1": i1, "I1_rel": r1, "I2": i2}) + "\n")
    print("min relative I1:", worst1)
    print("min I2:", worst2)
