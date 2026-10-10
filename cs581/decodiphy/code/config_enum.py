r"""Configuration-level (parameter-free) identifiability check, exhaustive over unlabeled tree shapes.

By Theorem 1 (see REPORT.md) two solutions with node measures m, m' give the same d iff
m' - m = g := -L alpha for some alpha supported on internal nodes (ybar' = ybar - sum alpha).
For a true edge set S (a matching) and a competitor S' (|S'| <= |S|, S' != S), there exist
parameters (p, x, ybar of the truth) making S' an exact alternative iff some alpha has
    g(u) = 0   for u outside V(S) u V(S'),
    g(u) > 0   for u in V(S') \ V(S),
    g(u) < 0   for u in V(S) \ V(S'),
(free on V(S) n V(S')), since alpha can be scaled down and ybar chosen large.
The strict inequalities are homogeneous, so we test feasibility with |g| >= 1. Also report whether
the true S alone has a continuum (kernel of L[V \ V(S), I] nontrivial).

usage: python config_enum.py NMAX OUT.md
"""
import sys, itertools, collections
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, __import__("os").path.dirname(__file__))
from pdd import all_topologies


def canon(n, edges):
    """canonical string of an unlabeled unrooted tree (AHU encoding over all centres)."""
    adj = collections.defaultdict(list)
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)

    def enc(u, p):
        return "(" + "".join(sorted(enc(v, u) for v in adj[u] if v != p)) + ")"
    return min(enc(r, -1) for r in adj if len(adj[r]) > 1)


def shapes(n):
    seen = {}
    for E in all_topologies(n):
        c = canon(n, E)
        if c not in seen:
            seen[c] = E
    return list(seen.values())


def laplacian(N, edges, lens):
    L = np.zeros((N, N))
    for (a, b), l in zip(edges, lens):
        L[a, a] += 1 / l; L[b, b] += 1 / l; L[a, b] -= 1 / l; L[b, a] -= 1 / l
    return L


def alt_exists(L, n, VS, VSp):
    N = L.shape[0]
    G = -L[:, n:]  # g = G @ alpha
    zero = [u for u in range(N) if u not in VS and u not in VSp]
    pos = [u for u in VSp if u not in VS]
    neg = [u for u in VS if u not in VSp]
    if not pos and not neg:
        return None
    A_ub = np.vstack([-G[pos], G[neg]]) if pos or neg else None
    b_ub = -np.ones(len(pos) + len(neg))
    A_eq = G[zero] if zero else None
    b_eq = np.zeros(len(zero)) if zero else None
    res = linprog(np.zeros(N - n), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=[(None, None)] * (N - n), method="highs")
    return res.status == 0


def main():
    NMAX = int(sys.argv[1]); OUT = sys.argv[2]
    rng = np.random.default_rng(7)
    lines = ["| n | shape | k | #S (matchings) | #S with continuum | #S with alternative k'<=k (unit lengths) | (random lengths) | claw-free S with alt | S with claw but no alt | claw centred in V(S) | alt <=> continuum? |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    detail = []
    for n in range(4, NMAX + 1):
        for si, E in enumerate(shapes(n)):
            N = 2 * n - 2
            m = len(E)
            nbr = collections.defaultdict(set)
            for a, b in E:
                nbr[a].add(b); nbr[b].add(a)
            Ls = {"unit": laplacian(N, E, np.ones(m)), "rand": laplacian(N, E, rng.exponential(1, m) + 0.05)}
            for k in ((1, 2, 3, 4) if n <= 8 else (1, 2, 3)):
                if 2 * k > N:
                    continue
                cnt = collections.Counter()
                for S in itertools.combinations(range(m), k):
                    VS = [x for e in S for x in E[e]]
                    if len(set(VS)) != len(VS):
                        continue  # only matchings (adjacent pairs always give a continuum)
                    VS = set(VS)
                    cnt["S"] += 1
                    claw = any(nbr[v] <= VS for v in range(n, N))
                    clawin = any(nbr[v] <= VS and v in VS for v in range(n, N))
                    cnt["clawin"] += clawin
                    outside = [u for u in range(N) if u not in VS]
                    cont = np.linalg.matrix_rank(Ls["unit"][outside][:, n:]) < N - n if outside else True
                    cnt["cont"] += cont
                    res = {}
                    for key, L in Ls.items():
                        found = False
                        for kp in range(1, k + 1):
                            for Sp in itertools.combinations(range(m), kp):
                                if set(Sp) == set(S):
                                    continue
                                VSp = {x for e in Sp for x in E[e]}
                                if alt_exists(L, n, VS, VSp):
                                    found = True; break
                            if found:
                                break
                        res[key] = found
                    cnt["alt_unit"] += res["unit"]; cnt["alt_rand"] += res["rand"]
                    cnt["noclaw_alt"] += (not claw) and (res["unit"] or res["rand"])
                    cnt["claw_noalt"] += claw and not (res["unit"] or res["rand"])
                    cnt["mismatch"] += (bool(cont) != bool(res["unit"]))
                    if (not claw) and (res["unit"] or res["rand"]) and len(detail) < 10:
                        detail.append(f"n={n} shape={si} edges={E} S={[E[e] for e in S]}")
                lines.append(f"| {n} | {si} | {k} | {cnt['S']} | {cnt['cont']} | {cnt['alt_unit']} | {cnt['alt_rand']} | {cnt['noclaw_alt']} | {cnt['claw_noalt']} | {cnt['clawin']} | {'yes' if cnt['mismatch'] == 0 else 'no (%d)' % cnt['mismatch']} |")
                print(lines[-1], flush=True)
    txt = "\n".join(lines) + "\n\nclaw-free counterexamples:\n" + "\n".join(detail) + "\n"
    open(OUT, "w").write(txt)


if __name__ == "__main__":
    main()
