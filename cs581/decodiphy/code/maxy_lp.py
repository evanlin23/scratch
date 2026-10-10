"""Test: does the single LP  max ybar  s.t.  F m + ybar 1 = d, 1'm = 1, m >= 0  (m on ALL nodes)
recover the true node measure m (hence the placements) from noise-free d?
Motivation: Theorem 1 - all solutions differ by Laplacian moves, and spreading mass to a node's
neighbours (spurious star split) always lowers ybar.
usage: python maxy_lp.py [ntrials]
"""
import sys, itertools
import numpy as np
from scipy.optimize import linprog
from pdd import Tree


def random_tree(n, rng):
    edges = [(0, n), (1, n), (2, n)]; nxt = n + 1
    for leaf in range(3, n):
        i = rng.integers(len(edges)); a, b = edges.pop(i)
        edges += [(a, nxt), (nxt, b), (nxt, leaf)]; nxt += 1
    return Tree(n, [(a, b, l) for (a, b), l in zip(edges, rng.exponential(1, len(edges)) + 0.05)])


def node_measure(t, S, p, x):
    m = np.zeros(t.N)
    for e, pq, xq in zip(S, p, x):
        a, b, _ = t.edges[e]
        m[a] += pq * (1 - xq); m[b] += pq * xq
    return m


def maxy(t, d):
    V = t.N
    F = t.D[: t.n, :V]
    A_eq = np.zeros((t.n + 1, V + 1)); A_eq[: t.n, :V] = F; A_eq[: t.n, V] = 1; A_eq[t.n, :V] = 1
    c = np.zeros(V + 1); c[V] = -1
    res = linprog(c, A_eq=A_eq, b_eq=np.concatenate([d, [1]]), bounds=[(0, None)] * (V + 1), method="highs")
    return res.x[:V], res.x[V]


def main():
    rng = np.random.default_rng(3)
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    import collections
    stats = collections.Counter()
    for trial in range(T):
        n = int(rng.integers(6, 41))
        t = random_tree(n, rng)
        k = int(rng.integers(1, 6))
        nbr = collections.defaultdict(set)
        for a, b, _ in t.edges:
            nbr[a].add(b); nbr[b].add(a)
        # random matching S of size k
        for _ in range(200):
            S = list(rng.choice(t.m, k, replace=False))
            VS = [v for e in S for v in t.edges[e][:2]]
            if len(set(VS)) == 2 * k:
                break
        else:
            continue
        closed = any(nbr[v] <= set(VS) and v in VS for v in range(t.n, t.N))
        opened = any(nbr[v] <= set(VS) and v not in VS for v in range(t.n, t.N))
        claw = "closed-claw" if closed else ("open-claw" if opened else "no-claw")
        p = rng.dirichlet(np.ones(k)); x = rng.uniform(0.05, 0.95, k); y = rng.uniform(0.0, 0.5)
        d = t.mixture(S, p, x, y)
        m_true = node_measure(t, S, p, x)
        m_hat, y_hat = maxy(t, d)
        ok = np.abs(m_hat - m_true).max() < 1e-6
        key = (k, claw)
        stats[key + ("N",)] += 1
        stats[key + ("recovered",)] += ok
        stats[key + ("ybar_hat>=ybar",)] += y_hat >= y - 1e-9
    for key in sorted({k[:2] for k in stats}):
        print(key, "recovered %d/%d" % (stats[key + ("recovered",)], stats[key + ("N",)]),
              "ybar_hat>=true %d" % stats[key + ("ybar_hat>=ybar",)])


if __name__ == "__main__":
    main()
