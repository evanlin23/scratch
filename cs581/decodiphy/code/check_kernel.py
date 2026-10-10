"""Numerically verify the kernel characterization:
   map (m in R^V, ybar) -> (F m + ybar*1, sum m), F = leaf-by-node distance matrix,
   has rank n+1 and kernel = span{(-L e_v, -1): v internal}, L = Laplacian with conductances 1/l."""
import numpy as np, random
from pdd import Tree, all_topologies
rng = np.random.default_rng(0)
worst = 0
for n in range(4, 13):
    for rep in range(20):
        # random topology by stepwise addition
        edges = [(0, n), (1, n), (2, n)]; nxt = n + 1
        for leaf in range(3, n):
            i = rng.integers(len(edges)); a, b = edges.pop(i)
            edges += [(a, nxt), (nxt, b), (nxt, leaf)]; nxt += 1
        t = Tree(n, [(a, b, l) for (a, b), l in zip(edges, rng.exponential(1, len(edges)) + 0.01)])
        V = t.N
        F = t.D[:n, :V]
        A = np.zeros((n + 1, V + 1)); A[:n, :V] = F; A[:n, V] = 1; A[n, :V] = 1
        r = np.linalg.matrix_rank(A)
        L = np.zeros((V, V))
        for a, b, l in t.edges:
            L[a, a] += 1 / l; L[b, b] += 1 / l; L[a, b] -= 1 / l; L[b, a] -= 1 / l
        K = np.array([np.concatenate([-L[:, v], [-1.0]]) for v in range(n, V)]).T
        assert r == n + 1, (n, r)
        assert np.linalg.matrix_rank(K) == V - n
        worst = max(worst, np.abs(A @ K).max())
print("rank n+1 and kernel = Laplacian star moves verified; max |A K| =", worst)
