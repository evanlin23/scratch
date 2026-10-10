"""Exact noise-free phylogenetic distance deconvolution (PDD) on small unrooted trees.

Model (Arasti et al. 2026, Eq. 1-2): queries q=1..k attach at points pi_q in the interior of
distinct edges of the reference tree R with abundances p (sum 1); the observed mixture distance
is d_r = sum_q p_q d_R(r, pi_q) + ybar.

For an edge e=(a,b) of length l and a point at distance t from a:
   d(r, point) = d(r,a) + t          if r is on the a side
               = d(r,b) + l - t      if r is on the b side
so with w = p*t the contribution is p*base_e(r) + sign_e(r)*w, sign = +1 (a side) / -1 (b side).
For a fixed edge set S the system is linear in (p, w, ybar):  M_S z = [d; 1].
"""
import itertools
import numpy as np
from scipy.optimize import linprog


class Tree:
    """Unrooted tree, leaves 0..n-1, internal nodes n.., edges list of (a, b, length)."""

    def __init__(self, n, edges):
        self.n = n
        self.edges = list(edges)
        self.N = 1 + max(max(a, b) for a, b, _ in self.edges)
        adj = [[] for _ in range(self.N)]
        for i, (a, b, l) in enumerate(self.edges):
            adj[a].append((b, l, i))
            adj[b].append((a, l, i))
        self.adj = adj
        D = np.zeros((self.N, self.N))
        for s in range(self.N):
            stack = [(s, -1, 0.0)]
            while stack:
                u, par, du = stack.pop()
                D[s, u] = du
                for v, l, _ in adj[u]:
                    if v != par:
                        stack.append((v, u, du + l))
        self.D = D
        # side[e][r] = True if leaf r is on the b side of edge e
        side = np.zeros((len(self.edges), n), dtype=bool)
        for i, (a, b, l) in enumerate(self.edges):
            side[i] = D[b, :n] < D[a, :n]
        self.side = side
        self.m = len(self.edges)
        base = np.zeros((self.m, n))
        sign = np.zeros((self.m, n))
        for i, (a, b, l) in enumerate(self.edges):
            bs = side[i]
            base[i] = np.where(bs, D[b, :n] + l, D[a, :n])
            sign[i] = np.where(bs, -1.0, 1.0)
        self.base, self.sign = base, sign
        self.length = np.array([l for _, _, l in self.edges])
        self.is_pendant = np.array([min(a, b) < n for a, b, _ in self.edges])

    def adjacent(self, e, f):
        a, b, _ = self.edges[e]
        c, d, _ = self.edges[f]
        return e != f and len({a, b} & {c, d}) > 0

    def system(self, S):
        """Matrix of the linear system for edge set S: columns p_1..p_k, w_1..w_k, ybar."""
        k = len(S)
        A = np.zeros((self.n + 1, 2 * k + 1))
        for j, e in enumerate(S):
            A[: self.n, j] = self.base[e]
            A[self.n, j] = 1.0
            A[: self.n, k + j] = self.sign[e]
        A[: self.n, 2 * k] = 1.0
        return A

    def mixture(self, S, p, x, ybar):
        """Noise-free d for queries on edges S at relative position x (from endpoint a)."""
        d = np.full(self.n, float(ybar))
        for e, pq, xq in zip(S, p, x):
            d += pq * (self.base[e] + self.sign[e] * xq * self.length[e])
        return d

    def newick(self):
        # root at internal node n (or leaf 0 for n<=2)
        root = self.n if self.N > self.n else 0

        def rec(u, par):
            ch = [(v, l) for v, l, _ in self.adj[u] if v != par]
            if not ch:
                return str(u)
            return "(" + ",".join(f"{rec(v, u)}:{l:g}" for v, l in ch) + ")" + ("" if u < self.n else f"i{u}")

        return rec(root, -1) + ";"


def all_topologies(n):
    """All unrooted binary topologies on leaves 0..n-1 (stepwise addition), as edge lists (a,b)."""
    out = []

    def grow(edges, nxt_leaf, nxt_int):
        if nxt_leaf == n:
            out.append(list(edges))
            return
        for i, (a, b) in enumerate(edges):
            u = nxt_int
            new = edges[:i] + edges[i + 1:] + [(a, u), (u, b), (u, nxt_leaf)]
            grow(new, nxt_leaf + 1, nxt_int + 1)

    if n == 3:
        return [[(0, 3), (1, 3), (2, 3)]]
    grow([(0, n), (1, n), (2, n)], 3, n + 1)
    return out


def strict_lp(A, rhs, S, length, tol_eq=1e-8):
    """Max slack s s.t. A z = rhs, p_q >= s, w_q >= s*l_q, p_q*l_q - w_q >= s*l_q, ybar >= 0.
    Returns (s*, z) or (None, None) if infeasible. s*>0 means a valid strictly interior solution."""
    k = len(S)
    nv = 2 * k + 2  # p, w, y, s
    c = np.zeros(nv)
    c[-1] = -1.0
    Aeq = np.hstack([A, np.zeros((A.shape[0], 1))])
    rows, b = [], []
    for j, e in enumerate(S):
        l = length[e]
        r = np.zeros(nv); r[j] = -1; r[-1] = 1; rows.append(r); b.append(0)          # s - p <= 0
        r = np.zeros(nv); r[k + j] = -1; r[-1] = l; rows.append(r); b.append(0)      # s*l - w <= 0
        r = np.zeros(nv); r[j] = -l; r[k + j] = 1; r[-1] = l; rows.append(r); b.append(0)  # w - p l + s l <= 0
    bounds = [(None, None)] * (2 * k) + [(0, None), (None, 1.0)]
    res = linprog(c, A_ub=np.array(rows), b_ub=np.array(b), A_eq=Aeq, b_eq=rhs,
                  bounds=bounds, method="highs")
    if res.status != 0:
        return None, None
    z = res.x
    if np.max(np.abs(A @ z[:-1] - rhs)) > tol_eq:
        return None, None
    return z[-1], z[:-1]


class Enumerator:
    """Precomputes residual projectors for all edge subsets of size <= K on one topology."""

    def __init__(self, tree, K):
        self.t = tree
        self.subsets = [S for k in range(1, K + 1) for S in itertools.combinations(range(tree.m), k)]
        self.P = []
        self.rank = []
        for S in self.subsets:
            A = tree.system(S)
            U, sv, _ = np.linalg.svd(A, full_matrices=False)
            r = int(np.sum(sv > 1e-9 * max(1.0, sv[0])))
            self.rank.append(r)
            Q = U[:, :r]
            self.P.append(Q)

    def solutions(self, d, Kmax=None, res_tol=1e-8, s_tol=1e-7):
        """All edge sets S' (|S'|<=Kmax) admitting a strictly valid exact solution for d."""
        rhs = np.concatenate([d, [1.0]])
        sc = max(1.0, np.abs(rhs).max())
        found = []
        for S, Q, r in zip(self.subsets, self.P, self.rank):
            if Kmax is not None and len(S) > Kmax:
                continue
            res = rhs - Q @ (Q.T @ rhs)
            if np.abs(res).max() > res_tol * sc:
                continue
            A = self.t.system(S)
            k = len(S)
            if r == 2 * k + 1:  # unique solution: check strict validity directly, no LP needed
                z = np.linalg.lstsq(A, rhs, rcond=None)[0]
                L = self.t.length[list(S)]
                p, w, y = z[:k], z[k:2 * k], z[2 * k]
                s = min(p.min(), (w / L).min(), ((p * L - w) / L).min(), 1.0)
                if y < -1e-9:
                    s = -1.0
            else:
                s, z = strict_lp(A, rhs, S, self.t.length)
            if s is not None and s > s_tol:
                found.append(dict(S=S, slack=s, z=z, rank=r, full_rank=(r == 2 * len(S) + 1)))
        return found
