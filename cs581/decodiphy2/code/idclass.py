"""Identifiability-aware post-processing of PDD / DecoDiPhy fits.

Given a reference tree R and a fit (edges, x, p, ybar), compute
  * the node measure m (Lemma 1),
  * flags: boundary placement (x ~ 0 or 1: the point is a node, so the edge is ambiguous; Claim 2 of
    the paper), adjacency (two placements share a node: (p, x) continuum, Thm C1), and the
    Laplacian-kernel test (B): L[V \\ V(S), I] rank-deficient (closed claw, split claw, ...),
  * the equivalence class on the fitted edge set (all (p', x', ybar') on the same edges that give the
    SAME fitted distance vector, with x' in [0, 1]):
        m' = m - L alpha,  m'_u = sum_q s_{u,q} (u in V(S)),  (L alpha)_u = 0 (u notin V(S)),
        ybar' = ybar - sum_v (deg v - 2) alpha_v >= 0,  s >= 0;   p'_q = s_{a_q,q} + s_{b_q,q}
    via 2k LPs (min / max of each p'_q), and a canonical representative (mean of those 2k vertex
    solutions, a central point of the class).
All members of the class have the same loss as the fit, so a solver's choice among them is arbitrary.

Trees: DecoDiPhy's rooted newick (degree-2 root) is converted to an unrooted tree by suppressing
the root; placements are re-expressed on the merged edge.
"""
import collections
import numpy as np
import scipy.sparse as sp
from scipy.optimize import linprog
import treeswift


class URTree:
    def __init__(self, newick):
        T = treeswift.read_tree_newick(newick)
        self.T = T
        lab = {}
        for i, nd in enumerate(T.traverse_preorder()):
            lab[nd] = i
        self.label_of = {nd.label: lab[nd] for nd in T.traverse_preorder() if nd.label is not None}
        adj = collections.defaultdict(dict)
        parent = {}
        for nd in T.traverse_preorder():
            if not nd.is_root():
                a, b = lab[nd], lab[nd.parent]
                adj[a][b] = adj[b][a] = float(nd.edge_length)
                parent[a] = b
        # suppress degree-2 nodes (the root); remember how old edges map to merged ones
        self.edge_map = {}  # (child, parent) old -> (u, v, offset_from_u, sign) merged
        for a, b in parent.items():
            self.edge_map[(a, b)] = None
        deg2 = [v for v in adj if len(adj[v]) == 2]
        merged = {}
        for v in deg2:
            (x, lx), (y, ly) = adj[v].items()
            del adj[x][v], adj[y][v], adj[v]
            adj[x][y] = adj[y][x] = lx + ly
            merged[v] = (x, lx, y, ly)
        self.merged = merged
        self.ids = sorted(adj)
        idx = {v: i for i, v in enumerate(self.ids)}
        self.idx = idx
        self.N = len(self.ids)
        self.adj = {idx[v]: {idx[w]: l for w, l in adj[v].items()} for v in adj}
        self.deg = np.array([len(self.adj[i]) for i in range(self.N)])
        self.internal = np.where(self.deg >= 3)[0]
        self.leaves = np.where(self.deg == 1)[0]
        rows, cols, vals = [], [], []
        for u in range(self.N):
            for w, l in self.adj[u].items():
                rows += [u, u]; cols += [u, w]; vals += [1 / l, -1 / l]
        self.L = sp.csr_matrix((vals, (rows, cols)), shape=(self.N, self.N))
        self.parent_old = parent
        self.lab = lab

    def placement(self, anchor_label, x):
        """anchor = child node of the placement edge (DecoDiPhy convention), x = relative distance from
        the anchor. Returns (u, w, t, L): merged edge endpoints (indices), distance t from u, length L."""
        a = self.label_of[str(anchor_label)] if str(anchor_label) in self.label_of else self.label_of[anchor_label]
        b = self.parent_old[a]
        if b in self.merged:
            x1, l1, y1, l2 = self.merged[b]
            other, lo = (y1, l2) if x1 == a else (x1, l1)
            la = l1 if x1 == a else l2
            return self.idx[a], self.idx[other], x * la, la + lo
        L = self.adj[self.idx[a]][self.idx[b]]
        return self.idx[a], self.idx[b], x * L, L

    def node_measure(self, pls, p):
        m = np.zeros(self.N)
        for (u, w, t, L), pq in zip(pls, p):
            m[u] += pq * (1 - t / L)
            m[w] += pq * t / L
        return m

    # ---------------- flags ----------------
    def flags(self, pls, x, tol_lo=1e-4, tol_hi=1e-3):
        W = set()
        adjacent = False
        for u, w, _, _ in pls:
            if u in W or w in W:
                adjacent = True
            W |= {u, w}
        boundary = any(xx < tol_lo or xx > 1 - tol_hi for xx in x)
        outside = [u for u in range(self.N) if u not in W]
        Lsub = self.L[outside][:, self.internal].toarray()
        kernel = (np.linalg.matrix_rank(Lsub, tol=1e-9) < len(self.internal)) if outside else True
        # closed claw: internal v with v and all neighbours in W
        claw = any(v in W and set(self.adj[v]) <= W for v in self.internal)
        return dict(adjacent=adjacent, boundary=boundary, kernel=bool(kernel), closed_claw=claw,
                    W=W)

    # ---------------- equivalence class on the fitted edges ----------------
    def eq_class(self, pls, p, ybar):
        k = len(pls)
        m = self.node_measure(pls, p)
        W = sorted({u for u, w, _, _ in pls} | {w for u, w, _, _ in pls})
        Wset = set(W)
        I = self.internal
        nI = len(I)
        LI = self.L[:, I]
        # variables: alpha (nI), s (2k: s_u,q then s_w,q)
        nv = nI + 2 * k
        Aeq_rows = []
        beq = []
        for u in range(self.N):
            row = sp.lil_matrix((1, nv))
            # m_u - (L alpha)_u - sum s = 0  ->  (L alpha)_u + sum s = m_u
            Lrow = LI[u]
            if Lrow.nnz == 0 and u not in Wset:
                continue
            row[0, :nI] = Lrow.toarray()
            if u in Wset:
                for q, (a, w, _, _) in enumerate(pls):
                    if a == u:
                        row[0, nI + 2 * q] = 1
                    if w == u:
                        row[0, nI + 2 * q + 1] = 1
                beq.append(m[u])
            else:
                beq.append(0.0)
            Aeq_rows.append(row)
        Aeq = sp.vstack(Aeq_rows).tocsr()
        Aub = sp.lil_matrix((1, nv)); Aub[0, :nI] = self.deg[I] - 2
        bub = [ybar]
        bounds = [(None, None)] * nI + [(0, None)] * (2 * k)
        lo, hi, sols = np.zeros(k), np.zeros(k), []
        for q in range(k):
            for sgn in (1, -1):
                c = np.zeros(nv); c[nI + 2 * q] = sgn; c[nI + 2 * q + 1] = sgn
                r = linprog(c, A_ub=Aub.tocsr(), b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=bounds, method="highs")
                if r.status != 0:
                    return None
                val = sgn * r.fun
                if sgn == 1:
                    lo[q] = val
                else:
                    hi[q] = val
                sols.append(r.x)
        Z = np.mean(sols, axis=0)
        s = Z[nI:].reshape(k, 2)
        alpha = Z[:nI]
        pc = s.sum(1)
        canon = []
        for q, (u, w, _, L) in enumerate(pls):
            t = L * s[q, 1] / pc[q] if pc[q] > 0 else 0.0
            canon.append((u, w, t, L))
        ybar_c = ybar - float((self.deg[I] - 2) @ alpha)
        verts = []
        for z in sols:
            sq = z[nI:].reshape(k, 2)
            pq = sq.sum(1)
            verts.append((pq, [(u, w, (L * sq[q, 1] / pq[q] if pq[q] > 0 else 0.0), L) for q, (u, w, _, L) in enumerate(pls)]))
        return dict(lo=lo, hi=hi, canon_p=pc, canon_pls=canon, canon_y=ybar_c, vertices=verts)

    # ---------------- earth mover's distance between point measures ----------------
    def emd(self, A, B):
        """A, B: lists of ((u, w, t, L), mass). Tree EMD (weighted UniFrac without normalisation by
        tree length), pendant lengths ignored."""
        root = int(self.internal[0]) if len(self.internal) else 0
        order, par = [], {root: -1}
        stack = [root]
        while stack:
            u = stack.pop(); order.append(u)
            for w in self.adj[u]:
                if w not in par:
                    par[w] = u; stack.append(w)
        on_edge = collections.defaultdict(list)  # child -> [(dist from child, signed mass)]
        for meas, sg in ((A, 1.0), (B, -1.0)):
            for (u, w, t, L), mass in meas:
                if par.get(u) == w:
                    on_edge[u].append((t, sg * mass))
                else:  # w is the child
                    on_edge[w].append((L - t, sg * mass))
        below = collections.defaultdict(float)
        total = 0.0
        for u in reversed(order):
            if par[u] == -1:
                continue
            L = self.adj[u][par[u]]
            cur = below[u]
            pos = 0.0
            for t, mass in sorted(on_edge[u]):
                total += abs(cur) * (t - pos); pos = t; cur += mass
            total += abs(cur) * (L - pos)
            below[par[u]] += cur
        return total
