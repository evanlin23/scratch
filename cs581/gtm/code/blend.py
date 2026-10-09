"""GTM-Blend: constrained SPR hill-climbing that starts from a merged tree (e.g. GTM)
and only accepts moves that keep every subset tree T_i induced (T|S_i == T_i).
Subsets may therefore interleave ("blend") in the output.

Score: Fitch parsimony on the alignment (exact incremental SPR deltas), or an
oracle score (RF to a given tree; used only to measure reachability).

Feasibility of an SPR move (prune X at p, regraft on target edge; path edges
e_1..e_m from p with far sides F_1 >= ... >= F_m): the new split set is the old one
minus {F_1..F_{m-1}} plus {F_2 u X, ..., F_m u X}. As T_i is binary, T'|S_i == T_i iff
every restriction (F_j u X) | S_i, j=2..m, is trivial or already a split of T_i.
"""
import sys
import numpy as np
from phylo import read_tree, bipartitions

ENC = np.full(256, 15, dtype=np.uint8)
for c, v in zip("ACGTU", (1, 2, 4, 8, 8)):
    ENC[ord(c)] = v
    ENC[ord(c.lower())] = v
for c, v in (("R", 5), ("Y", 10), ("S", 6), ("W", 9), ("K", 12), ("M", 3)):
    ENC[ord(c)] = v


def compress(seqs, names):
    """Parsimony-informative unique columns + weights. seqs: name->str."""
    A = np.vstack([ENC[np.frombuffer(seqs[n].encode(), dtype=np.uint8)] for n in names])
    keep = []
    for j in range(A.shape[1]):
        col = A[:, j]
        cnt = [(col == s).sum() for s in (1, 2, 4, 8)]
        if sum(c >= 2 for c in cnt) >= 2:
            keep.append(j)
    A = A[:, keep]
    u, inv, w = np.unique(A, axis=1, return_inverse=True, return_counts=True)
    return np.ascontiguousarray(u), w.astype(np.int64)


def fitch(a, b):
    i = a & b
    return np.where(i == 0, a | b, i)


def resolve(tree):
    """Copy of phylo.Tree with every node of degree > 3 resolved (caterpillar)."""
    t = tree.copy()
    for v in list(t.adj):
        while len(t.adj[v]) > 3:
            nb = list(t.adj[v])
            a, b = nb[0], nb[1]
            z = t.new()
            t.disconnect(v, a)
            t.disconnect(v, b)
            t.connect(z, a)
            t.connect(z, b)
            t.connect(v, z)
    return t


class BlendSearch:
    def __init__(self, tree, names, states, weights, subset_names, radius=6, oracle=None):
        """tree: phylo.Tree; names: taxon order (index = bit); states: (n, L) uint8;
        subset_names: list of lists of leaf names; oracle: set of true splits (optional)."""
        self.names = names
        self.idx = {n: i for i, n in enumerate(names)}
        self.n = len(names)
        self.S = states
        self.w = weights
        self.radius = radius
        self.oracle = oracle
        self.ALL = (1 << self.n) - 1
        # adjacency with leaves 0..n-1, internal n..2n-3
        tree = resolve(tree)
        old2new = {}
        for v, nm in tree.label.items():
            old2new[v] = self.idx[nm]
        k = self.n
        for v in tree.adj:
            if v not in old2new:
                old2new[v] = k
                k += 1
        self.N = k
        self.adj = [[] for _ in range(k)]
        for v, nb in tree.adj.items():
            self.adj[old2new[v]] = [old2new[u] for u in nb]
        # subsets: masks and normalized split sets (from the CURRENT tree, which must
        # already satisfy the constraints)
        self.masks = []
        for sub in subset_names:
            m = 0
            for nm in sub:
                m |= 1 << self.idx[nm]
            self.masks.append(m)
        self.recompute()
        self.splits = []
        allsp = self.all_splits()
        for m in self.masks:
            self.splits.append({r for b in allsp if (r := self.rnorm(b, m)) is not None})

    # ---------- bitset helpers
    @staticmethod
    def rnorm(b, mask):
        r = b & mask
        c = r.bit_count()
        if c <= 1 or c >= mask.bit_count() - 1:
            return None
        low = mask & -mask
        return r if not (r & low) else mask ^ r

    def far(self, u, v):
        """taxa on v's side of edge (u,v)."""
        return self.C[v] if self.par[v] == u else self.ALL ^ self.C[u]

    def D(self, u, v):
        """Fitch set of v's side of edge (u,v)."""
        return self.down[v] if self.par[v] == u else self.up[u]

    def all_splits(self):
        out = []
        for v in range(self.n, self.N):
            if self.par[v] is not None and self.par[v] >= self.n:
                b = self.C[v]
                out.append(b if not (b & 1) else self.ALL ^ b)
        return out

    # ---------- full recomputation (rooted at leaf 0)
    def recompute(self):
        N = self.N
        par = [None] * N
        order = []
        stack = [0]
        par[0] = -1
        while stack:
            v = stack.pop()
            order.append(v)
            for u in self.adj[v]:
                if par[u] is None:
                    par[u] = v
                    stack.append(u)
        par[0] = None
        self.par = par
        self.order = order
        C = [0] * N
        down = [None] * N
        length = 0
        for v in reversed(order):
            ch = [u for u in self.adj[v] if u != par[v]]
            if v < self.n:
                C[v] = 1 << v
                down[v] = self.S[v]
                if v != 0:
                    continue
            if v == 0:
                continue
            a, b = ch
            C[v] = C[a] | C[b]
            i = down[a] & down[b]
            e = i == 0
            length += int(self.w[e].sum())
            down[v] = np.where(e, down[a] | down[b], i)
        r1 = self.adj[0][0]
        length += int(self.w[(down[r1] & self.S[0]) == 0].sum())
        C[0] = self.ALL
        self.C = C
        self.down = down
        up = [None] * N
        up[r1] = self.S[0]
        for v in order[1:]:
            if v < self.n:
                continue
            ch = [u for u in self.adj[v] if u != par[v]]
            a, b = ch
            up[a] = fitch(up[v], down[b])
            up[b] = fitch(up[v], down[a])
        self.up = up
        self.length = length

    def rf_to_oracle(self):
        return len(self.oracle - set(self.all_splits()))

    # ---------- moves
    def feasible(self, X, fars):
        for m, sp in zip(self.masks, self.splits):
            if not (X & m) or (m & ~X) == 0:
                continue
            for F in fars:
                r = self.rnorm(F | X, m)
                if r is not None and r not in sp:
                    return False
        return True

    def candidates(self):
        """Yield (delta, p, x, path) for improving parsimony moves.
        path = list of nodes [p, a, ..., c, d]; target edge (c, d)."""
        w = self.w
        out = []
        for p in range(self.n, self.N):
            nb = self.adj[p]
            for x in nb:
                DX = self.D(p, x)
                a_, b_ = [u for u in nb if u != x]
                base = None
                for s, o in ((a_, b_), (b_, a_)):
                    Ds, Do = self.D(p, s), self.D(p, o)
                    if base is None:
                        base = int(w[(fitch(Ds, Do) & DX) == 0].sum())
                    # walk from s away from p; N_in = set of near side (o side)
                    stack = [(p, s, Do, 1, [p, s])]
                    while stack:
                        uprev, u, Nin, depth, path = stack.pop()
                        if u < self.n or depth > self.radius:
                            continue
                        others = [v for v in self.adj[u] if v != uprev]
                        for t in others:
                            t2 = others[0] if others[1] == t else others[1]
                            near = fitch(Nin, self.D(u, t2))
                            c = int(w[(fitch(near, self.D(u, t)) & DX) == 0].sum())
                            if c < base:
                                out.append((c - base, p, x, path + [t]))
                            stack.append((u, t, near, depth + 1, path + [t]))
        out.sort(key=lambda z: z[0])
        return out

    def oracle_candidates(self):
        """All feasible-or-not moves within radius scored by RF to the oracle: here we
        just enumerate and let the caller score (slow; used on small budgets)."""
        out = []
        for p in range(self.n, self.N):
            nb = self.adj[p]
            for x in nb:
                a_, b_ = [u for u in nb if u != x]
                for s in (a_, b_):
                    stack = [(p, s, 1, [p, s])]
                    while stack:
                        uprev, u, depth, path = stack.pop()
                        if u < self.n or depth > self.radius:
                            continue
                        for t in self.adj[u]:
                            if t == uprev:
                                continue
                            out.append((p, x, path + [t]))
                            stack.append((u, t, depth + 1, path + [t]))
        return out

    def path_fars(self, path):
        # far sides of e_2..e_m, i.e. edges (path[k], path[k+1]) for k>=1
        return [self.far(path[k], path[k + 1]) for k in range(1, len(path) - 1)]

    def apply(self, p, x, path):
        a_, b_ = [u for u in self.adj[p] if u != x]
        c, d = path[-2], path[-1]
        # remove p
        self.adj[p] = [x]
        self.adj[a_].remove(p)
        self.adj[b_].remove(p)
        self.adj[a_].append(b_)
        self.adj[b_].append(a_)
        # regraft on (c, d)
        self.adj[c].remove(d)
        self.adj[d].remove(c)
        self.adj[c].append(p)
        self.adj[d].append(p)
        self.adj[p] += [c, d]

    def new_splits_if(self, p, x, path):
        """Split set after the move (for oracle scoring)."""
        X = self.far(p, x)
        fars = [self.far(path[k], path[k + 1]) for k in range(0, len(path) - 1)]
        cur = set(self.all_splits())
        norm = lambda b: b if not (b & 1) else self.ALL ^ b
        for F in fars[:-1]:
            cur.discard(norm(F))
        for F in fars[1:]:
            cur.add(norm(F | X))
        return cur

    def run_parsimony(self, max_moves=10 ** 6, log=None):
        moves = 0
        while moves < max_moves:
            done = False
            for delta, p, x, path in self.candidates():
                X = self.far(p, x)
                if self.feasible(X, self.path_fars(path)):
                    before = self.length
                    self.apply(p, x, path)
                    self.recompute()
                    assert self.length - before == delta, (self.length, before, delta)
                    moves += 1
                    done = True
                    if log:
                        log(moves, self.length, delta)
                    break
            if not done:
                break
        return moves

    def to_newick(self):
        def rec(v, p):
            if v < self.n:
                return self.names[v]
            return "(" + ",".join(rec(u, v) for u in self.adj[v] if u != p) + ")"
        r = self.adj[0][0]
        return "(" + ",".join([self.names[0]] + [rec(u, r) for u in self.adj[r] if u != 0]) + ");"


def run_oracle(bs, max_moves=10 ** 6, log=None):
    """Hill-climb on #true splits recovered (oracle; measures reachability only)."""
    norm = lambda b: b if not (b & 1) else bs.ALL ^ b
    moves = 0
    while moves < max_moves:
        cur = set(bs.all_splits())
        best = []
        for p, x, path in bs.oracle_candidates():
            X = bs.far(p, x)
            fars = [bs.far(path[k], path[k + 1]) for k in range(len(path) - 1)]
            removed = {norm(F) for F in fars[:-1]}
            added = {norm(F | X) for F in fars[1:]}
            after = (cur - removed) | added
            d = len(bs.oracle & cur) - len(bs.oracle & after)  # negative = better
            if d < 0:
                best.append((d, p, x, path))
        best.sort(key=lambda z: z[0])
        done = False
        for d, p, x, path in best:
            if bs.feasible(bs.far(p, x), bs.path_fars(path)):
                bs.apply(p, x, path)
                bs.recompute()
                moves += 1
                done = True
                if log:
                    log(moves, bs.rf_to_oracle(), d)
                break
        if not done:
            break
    return moves
