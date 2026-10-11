"""Constrained insertion merger (blending DTM).

Start with the largest subset tree; insert every other taxon one at a time on an
edge that keeps T*|S_i == T_i|(inserted part of S_i) for its own subset i. Among the
feasible edges choose by
  - 'pars'  : Fitch insertion cost on the alignment (exact length increase), ties
              broken by the guide tree (fewest new conflicts with guide splits);
  - 'guide' : best match to the guide tree's attachment split (data-free);
  - 'oracle': best match to the TRUE tree's attachment split (measures headroom).
Subsets therefore interleave freely ("blending").
"""
import numpy as np
from phylo import read_tree


def fitch(a, b):
    i = a & b
    return np.where(i == 0, a | b, i)


class Grow:
    def __init__(self, names, states, weights):
        self.names = names
        self.idx = {n: i for i, n in enumerate(names)}
        self.n = len(names)
        self.S = states
        self.w = weights
        self.adj = {}
        self.next = self.n  # internal node ids

    def init_from(self, tree):
        """tree: phylo.Tree on a subset of names."""
        m = {}
        for v in tree.adj:
            if v in tree.label:
                m[v] = self.idx[tree.label[v]]
            else:
                m[v] = self.next
                self.next += 1
        for v, nb in tree.adj.items():
            self.adj[m[v]] = [m[u] for u in nb]

    def leaves(self):
        return [v for v in self.adj if v < self.n]

    def recompute(self, with_fitch=True):
        root = min(self.leaves())
        par = {root: None}
        order = []
        stack = [root]
        while stack:
            v = stack.pop()
            order.append(v)
            for u in self.adj[v]:
                if u not in par:
                    par[u] = v
                    stack.append(u)
        self.par, self.order, self.root = par, order, root
        C = {}
        for v in reversed(order):
            if v < self.n:
                C[v] = 1 << v
            else:
                c = 0
                for u in self.adj[v]:
                    if u != par[v]:
                        c |= C[u]
                C[v] = c
        allm = 0
        for v in self.leaves():
            allm |= 1 << v
        C[root] = allm
        self.C, self.ALL = C, allm
        if not with_fitch:
            return
        down = {}
        for v in reversed(order):
            if v < self.n:
                down[v] = self.S[v]
            if v >= self.n:
                a, b = [u for u in self.adj[v] if u != par[v]]
                down[v] = fitch(down[a], down[b])
        up = {}
        r1 = self.adj[root][0]
        up[r1] = self.S[root]
        for v in order:
            if v >= self.n:
                a, b = [u for u in self.adj[v] if u != par[v]]
                up[a] = fitch(up[v], down[b])
                up[b] = fitch(up[v], down[a])
        self.down, self.up = down, up

    def edges(self):
        """(parent, child) pairs; far side of child = C[child]."""
        return [(self.par[v], v) for v in self.order if self.par[v] is not None]

    def insert(self, leaf, p, c):
        z = self.next
        self.next += 1
        self.adj[p].remove(c)
        self.adj[c].remove(p)
        self.adj[p].append(z)
        self.adj[c].append(z)
        self.adj[z] = [p, c, leaf]
        self.adj[leaf] = [z]

    def to_newick(self):
        def rec(v, p):
            if v < self.n:
                return self.names[v]
            return "(" + ",".join(rec(u, v) for u in self.adj[v] if u != p) + ")"
        r = self.root
        z = self.adj[r][0]
        return "(" + ",".join([self.names[r]] + [rec(u, z) for u in self.adj[z] if u != r]) + ");"


def norm(b, mask):
    low = mask & -mask
    return b if not (b & low) else mask ^ b


def sibling_sides(tree, leafname, keep, idx):
    """Restrict phylo.Tree `tree` to keep+{leafname}; return (Y1, Y2) bitsets of the
    two clades adjacent to leafname (Y1|Y2 partition keep)."""
    t = tree.copy().restrict(set(keep) | {leafname})
    lv = [v for v in t.label if t.label[v] == leafname][0]
    z = next(iter(t.adj[lv]))
    out = []
    for u in t.adj[z]:
        if u == lv:
            continue
        b = 0
        stack = [(u, z)]
        while stack:
            v, p = stack.pop()
            if v in t.label:
                b |= 1 << idx[t.label[v]]
            for y in t.adj[v]:
                if y != p:
                    stack.append((y, v))
        out.append(b)
    return out


def feasible_edges(g, Q, tau):
    """Edges (p,c) of g where a new leaf can go so that its subset's induced tree has
    the leaf on the T_i edge with Q-split tau (bitset, normalized within Q)."""
    E = g.edges()
    if Q.bit_count() <= 2:
        return E
    path = [(p, c) for p, c in E if norm(g.C[c] & Q, Q) == tau]
    if not path:
        raise RuntimeError("no path edges")
    pe = set(path)
    cnt = {}
    for p, c in path:
        cnt[p] = cnt.get(p, 0) + 1
        cnt[c] = cnt.get(c, 0) + 1
    out = list(path)
    for v, k in cnt.items():
        if k == 2:  # interior node of the path: add Q-free hanging subtree
            for u in g.adj[v]:
                if (v, u) in pe or (u, v) in pe:
                    continue
                stack = [(u, v)]
                while stack:
                    x, px = stack.pop()
                    out.append((px, x) if g.par.get(x) == px else (x, px))
                    for y in g.adj[x]:
                        if y != px:
                            stack.append((y, x))
    return out


def run(names, states, weights, subsets, mode, ref_tree=None, order_key=None, log=None):
    """subsets: list of phylo.Tree. ref_tree: guide (mode guide/pars tiebreak) or
    true tree (mode oracle)."""
    g = Grow(names, states, weights)
    idx = g.idx
    subsets = sorted(subsets, key=lambda t: -len(t.leaves()))
    g.init_from(subsets[0])
    done = {g.names[v] for v in g.leaves()}
    sub_of = {}
    for i, t in enumerate(subsets):
        for v in t.leaves():
            sub_of[t.label[v]] = i
    rest = [n for n in names if n not in done]
    if order_key:
        rest.sort(key=order_key)
    inserted = {i: [] for i in range(len(subsets))}
    for ln in rest:
        i = sub_of[ln]
        g.recompute(with_fitch=(mode == "pars"))
        Qn = inserted[i]
        Q = 0
        for q in Qn:
            Q |= 1 << idx[q]
        if len(Qn) >= 3:
            y1, y2 = sibling_sides(subsets[i], ln, Qn, idx)
            tau = norm(y1, Q)
        else:
            tau = None
        cand = feasible_edges(g, Q, tau)
        leaf = idx[ln]
        if True:
            if ref_tree is not None:
                P = list(done)
                sides = sibling_sides(ref_tree, ln, P, idx)  # >2 if ref has a polytomy

                def gscore(pc):
                    F = g.C[pc[1]]
                    return min((F ^ r).bit_count() for r in sides)
        if mode == "pars":
            sl = g.S[leaf]
            best = None
            for p, c in cand:
                dn = g.down[c]
                upc = g.up[c] if c in g.up else None
                if upc is None:
                    continue
                cost = int(g.w[(fitch(dn, upc) & sl) == 0].sum())
                key = (cost, gscore((p, c)) if ref_tree is not None else 0)
                if best is None or key < best[0]:
                    best = (key, p, c)
            _, p, c = best
        else:
            p, c = min(cand, key=gscore)
        g.insert(leaf, p, c)
        done.add(ln)
        inserted[i].append(ln)
        if log and len(done) % 100 == 0:
            log(len(done))
    g.recompute(with_fitch=False)
    return g
