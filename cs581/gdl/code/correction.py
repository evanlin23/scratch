"""Survival-reweighted ASTRID-Pro ("ASTRID-Pro-S"): a candidate correction for supercritical
branches.

Theory (results/theory.md). For an orthologous pair the expected speciation-only distance is
1 + sum_v s(off(v)). If every visible speciation node v is counted with weight 1/s(off(v)),
the expectation becomes 1 + (#internal species nodes strictly inside the path), i.e. the
topological distance in the species tree, which is additive.

s(c) is estimated by mapping gene trees onto a rooted first-pass species tree S.
For every maximal gene subtree X whose species lie inside clade(w) (and whose parent's do not),
the copy that entered branch w from u = parent(w) survived on the w side. Its sibling daughter
(entering branch sib(w)) survived iff X's parent is a speciation node mapped to u.
So  s_hat(sib(w)) = #both / #observed.  Independence of the two daughters makes this unbiased
when S is correct.
"""
import numpy as np

from phylo import Tree, parse_newick


class SpTree:
    """Rooted species tree with clade bitmasks."""

    def __init__(self, t, sp_index):
        self.t = t
        n = len(t.parent)
        self.mask = [0] * n
        for v in t.postorder():
            if t.children[v]:
                m = 0
                for c in t.children[v]:
                    m |= self.mask[c]
                self.mask[v] = m
            else:
                self.mask[v] = 1 << sp_index[t.label[v]]
        # nodes by increasing clade size, for LCA lookup
        self.order = sorted(range(n), key=lambda v: bin(self.mask[v]).count("1"))
        self.sib = {}
        for v in range(n):
            ch = t.children[v]
            if len(ch) == 2:
                self.sib[ch[0]], self.sib[ch[1]] = ch[1], ch[0]
        self._cache = {}

    def lca(self, m):
        r = self._cache.get(m)
        if r is None:
            for v in self.order:
                if self.mask[v] & m == m:
                    r = v
                    break
            self._cache[m] = r
        return r


def reroot_species_tree(unrooted, outgroup_species):
    """Root an unrooted Newick tree on the branch above the given leaf (simple outgroup rooting)."""
    t = unrooted
    nb = [[] for _ in t.parent]
    for v, p in enumerate(t.parent):
        if p >= 0:
            nb[v].append(p); nb[p].append(v)
    leaf = [v for v in t.leaves() if t.label[v] == outgroup_species][0]
    other = nb[leaf][0]
    rt = Tree()
    root = rt.add(-1)
    rt.root = root
    stack = [(leaf, other, root), (other, leaf, root)]
    while stack:
        v, frm, p = stack.pop()
        i = rt.add(p)
        rt.label[i] = t.label[v] if not [q for q in nb[v] if q != frm] else None
        for q in nb[v]:
            if q != frm:
                stack.append((q, v, i))
    from phylo import compact
    return compact(rt)


def estimate_survival(tagged, S):
    """tagged: list of (rt, tags, leafsp) with rt._S species bitmasks. Returns dict c -> s_hat."""
    both, tot = {}, {}
    sroot = S.t.root
    for rt, tags, _ in tagged:
        sm = rt._S
        for x in range(len(rt.parent)):
            p = rt.parent[x]
            mx = S.lca(sm[x])
            mp = S.lca(sm[p]) if p >= 0 else None
            # walk w from mx up to (excluding) mp (or the species root if x is the gene root)
            w = mx
            stop = mp if p >= 0 else sroot
            while w is not None and w != stop and w != sroot:
                u = S.t.parent[w]
                if p >= 0 and u == stop:
                    # last step: parent gene node is mapped to u
                    is_both = tags[p] is False and S.lca(sm[p]) == u
                else:
                    is_both = False
                if p >= 0 and S.mask[w] & sm[p] == sm[p]:
                    break  # parent still inside clade(w): x is not maximal for w
                s = S.sib.get(w)
                if s is not None:
                    tot[s] = tot.get(s, 0) + 1
                    if is_both:
                        both[s] = both.get(s, 0) + 1
                w = u
    return {c: both.get(c, 0) / tot[c] for c in tot if tot[c] > 0}


def corrected_distances(rt, tags, leafsp, nsp, S, shat, floor=0.02):
    """Orthologous pairs only; LCA counts 1; each other visible speciation node v on the path
    counts 1/s_hat(off-path child clade of map(v))."""
    sm = rt._S
    n = len(rt.parent)
    E = [0.0] * n
    for v in rt.postorder()[::-1]:
        p = rt.parent[v]
        if p < 0:
            continue
        add = 0.0
        if not tags[p] and len(rt.children[p]) == 2:
            off = [c for c in rt.children[p] if c != v][0]
            u = S.lca(sm[p])
            # S child of u containing the off side
            coff = None
            for c in S.t.children[u]:
                if S.mask[c] & sm[off] == sm[off]:
                    coff = c
            sh = shat.get(coff) if coff is not None else None
            add = 1.0 / max(sh, floor) if sh else 1.0
        E[v] = E[p] + add
    tot = np.zeros((nsp, nsp)); cnt = np.zeros((nsp, nsp))
    leaves_under = [None] * n
    for v in rt.postorder():
        ch = rt.children[v]
        if not ch:
            leaves_under[v] = np.array([v]); continue
        L = [leaves_under[c] for c in ch]
        leaves_under[v] = np.concatenate(L)
        if tags[v]:
            continue
        for i in range(len(ch)):
            for j in range(i + 1, len(ch)):
                a, b = L[i], L[j]
                ea = np.array([E[x] for x in a]) - E[ch[i]]
                eb = np.array([E[x] for x in b]) - E[ch[j]]
                d = ea[:, None] + eb[None, :] + 1.0
                sa = np.array([leafsp[x] for x in a]); sb = np.array([leafsp[x] for x in b])
                SA = np.broadcast_to(sa[:, None], d.shape).ravel(); SB = np.broadcast_to(sb[None, :], d.shape).ravel()
                dd = d.ravel(); keep = SA != SB
                np.add.at(tot, (SA[keep], SB[keep]), dd[keep]); np.add.at(cnt, (SA[keep], SB[keep]), 1)
    tot = tot + tot.T; cnt = cnt + cnt.T
    has = cnt > 0
    out = np.zeros((nsp, nsp)); out[has] = tot[has] / cnt[has]
    return out, has.astype(float)
