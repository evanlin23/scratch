"""Minimal tree utilities: Newick I/O, bipartitions, RF error.

Trees are stored as flat arrays (parent, children, label) so that the distance
code in methods.py can work on integer node ids.
"""
import re

_TOK = re.compile(r"\s*([(),;:]|[^(),;:\s]+)")


class Tree:
    __slots__ = ("parent", "children", "label", "length", "root", "_S")

    def __init__(self):
        self.parent, self.children, self.label, self.length = [], [], [], []
        self.root = 0
        self._S = None

    def add(self, par):
        i = len(self.parent)
        self.parent.append(par)
        self.children.append([])
        self.label.append(None)
        self.length.append(None)
        if par >= 0:
            self.children[par].append(i)
        return i

    def is_leaf(self, v):
        return not self.children[v]

    def leaves(self):
        return [v for v in range(len(self.parent)) if not self.children[v]]

    def postorder(self, v=None):
        v = self.root if v is None else v
        out, stack = [], [(v, False)]
        while stack:
            u, done = stack.pop()
            if done:
                out.append(u)
            else:
                stack.append((u, True))
                for c in self.children[u]:
                    stack.append((c, False))
        return out

    def newick(self, v=None, lengths=False):
        v = self.root if v is None else v
        def rec(u):
            if not self.children[u]:
                s = self.label[u]
            else:
                s = "(" + ",".join(rec(c) for c in self.children[u]) + ")"
            if lengths and self.length[u] is not None and u != v:
                s += ":%g" % self.length[u]
            return s
        return rec(v) + ";"


def parse_newick(s):
    t = Tree()
    cur = t.add(-1)
    t.root = cur
    expect_len = False
    for tok in _TOK.findall(s):
        if tok == "(":
            cur = t.add(cur)
        elif tok == ",":
            cur = t.add(t.parent[cur])
        elif tok == ")":
            cur = t.parent[cur]
        elif tok == ":":
            expect_len = True
            continue
        elif tok == ";":
            break
        else:
            if expect_len:
                try:
                    t.length[cur] = float(tok)
                except ValueError:
                    pass
            else:
                t.label[cur] = tok.strip("'\"")
        expect_len = False
    return compact(t)


def compact(t):
    """Suppress unifurcations (degree-2 non-root nodes) and renumber."""
    new = Tree()
    def rec(u, par):
        while len(t.children[u]) == 1:
            u = t.children[u][0]
        i = new.add(par)
        new.label[i] = t.label[u]
        new.length[i] = t.length[u]
        for c in t.children[u]:
            rec(c, i)
        return i
    import sys
    sys.setrecursionlimit(100000)
    new.root = rec(t.root, -1)
    return new


def read_trees(path):
    with open(path) as f:
        return [parse_newick(line) for line in f if ";" in line]


def bipartitions(t, taxa_index, leaf_species=None):
    """Nontrivial bipartitions of an unrooted single-copy tree as bitmask ints
    normalised so the mask never contains taxon 0."""
    n = len(taxa_index)
    full = (1 << n) - 1
    mask = [0] * len(t.parent)
    for v in t.postorder():
        if t.children[v]:
            m = 0
            for c in t.children[v]:
                m |= mask[c]
            mask[v] = m
        else:
            lab = t.label[v] if leaf_species is None else leaf_species(t.label[v])
            mask[v] = 1 << taxa_index[lab]
    out = set()
    for v in range(len(t.parent)):
        m = mask[v]
        if m & 1:
            m = full ^ m
        k = bin(m).count("1")
        if 2 <= k <= n - 2:
            out.add(m)
    return out


def rf_error(est, true):
    """Returns (FN, FP, nI_true, nI_est) on the shared leaf set (species trees)."""
    taxa = sorted(set(est.label[v] for v in est.leaves()))
    tt = sorted(set(true.label[v] for v in true.leaves()))
    assert taxa == tt, "leaf sets differ"
    idx = {x: i for i, x in enumerate(taxa)}
    be, bt = bipartitions(est, idx), bipartitions(true, idx)
    return len(bt - be), len(be - bt), len(bt), len(be)
