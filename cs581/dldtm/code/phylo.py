"""Minimal unrooted-tree utilities: Newick I/O, bipartitions as int bitsets,
FN/FP/RF error (Park et al. 2021 criterion), induced-subtree checks.

A tree is stored as an undirected adjacency dict {node_id: set(node_id)} plus
leaf labels {node_id: name}. Leaves have degree 1.
"""
import re
import sys

sys.setrecursionlimit(100000)

_tok = re.compile(r"\s*('(?:[^']|'')*'|[(),;:]|[^(),;:\s]+)")


class Tree:
    def __init__(self):
        self.adj = {}
        self.label = {}   # leaf id -> name
        self.n = 0

    def new(self):
        i = self.n
        self.n += 1
        self.adj[i] = set()
        return i

    def connect(self, a, b):
        self.adj[a].add(b)
        self.adj[b].add(a)

    def disconnect(self, a, b):
        self.adj[a].discard(b)
        self.adj[b].discard(a)

    def leaves(self):
        return [v for v in self.adj if len(self.adj[v]) == 1 and v in self.label]

    def copy(self):
        t = Tree()
        t.adj = {k: set(v) for k, v in self.adj.items()}
        t.label = dict(self.label)
        t.n = self.n
        return t

    def suppress_deg2(self):
        for v in list(self.adj):
            if v in self.adj and len(self.adj[v]) == 2 and v not in self.label:
                a, b = self.adj[v]
                self.disconnect(v, a)
                self.disconnect(v, b)
                self.connect(a, b)
                del self.adj[v]

    def restrict(self, keep):
        """Induced subtree on leaf names in `keep` (in place)."""
        keep = set(keep)
        for v, nm in list(self.label.items()):
            if nm not in keep:
                for u in list(self.adj[v]):
                    self.disconnect(v, u)
                del self.adj[v]
                del self.label[v]
        # prune dangling internal nodes
        stack = [v for v in self.adj if len(self.adj[v]) <= 1 and v not in self.label]
        while stack:
            v = stack.pop()
            if v not in self.adj:
                continue
            for u in list(self.adj[v]):
                self.disconnect(v, u)
                if len(self.adj[u]) <= 1 and u not in self.label:
                    stack.append(u)
            del self.adj[v]
        self.suppress_deg2()
        return self

    def to_newick(self):
        lv = self.leaves()
        # root at an internal node next to the first leaf
        root = next(iter(self.adj[lv[0]])) if len(self.adj) > 2 else lv[0]
        out = []

        def rec(v, p):
            ch = [u for u in self.adj[v] if u != p]
            if not ch:
                out.append(self.label[v])
                return
            out.append("(")
            for i, u in enumerate(ch):
                if i:
                    out.append(",")
                rec(u, v)
            out.append(")")
        # iterative-safe recursion is fine for n<=1e4 with raised limit
        rec(root, None)
        return "".join(out) + ";"


def parse_newick(s):
    t = Tree()
    stack = []
    last = None   # node just closed by ')', so a following name is its label
    s = re.sub(r"\[[^\]]*\]", "", s)  # strip [&R] and other comments
    toks = [m.group(1) for m in _tok.finditer(s.strip())]
    i = 0
    while i < len(toks):
        x = toks[i]
        if x == "(":
            c = t.new()
            if stack:
                t.connect(stack[-1], c)
            stack.append(c)
            last = None
        elif x == ",":
            last = None
        elif x == ")":
            last = stack.pop()
        elif x == ":":
            i += 1  # skip branch length
        elif x == ";":
            break
        elif last is None:
            nm = x[1:-1].replace("''", "'") if x.startswith("'") else x
            c = t.new()
            t.connect(stack[-1], c)
            t.label[c] = nm
        # else: internal node label / support value, ignored
        i += 1
    t.suppress_deg2()
    return t


def read_tree(path):
    with open(path) as f:
        return parse_newick(f.read().strip().split(";")[0] + ";")


def bipartitions(t, index):
    """Set of nontrivial bipartitions as int bitsets, normalized to exclude
    the taxon with index 0. `index` maps leaf name -> bit."""
    lv = t.leaves()
    full = 0
    for v in lv:
        full |= 1 << index[t.label[v]]
    root = lv[0]
    bits = {}
    order = []
    parent = {root: None}
    stack = [root]
    while stack:
        v = stack.pop()
        order.append(v)
        for u in t.adj[v]:
            if u != parent[v]:
                parent[u] = v
                stack.append(u)
    out = set()
    nl = len(lv)
    for v in reversed(order):
        b = (1 << index[t.label[v]]) if v in t.label else 0
        for u in t.adj[v]:
            if u != parent[v]:
                b |= bits[u]
        bits[v] = b
        if v != root and v not in t.label:
            c = bin(b).count("1")
            if 1 < c < nl - 1:
                out.add(b if not (b & 1) else full ^ b)
    return out


def fn_fp(est, true):
    """Return (FN rate, FP rate, nFN, nFP, n_int_true, n_int_est) on the
    common leaf set."""
    a = {est.label[v] for v in est.leaves()}
    b = {true.label[v] for v in true.leaves()}
    common = sorted(a & b)
    if a != b:
        est = est.copy().restrict(common)
        true = true.copy().restrict(common)
    idx = {n: i for i, n in enumerate(common)}
    be = bipartitions(est, idx)
    bt = bipartitions(true, idx)
    fn = len(bt - be)
    fp = len(be - bt)
    return fn / len(bt), (fp / len(be) if be else 0.0), fn, fp, len(bt), len(be)


def is_induced(big, small):
    """True iff big restricted to leaves(small) equals small (topologically)."""
    names = [small.label[v] for v in small.leaves()]
    r = big.copy().restrict(names)
    idx = {n: i for i, n in enumerate(sorted(names))}
    return bipartitions(r, idx) == bipartitions(small, idx)


if __name__ == "__main__":
    e, tr = read_tree(sys.argv[1]), read_tree(sys.argv[2])
    print("FN=%.4f FP=%.4f nFN=%d nFP=%d intT=%d intE=%d" % fn_fp(e, tr))
