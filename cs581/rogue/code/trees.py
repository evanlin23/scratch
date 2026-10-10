"""Minimal Newick trees: parse, write, restrict, bipartitions, FN/FP rates, attach leaves."""

import random


class Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self, name=None, length=None):
        self.name, self.length, self.children, self.parent = name, length, [], None

    def add(self, child):
        child.parent = self
        self.children.append(child)
        return child

    def leaves(self):
        out, stack = [], [self]
        while stack:
            n = stack.pop()
            if n.children:
                stack.extend(n.children)
            else:
                out.append(n)
        return out

    def postorder(self):
        out, stack = [], [(self, False)]
        while stack:
            n, seen = stack.pop()
            if seen or not n.children:
                out.append(n)
            else:
                stack.append((n, True))
                stack.extend((c, False) for c in n.children)
        return out


def parse(text):
    """Recursive-descent Newick parser (names, branch lengths; internal labels kept as names)."""
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]
    pos = [0]

    def label():
        j = pos[0]
        while j < len(s) and s[j] not in ",():;":
            j += 1
        tok = s[pos[0]:j].strip().strip("'")
        pos[0] = j
        return tok or None

    def length():
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            j = pos[0]
            while j < len(s) and s[j] not in ",();":
                j += 1
            v = s[pos[0]:j].strip()
            pos[0] = j
            return float(v) if v else None
        return None

    def node():
        n = Node()
        if s[pos[0]] == "(":
            pos[0] += 1
            while True:
                n.add(node())
                if s[pos[0]] == ",":
                    pos[0] += 1
                    continue
                assert s[pos[0]] == ")", s[pos[0]:pos[0] + 20]
                pos[0] += 1
                break
            n.name = label()  # internal label (or support value)
        else:
            n.name = label()
        n.length = length()
        return n

    return node()


def read(path):
    return parse(open(path).read())


def write(root, lengths=True):
    def rec(n):
        lab = "(" + ",".join(rec(c) for c in n.children) + ")" if n.children else (n.name or "")
        if lengths and n.length is not None and n.parent is not None:
            lab += ":%.6g" % n.length
        return lab
    return rec(root) + ";"


def suppress_unifurcations(root):
    for n in root.postorder():
        if n.children and len(n.children) == 1:
            c = n.children[0]
            if n.parent is None:
                if c.children:
                    root.children = c.children
                    for g in c.children:
                        g.parent = root
                continue
            p = n.parent
            idx = p.children.index(n)
            c.length = (c.length or 0) + (n.length or 0)
            c.parent = p
            p.children[idx] = c
    while root.children and len(root.children) == 1 and root.children[0].children:
        c = root.children[0]
        root.children = c.children
        for g in c.children:
            g.parent = root
    return root


def restrict(root, keep):
    """Copy of the tree induced on leaf names `keep`."""
    keep = set(keep)
    def rec(n):
        if not n.children:
            return Node(n.name, n.length) if n.name in keep else None
        kids = [k for k in (rec(c) for c in n.children) if k is not None]
        if not kids:
            return None
        m = Node(None, n.length)
        for k in kids:
            m.add(k)
        return m
    r = rec(root)
    return suppress_unifurcations(r)


def bipartitions(root, taxa=None):
    """Set of frozensets (the side not containing a fixed reference leaf) of nontrivial splits."""
    leaves = sorted(l.name for l in root.leaves())
    taxa = set(leaves) if taxa is None else set(taxa)
    ref = min(taxa)
    n = len(taxa)
    clade = {}
    out = set()
    for node in root.postorder():
        if not node.children:
            clade[node] = frozenset([node.name]) if node.name in taxa else frozenset()
        else:
            s = frozenset().union(*(clade[c] for c in node.children))
            clade[node] = s
            if node.parent is not None:
                side = s if ref not in s else frozenset(taxa - s)
                if 2 <= len(side) <= n - 2:
                    out.add(side)
    return out


def fn_fp(true_root, est_root, taxa=None):
    if taxa is None:
        taxa = {l.name for l in true_root.leaves()} & {l.name for l in est_root.leaves()}
    t = bipartitions(restrict(true_root, taxa), taxa)
    e = bipartitions(restrict(est_root, taxa), taxa)
    fn = len(t - e) / max(1, len(t))
    fp = len(e - t) / max(1, len(e))
    return fn, fp


def attach_leaves(root, names, pendant_lengths, rng):
    """Attach new leaves at uniformly random edges (by count), splitting each edge at a random point."""
    for name, plen in zip(names, pendant_lengths):
        edges = [n for n in root.postorder() if n.parent is not None]
        e = rng.choice(edges)
        p = e.parent
        L = e.length or 0.0
        x = rng.random() * L
        mid = Node(None, L - x)
        idx = p.children.index(e)
        p.children[idx] = mid
        mid.parent = p
        e.length = x
        mid.add(e)
        mid.add(Node(name, plen))
    return root


def root_to_tip(root):
    d = {root: 0.0}
    out = {}
    stack = [root]
    while stack:
        n = stack.pop()
        for c in n.children:
            d[c] = d[n] + (c.length or 0.0)
            stack.append(c)
        if not n.children:
            out[n.name] = d[n]
    return out
