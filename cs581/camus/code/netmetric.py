"""Rooted phylogenetic network utilities: extended-newick parsing, displayed trees,
softwired / hardwired clusters, and cluster FN/FP rates.

Clusters are leaf sets below a node. Trivial clusters (singletons and the full leaf
set) are excluded. Softwired clusters = union of the clusters of all displayed trees
(2^r switchings; fine for the r <= 6 networks used here).
"""
import itertools
import re
import sys


class Node:
    __slots__ = ("label", "children", "hyb", "gamma")

    def __init__(self, label=None):
        self.label = label
        self.children = []
        self.hyb = None  # e.g. "H46" for a reticulation node


def parse_enewick(s):
    """Return (root, hybrid_dict). Hybrid nodes with the same #H name are merged into one
    node object (the occurrence that carries children keeps them)."""
    s = s.strip().rstrip(";")
    i = 0
    hyb = {}
    last_gamma = [None]
    gamma = {}  # (id(parent), id(hybrid)) -> inheritance probability

    def read_label():
        nonlocal i
        j = i
        while i < len(s) and s[i] not in ",():;[":
            i += 1
        return s[j:i]

    def skip_meta():
        nonlocal i
        # branch length / support / gamma: everything after ':' up to , or )
        j = i
        while i < len(s) and s[i] == ":":
            i += 1
            while i < len(s) and s[i] not in ",():;[":
                i += 1
        if i < len(s) and s[i] == "[":
            i = s.index("]", i) + 1
        parts = s[j:i].split(":")
        return float(parts[3]) if len(parts) > 3 and parts[3] else None

    def node():
        nonlocal i
        n = Node()
        if s[i] == "(":
            i += 1
            while True:
                last_gamma[0] = None
                c = node()
                if last_gamma[0] is not None:
                    gamma[(id(n), id(c))] = last_gamma[0]
                n.children.append(c)
                if s[i] != ",":
                    break
                i += 1
            assert s[i] == ")", s[i - 20:i + 20]
            i += 1
        lab = read_label()
        g = skip_meta()
        m = re.match(r"^(.*?)#([A-Za-z]*\d+)$", lab)
        if m:
            name = m.group(2)
            last_gamma[0] = g
            if name in hyb:
                h = hyb[name]
                if n.children:
                    h.children = n.children
                return h
            n.hyb = name
            n.label = m.group(1) or None
            hyb[name] = n
            return n
        n.label = lab if lab != "" else None
        last_gamma[0] = None
        return n

    root = node()
    root.gamma = gamma
    return root, hyb


def _parents(root):
    par = {}
    seen = set()
    stack = [root]
    while stack:
        v = stack.pop()
        if id(v) in seen:
            continue
        seen.add(id(v))
        for c in v.children:
            par.setdefault(id(c), []).append(v)
            stack.append(c)
    return par


def leaves(root):
    out = set()
    seen = set()
    stack = [root]
    while stack:
        v = stack.pop()
        if id(v) in seen:
            continue
        seen.add(id(v))
        if not v.children:
            out.add(v.label)
        stack.extend(v.children)
    return frozenset(out)


def hardwired_clusters(root):
    memo = {}

    def rec(v):
        if id(v) in memo:
            return memo[id(v)]
        if not v.children:
            r = frozenset([v.label])
        else:
            r = frozenset().union(*[rec(c) for c in v.children])
        memo[id(v)] = r
        return r

    rec(root)
    allL = memo[id(root)]
    return {c for c in memo.values() if 1 < len(c) < len(allL)}


def displayed_cluster_sets(root):
    """List of cluster sets, one per displayed tree (one per choice of parent edge
    for every reticulation node)."""
    par = _parents(root)
    hnodes = {}
    stack, seen = [root], set()
    while stack:
        v = stack.pop()
        if id(v) in seen:
            continue
        seen.add(id(v))
        if len(par.get(id(v), [])) > 1:
            hnodes[id(v)] = par[id(v)]
        stack.extend(v.children)
    keys = list(hnodes)
    allL = leaves(root)
    out = []
    for choice in itertools.product(*[range(len(hnodes[k])) for k in keys]):
        keep = {k: id(hnodes[k][c]) for k, c in zip(keys, choice)}
        memo = {}

        def rec(v):
            if id(v) in memo:
                return memo[id(v)]
            if not v.children:
                r = frozenset([v.label])
            else:
                parts = []
                for c in v.children:
                    if id(c) in keep and keep[id(c)] != id(v):
                        continue
                    parts.append(rec(c))
                r = frozenset().union(*parts) if parts else frozenset()
            memo[id(v)] = r
            return r

        rec(root)
        out.append({c for c in memo.values() if 1 < len(c) < len(allL)})
    return out


def softwired_clusters(root):
    s = set()
    for cs in displayed_cluster_sets(root):
        s |= cs
    return s


def num_reticulations(root):
    par = _parents(root)
    return sum(len(p) - 1 for p in par.values() if len(p) > 1)


def major_tree_newick(root):
    """Displayed tree keeping, at each reticulation, the parent edge with the larger gamma."""
    par = _parents(root)
    g = getattr(root, "gamma", {}) or {}
    best = {}
    for cid, ps in par.items():
        if len(ps) > 1:
            best[cid] = max(ps, key=lambda p: g.get((id(p), cid), 0.5))
    choice_keep = {cid: id(p) for cid, p in best.items()}
    return displayed_tree_newicks(root, choice_keep)[0]


def displayed_tree_newicks(root, fixed=None):
    """Newick strings (topology only) of all displayed trees, degree-2 nodes suppressed."""
    par = _parents(root)
    hnodes = {}
    stack, seen = [root], set()
    while stack:
        v = stack.pop()
        if id(v) in seen:
            continue
        seen.add(id(v))
        if len(par.get(id(v), [])) > 1:
            hnodes[id(v)] = par[id(v)]
        stack.extend(v.children)
    keys = list(hnodes)
    res = []
    choices = [None] if fixed is not None else itertools.product(*[range(len(hnodes[k])) for k in keys])
    for choice in choices:
        keep = fixed if fixed is not None else {k: id(hnodes[k][c]) for k, c in zip(keys, choice)}

        def rec(v):
            if not v.children:
                return v.label
            parts = [rec(c) for c in v.children if not (id(c) in keep and keep[id(c)] != id(v))]
            parts = [p for p in parts if p is not None]
            if not parts:
                return None
            if len(parts) == 1:
                return parts[0]
            return "(" + ",".join(parts) + ")"

        res.append(rec(root) + ";")
    return res


def topology_enewick(s):
    """Strip branch lengths / gammas / supports (for PhyloNet)."""
    s = s.strip()
    s = re.sub(r"\[[^\]]*\]", "", s)
    s = re.sub(r":[^,();]*", "", s)
    s = re.sub(r"\)[0-9.eE+-]+", ")", s)  # internal support labels
    return s


def cluster_error(true_cl, est_cl):
    fn = len(true_cl - est_cl) / len(true_cl) if true_cl else 0.0
    fp = len(est_cl - true_cl) / len(est_cl) if est_cl else 0.0
    return fn, fp


def compare(true_nwk, est_nwk, mode="soft"):
    t, _ = parse_enewick(true_nwk)
    e, _ = parse_enewick(est_nwk)
    assert leaves(t) == leaves(e), (sorted(leaves(t) ^ leaves(e)))
    f = softwired_clusters if mode == "soft" else hardwired_clusters
    return cluster_error(f(t), f(e))


if __name__ == "__main__":
    a, b = open(sys.argv[1]).readline(), open(sys.argv[2]).readline()
    print("soft FN FP", compare(a, b, "soft"), "hard FN FP", compare(a, b, "hard"))
