"""Small tree helpers: bipartitions, RF (FN/FP), ASTRID internode-distance matrix."""
import itertools
import numpy as np
import treeswift


def read_newick_list(fn):
    txt = open(fn).read()
    return [t.strip() + ";" for t in txt.split(";") if t.strip()]


def _tree(nwk):
    t = treeswift.read_tree_newick(nwk)
    t.suppress_unifurcations()
    return t


def bipartitions(nwk):
    """Set of nontrivial bipartitions as frozensets of the side not containing the min label."""
    t = _tree(nwk)
    leaves = sorted(str(l.label) for l in t.traverse_leaves())
    ref = leaves[0]
    full = frozenset(leaves)
    out = set()
    below = {}
    for n in t.traverse_postorder():
        if n.is_leaf():
            below[n] = frozenset([str(n.label)])
        else:
            below[n] = frozenset().union(*(below[c] for c in n.children))
        s = below[n]
        if 1 < len(s) < len(full) - 1:
            side = s if ref not in s else full - s
            out.add(side)
    return out


def rf(true_bp, est_bp):
    """(false negatives, false positives, #true internal edges)."""
    return len(true_bp - est_bp), len(est_bp - true_bp), len(true_bp)


def astrid_matrix(genes):
    """ASTRID / NJst: average topological internode distance (missing pairs ignored)."""
    names = sorted({str(l.label) for g in genes for l in _tree(g).traverse_leaves()})
    idx = {n: i for i, n in enumerate(names)}
    S = np.zeros((len(names), len(names)))
    C = np.zeros_like(S)
    for g in genes:
        t = _tree(g)
        for e in t.traverse_preorder():
            e.edge_length = 1.0
        dm = t.distance_matrix(leaf_labels=True)
        for a, row in dm.items():
            for b, v in row.items():
                S[idx[str(a)], idx[str(b)]] += v
                C[idx[str(a)], idx[str(b)]] += 1
    M = np.where(C > 0, S / np.maximum(C, 1), 0)
    return names, M
