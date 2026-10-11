"""Quartet topology counts from gene trees (numpy, four-point condition on topological
distances) and the CAMUS 'max' objective (number of filtered gene-tree quartets displayed
by a network) for comparing networks built on different base trees."""
import itertools
import numpy as np
import treeswift as ts


def tree_dist(newick, taxa):
    """Topological (edge-count) distance matrix over `taxa` (unrooted)."""
    t = ts.read_tree_newick(newick)
    for n in t.traverse_preorder():
        n.edge_length = 1.0
    d = t.distance_matrix(leaf_labels=True)
    idx = {x: i for i, x in enumerate(taxa)}
    D = np.full((len(taxa), len(taxa)), np.nan)
    for a, row in d.items():
        for b, v in row.items():
            D[idx[a], idx[b]] = v
    np.fill_diagonal(D, 0)
    return D


class QuartetTable:
    def __init__(self, taxa):
        self.taxa = list(taxa)
        self.sets = np.array(list(itertools.combinations(range(len(taxa)), 4)), dtype=np.int32)
        a, b, c, d = self.sets.T
        self.abcd = (a, b, c, d)

    def topo(self, D):
        """0: ab|cd, 1: ac|bd, 2: ad|bc, -1 unresolved/missing."""
        a, b, c, d = self.abcd
        s = np.stack([D[a, b] + D[c, d], D[a, c] + D[b, d], D[a, d] + D[b, c]], axis=1)
        top = np.argmin(s, axis=1)
        srt = np.sort(s, axis=1)
        bad = ~(srt[:, 0] < srt[:, 1]) | np.isnan(srt).any(axis=1)
        top[bad] = -1
        return top

    def counts(self, newicks):
        C = np.zeros((len(self.sets), 3), dtype=np.int64)
        r = np.arange(len(self.sets))
        for nw in newicks:
            tp = self.topo(tree_dist(nw, self.taxa))
            ok = tp >= 0
            np.add.at(C, (r[ok], tp[ok]), 1)
        return C


def filter_mask(C, t=0.5, z=0.0, zmode="and"):
    """Boolean mask [n,3] of quartet topologies kept by CAMUS' filter (-q 2):
    dominant always; second if (c2-c3) > t*(c2+c3) [and/or/only a z-test]."""
    order = np.argsort(-C, axis=1, kind="stable")
    srt = np.take_along_axis(C, order, axis=1)
    c1, c2, c3 = srt[:, 0], srt[:, 1], srt[:, 2]
    s = c2 + c3
    ratio = np.floor(t * s) < (c2 - c3)
    if z > 0:
        sig = (s > 0) & ((c2 - c3 - 1) / np.sqrt(np.maximum(s, 1)) > z)
        keep2 = {"only": sig, "or": ratio | sig}.get(zmode, ratio & sig)
    else:
        keep2 = ratio
    M = np.zeros_like(C, dtype=bool)
    r = np.arange(len(C))
    M[r, order[:, 0]] = c1 > 0
    M[r[keep2], order[keep2, 1]] = True
    return M


def network_score(qt, C, M, displayed_newicks):
    """Sum of kept gene-tree quartet counts whose topology is displayed by any displayed tree."""
    disp = np.zeros_like(M)
    r = np.arange(len(C))
    for nw in displayed_newicks:
        tp = qt.topo(tree_dist(nw, qt.taxa))
        ok = tp >= 0
        disp[r[ok], tp[ok]] = True
    return int((C * (disp & M)).sum())
