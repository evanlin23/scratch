"""Minimal multispecies-coalescent gene-tree simulator (one lineage per species).

Species tree: rooted binary, branch lengths in coalescent units. Gene trees are returned
as newick topologies (no branch lengths needed for summary methods).
"""
import numpy as np


class Node:
    __slots__ = ("children", "length", "label")

    def __init__(self, children=(), length=0.0, label=None):
        self.children = list(children)
        self.length = length
        self.label = label


def caterpillar(n, f, pendant=None):
    """((((t0,t1),t2),t3)...) with every internal branch = f. Leaf ages irrelevant (topology only)."""
    node = Node([Node(label="t0"), Node(label="t1")], length=f)
    for i in range(2, n):
        node = Node([node, Node(label=f"t{i}")], length=f)
    node.length = np.inf
    return node


def balanced(n, f):
    """Perfectly balanced tree on n = 2^h leaves, every internal branch = f."""
    leaves = [Node(label=f"t{i}") for i in range(n)]
    level = leaves
    while len(level) > 1:
        level = [Node([level[i], level[i + 1]], length=f) for i in range(0, len(level), 2)]
    level[0].length = np.inf
    return level[0]


def to_newick(node):
    if not node.children:
        return node.label
    return "(" + ",".join(to_newick(c) for c in node.children) + ")"


def species_newick(root):
    return to_newick(root) + ";"


def sim_gene(root, rng):
    """Return newick of one gene tree under the MSC."""
    def rec(node):
        if not node.children:
            lin = [node.label]
        else:
            lin = []
            for c in node.children:
                lin.extend(rec(c))
        # coalesce within this branch of length node.length
        t = 0.0
        while len(lin) > 1:
            k = len(lin)
            t += rng.exponential(2.0 / (k * (k - 1)))
            if t > node.length:
                break
            i, j = rng.choice(k, 2, replace=False)
            a, b = lin[i], lin[j]
            lin = [x for idx, x in enumerate(lin) if idx not in (i, j)] + [f"({a},{b})"]
        return lin
    out = rec(root)
    assert len(out) == 1
    return out[0] + ";"


def sim_genes(root, k, seed):
    rng = np.random.default_rng(seed)
    return [sim_gene(root, rng) for _ in range(k)]


if __name__ == "__main__":
    # sanity check: quartet ((a,b),c),d with internal branch f -> P(gene topology = species) = 1 - 2/3 e^-f
    import sys
    for f in [0.05, 0.2, 1.0]:
        sp = caterpillar(4, f)
        # 4-taxon caterpillar has 2 internal branches of length f; unrooted quartet topology ab|cd
        # is affected only by the branch above (t0,t1): P(ab|cd) = 1 - 2/3 exp(-f)
        g = sim_genes(sp, 20000, 1)
        import treeswift as ts
        hit = 0
        for nw in g:
            t = ts.read_tree_newick(nw)
            ok = False
            for nd in t.traverse_internal():
                s = {l.label for l in nd.traverse_leaves()}
                if s in ({"t0", "t1"}, {"t2", "t3"}) or s in ({"t0", "t1", "t2", "t3"} - {"t0", "t1"},):
                    ok = True
            hit += ok
        print(f, hit / len(g), 1 - 2 / 3 * np.exp(-f))
