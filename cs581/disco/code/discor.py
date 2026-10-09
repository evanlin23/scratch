"""DISCO-R: species-tree-guided rooting/tagging of gene-family trees before DISCO decomposition.

For each unrooted gene-family tree, every rooting is scored by duplication-loss parsimony under
LCA mapping against a reference species tree S0 (NOTUNG-style, O(n) rerooting DP). The best rooting
(ties broken by DISCO's own MinDL score) is then tagged either with DISCO's species-overlap rule
("overlap") or with the LCA-reconciliation rule ("lca"), and decomposed with DISCO's decompose().

Usage (library): see run_pipeline.py. CLI:
  python discor.py -i genes.tre -s S0.tre -o out.tre [--tag overlap|lca] [--no-decomp] [--dup 1.5 --loss 1]
"""
import argparse
import sys
import warnings

import treeswift

sys.path.insert(0, "/opt/tools/DISCO")
import disco  # noqa: E402  (DISCO v1.4.1)

sys.setrecursionlimit(100000)


def species_of(label):
    return label.split("_")[0]


class SpeciesTree:
    """Rooted species tree with depth + LCA via parent pointers (small trees)."""

    def __init__(self, newick, contract_below=None):
        t = treeswift.read_tree_newick(newick)
        if contract_below is not None:
            for n in list(t.traverse_postorder(leaves=False)):
                if n.is_root():
                    continue
                try:
                    sup = float(n.label)
                except (TypeError, ValueError):
                    continue
                if sup < contract_below:
                    n.contract()
        self.tree = t
        self.parent, self.depth, self.leaf = {}, {}, {}
        for i, n in enumerate(t.traverse_preorder()):
            n.idx = i
            self.parent[i] = None if n.is_root() else n.get_parent().idx
            self.depth[i] = 0 if n.is_root() else self.depth[self.parent[i]] + 1
            if n.is_leaf():
                self.leaf[n.label] = i
        self.root = t.root.idx
        self._cache = {}

    def lca(self, a, b):
        if a == b:
            return a
        key = (a, b) if a < b else (b, a)
        r = self._cache.get(key)
        if r is not None:
            return r
        x, y = a, b
        while self.depth[x] > self.depth[y]:
            x = self.parent[x]
        while self.depth[y] > self.depth[x]:
            y = self.parent[y]
        while x != y:
            x, y = self.parent[x], self.parent[y]
        self._cache[key] = x
        return x


def local_cost(S, ma, mb, wd, wl):
    """DL cost of joining child mappings ma, mb at a new node; returns (mapping, cost, is_dup)."""
    m = S.lca(ma, mb)
    dup = m == ma or m == mb
    d = S.depth[m]
    off = 0 if dup else 1
    losses = (S.depth[ma] - d - off) + (S.depth[mb] - d - off)
    return m, wd * dup + wl * losses, dup


def best_dl_root(tree, S, wd=1.5, wl=1.0):
    """Score all rootings of `tree` (treeswift, any rooting) by DL parsimony against S.
    Returns (best node to reroot above, best cost, list of tied nodes)."""
    if tree.root.num_children() != 2:
        disco.reroot_on_edge(tree, tree.root.child_nodes()[0])
    tree.resolve_polytomies()
    for n in tree.traverse_postorder():
        if n.is_leaf():
            n.dm = S.leaf[species_of(n.label)]
            n.dc = 0.0
        else:
            a, b = n.child_nodes()
            m, c, _ = local_cost(S, a.dm, b.dm, wd, wl)
            n.dm, n.dc = m, a.dc + b.dc + c
    root = tree.root
    left, right = root.child_nodes()
    left.um, left.uc = right.dm, right.dc
    right.um, right.uc = left.dm, left.dc
    best, best_cost, ties = None, float("inf"), []
    for n in tree.traverse_preorder():
        if n.is_root():
            continue
        p = n.get_parent()
        if not p.is_root():
            sib = [c for c in p.child_nodes() if c is not n][0]
            m, c, _ = local_cost(S, p.um, sib.dm, wd, wl)
            n.um, n.uc = m, p.uc + sib.dc + c
        # rooting on edge (n, parent): children = subtree(n), rest
        if p.is_root() and n is right:
            continue  # same edge as root-left
        _, c, _ = local_cost(S, n.dm, n.um, wd, wl)
        tot = n.dc + n.uc + c
        if tot < best_cost - 1e-9:
            best, best_cost, ties = n, tot, [n]
        elif abs(tot - best_cost) <= 1e-9:
            ties.append(n)
    return best, best_cost, ties


def disco_scores(tree):
    """DISCO MinDL score of each rooting, keyed by id(node) of the node above which we root."""
    _, _, _ = disco.get_min_root(tree, species_of)
    sc = {}

    def score(total, s1, s2):
        if s1 & s2:
            if total == s1 or total == s2:
                return 1 if s1 == s2 else 2
            return 3
        return 0

    for n in tree.traverse_preorder():
        if n.is_root() or not hasattr(n, "up"):
            continue
        sc[id(n)] = n.u_score + n.d_score + score(n.up | n.down, n.up, n.down)
    return sc


def lca_tag(tree, S):
    """Tag rooted binary tree by LCA reconciliation (dup iff M(v) == M(child)); sets node.s like DISCO."""
    tree.suppress_unifurcations()
    tree.resolve_polytomies()
    for n in tree.traverse_postorder():
        if n.is_leaf():
            n.s = {species_of(n.label)}
            n.m = S.leaf[species_of(n.label)]
        else:
            a, b = n.child_nodes()
            n.s = a.s | b.s
            n.m, _, dup = local_cost(S, a.m, b.m, 1, 1)
            n.tag = "D" if dup else "S"


def root_and_tag(tree, mode, S=None, tagmode="overlap", wd=1.5, wl=1.0):
    """mode: 'disco' (MinDL) or 'dl' (DL parsimony vs S, ties -> DISCO score)."""
    if tree.root.num_children() == 0:
        return tree
    if mode == "disco":
        r, _, _ = disco.get_min_root(tree, species_of)
    else:
        r, _, ties = best_dl_root(tree, S, wd, wl)
        if len(ties) > 1:
            sc = disco_scores(tree)
            r = min(ties, key=lambda n: (sc.get(id(n), 99999), 0))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        disco.reroot_on_edge(tree, r)
    if tagmode == "lca":
        lca_tag(tree, S)
    else:
        disco.tag(tree, species_of)
    return tree


def decompose_trees(trees, minimum=4):
    out = []
    for t in trees:
        if t.root.num_children() == 0:
            continue
        for d in disco.decompose(t):
            if d.num_nodes(internal=False) >= minimum:
                disco.unroot(d)
                disco.relabel(d, species_of)
                d.suppress_unifurcations()
                out.append(d.newick())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", required=True)
    ap.add_argument("-s", help="reference species tree (omit -> plain DISCO rooting)")
    ap.add_argument("-o", required=True)
    ap.add_argument("--tag", default="overlap", choices=["overlap", "lca"])
    ap.add_argument("--dup", type=float, default=1.5)
    ap.add_argument("--loss", type=float, default=1.0)
    ap.add_argument("--contract", type=float, default=None)
    ap.add_argument("--no-decomp", action="store_true")
    a = ap.parse_args()
    S = SpeciesTree(open(a.s).readline(), a.contract) if a.s else None
    trees = []
    for line in open(a.i):
        if line.strip():
            t = treeswift.read_tree_newick(line)
            root_and_tag(t, "dl" if S else "disco", S, a.tag, a.dup, a.loss)
            trees.append(t)
    with open(a.o, "w") as fo:
        if a.no_decomp:
            for t in trees:
                fo.write(t.newick() + "\n")
        else:
            for nw in decompose_trees(trees):
                fo.write(nw + "\n")


if __name__ == "__main__":
    main()
