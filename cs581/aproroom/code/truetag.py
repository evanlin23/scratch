"""Label SimPhy true gene trees with their TRUE duplication/speciation tags.

Usage: python truetag.py S_TREE L_TREES G_TREES OUT [--genes IDS.txt]

A locus-tree node is a speciation iff its time before present equals that of a species-tree node (SimPhy
places speciations exactly at species-tree node times; duplications fall at continuous random times).
A gene-tree node g is tagged D iff the locus-tree LCA of the loci of g's leaves is a duplication node
(orthology defined by the locus tree; exact without ILS between gene and locus tree). The true root is kept.
Output: one rooted gene tree per line, internal nodes labelled D or S, leaf labels unchanged
(species_locus_individual). With --genes, only lines whose 1-based index is listed are written.
Prints, on stderr, a check: fraction of locus-tree S nodes whose children have disjoint species sets.
"""
import sys
import treeswift


def heights(t):
    """time before present for every node (present = deepest extant leaf)."""
    d = {}
    for n in t.traverse_preorder():
        d[n] = (0.0 if n.is_root() else d[n.parent] + (n.edge_length or 0.0))
    top = max(d[n] for n in t.traverse_leaves() if not str(n.label).startswith("Lost"))
    return {n: top - x for n, x in d.items()}


def main():
    a = sys.argv[1:]
    sfile, lfile, gfile, out = a[:4]
    keep = None
    if "--genes" in a:
        keep = set(int(x) for x in open(a[a.index("--genes") + 1]).read().split())
    st = treeswift.read_tree_newick(open(sfile).read().strip())
    sh = sorted(h for n, h in heights(st).items() if not n.is_leaf())
    import bisect

    def is_spec_time(h):
        i = bisect.bisect_left(sh, h)
        for j in (i - 1, i):
            if 0 <= j < len(sh) and abs(sh[j] - h) <= 1e-6 * max(1.0, sh[-1]):
                return True
        return False
    ok = bad = 0
    ndup = nnode = 0
    with open(lfile) as lf, open(gfile) as gf, open(out, "w") as fo:
        for i, (ll, gl) in enumerate(zip(lf, gf), 1):
            if keep is not None and i not in keep:
                continue
            lt = treeswift.read_tree_newick(ll.strip())
            lh = heights(lt)
            dupnode = {}
            for n in lt.traverse_postorder():
                if n.is_leaf():
                    n.sp = set() if str(n.label).startswith("Lost") else {str(n.label).split("_")[0]}
                    continue
                n.sp = set().union(*(c.sp for c in n.children))
                dupnode[n] = not is_spec_time(lh[n])
                if not dupnode[n] and len(n.children) == 2:
                    a1, a2 = n.children
                    ok += not (a1.sp & a2.sp)
                    bad += bool(a1.sp & a2.sp)
            leaf = {str(n.label): n for n in lt.traverse_leaves()}
            # locus-tree LCA via depth/parent walk
            depth = {}
            for n in lt.traverse_preorder():
                depth[n] = 0 if n.is_root() else depth[n.parent] + 1

            def lca(x, y):
                if x is None:
                    return y
                while depth[x] > depth[y]:
                    x = x.parent
                while depth[y] > depth[x]:
                    y = y.parent
                while x is not y:
                    x, y = x.parent, y.parent
                return x
            gt = treeswift.read_tree_newick(gl.strip())
            gt.suppress_unifurcations()
            for n in gt.traverse_postorder():
                if n.is_leaf():
                    lab = str(n.label)
                    n.L = leaf["_".join(lab.split("_")[:2])]
                    continue
                m = None
                for c in n.children:
                    m = lca(m, c.L)
                n.L = m
                n.label = "D" if dupnode.get(m, False) else "S"
                ndup += n.label == "D"
                nnode += 1
            fo.write(gt.newick().replace("[&R] ", "") + "\n")
    sys.stderr.write("locus S-node disjointness ok=%d bad=%d; gene nodes %d, D=%d\n" % (ok, bad, nnode, ndup))


if __name__ == "__main__":
    main()
