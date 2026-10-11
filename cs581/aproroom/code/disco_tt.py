"""DISCO decomposition with GIVEN root and D/S tags (internal labels from truetag.py), using DISCO v1.4.1's
own decompose(). Usage: python disco_tt.py TAGGED.trees OUT.trees   (output leaves = species, >= 4 leaves)"""
import sys
import treeswift
sys.path.insert(0, "/opt/src/DISCO")
from disco import decompose, unroot  # noqa: E402

sp = lambda x: x.split("_")[0]
with open(sys.argv[1]) as fi, open(sys.argv[2], "w") as fo:
    for line in fi:
        if ";" not in line:
            continue
        t = treeswift.read_tree_newick(line.strip())
        t.suppress_unifurcations()
        for n in t.traverse_postorder():
            if n.is_leaf():
                n.s = {sp(n.label)}
            else:
                n.s = set().union(*(c.s for c in n.children))
                n.tag = n.label
                n.label = None
        for o in decompose(t):
            if o.num_nodes(internal=False) < 4:
                continue
            unroot(o)
            for l in o.traverse_leaves():
                l.label = sp(l.label)
            o.suppress_unifurcations()
            fo.write(o.newick() + "\n")
