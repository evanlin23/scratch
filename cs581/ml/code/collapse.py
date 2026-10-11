"""Brown et al. (2017) baseline: the single tree with low-support edges collapsed.
FastTree trees carry SH-like local supports; collapse edges with support < S and score FN/FP/RF.

    python collapse.py DATASET REP [S ...]
"""
import os
import sys

import dendropy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runtrees as rt  # noqa: E402
import treeerr  # noqa: E402

ds, rep = sys.argv[1], sys.argv[2]
d = os.path.join(rt.MLDATA, ds, "R%s" % rep)
src = os.path.join(d, "trees", "true_align.fasttree.tre")
for S in [float(x) for x in sys.argv[3:]] or [0.5, 0.7, 0.8, 0.9, 0.95]:
    t = dendropy.Tree.get(path=src, schema="newick", preserve_underscores=True)
    for nd in list(t.postorder_internal_node_iter(exclude_seed_node=True)):
        if nd.label is not None and float(nd.label) < S:
            nd.edge.collapse()
    out = os.path.join(d, "trees", "true_align.fasttree_collapse%g.tre" % S)
    t.write(path=out, schema="newick", suppress_rooting=True, suppress_internal_node_labels=True)
    e = treeerr.error(os.path.join(d, "true_tree.tre"), out)
    print("%s R%s FastTree collapsed <%.2f: FN %.1f FP %.1f RF %.1f edges %d" % (
        ds, rep, S, 100 * e["fn_rate"], 100 * e["fp_rate"], 100 * e["rf_rate"], e["est_int"]))
