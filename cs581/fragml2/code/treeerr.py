"""Missing-branch (FN), false-positive (FP) and RF error of estimated trees vs a reference tree.

    python treeerr.py TRUE_TREE EST_TREE [EST_TREE ...]

Both trees are restricted to their common leaf set and treated as unrooted. Rates are
normalised by the number of internal edges of the (binary) reference: n - 3.
"""
import sys

import dendropy


def load(path, tns):
    t = dendropy.Tree.get(path=path, schema="newick", taxon_namespace=tns,
                          preserve_underscores=True, rooting="force-unrooted")
    return t


def error(true_path, est_path):
    tns = dendropy.TaxonNamespace()
    t = load(true_path, tns)
    e = load(est_path, tns)
    common = {x.taxon.label for x in t.leaf_node_iter()} & {x.taxon.label for x in e.leaf_node_iter()}
    for tree in (t, e):
        tree.retain_taxa_with_labels(common)
        tree.encode_bipartitions()
    fp, fn = dendropy.calculate.treecompare.false_positives_and_negatives(t, e, is_bipartitions_updated=True)
    ti = len([b for b in t.bipartition_encoding if not b.is_trivial()])
    ei = len([b for b in e.bipartition_encoding if not b.is_trivial()])
    n = len(common)
    return {"n": n, "fn": fn, "fp": fp, "true_int": ti, "est_int": ei,
            "fn_rate": fn / ti if ti else 0.0, "fp_rate": fp / ei if ei else 0.0,
            "rf_rate": (fn + fp) / (2.0 * (n - 3))}


if __name__ == "__main__":
    for p in sys.argv[2:]:
        r = error(sys.argv[1], p)
        print("%s\tFN=%.4f\tFP=%.4f\tRF=%.4f\t(n=%d, true internal edges=%d, est=%d)" %
              (p, r["fn_rate"], r["fp_rate"], r["rf_rate"], r["n"], r["true_int"], r["est_int"]))
