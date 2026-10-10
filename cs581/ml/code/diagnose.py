"""Where does the error on fragmentary data come from?

    python diagnose.py DATASET REP METHOD [METHOD ...]

For each method's tree (trees/true_align.<METHOD>.tre) reports
  FN overall, FN restricted to the backbone (full-length) taxa, and
  the mean fragment placement error: for every fragment f, prune all other fragments,
  and count the true backbone+f bipartitions that are missing beyond the backbone-only
  error ("delta error" per fragment, as in phylogenetic placement papers).
"""
import os
import statistics
import sys

import dendropy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runtrees as rt  # noqa: E402

MLDATA = rt.MLDATA


def splits(tree, index):
    """All edge splits of `tree` as int bitmasks over `index` (label -> bit position)."""
    out, below = [], {}
    for nd in tree.postorder_node_iter():
        if nd.is_leaf():
            below[nd] = 1 << index[nd.taxon.label]
        else:
            below[nd] = 0
            for c in nd.child_node_iter():
                below[nd] |= below[c]
        if nd is not tree.seed_node:
            out.append(below[nd])
    return out


def restrict(split_list, mask):
    """Nontrivial splits of the tree induced on the taxa in `mask`, canonicalised."""
    ref = mask & -mask
    n = bin(mask).count("1")
    out = set()
    for s in split_list:
        x = s & mask
        if x & ref:
            x = mask & ~x
        k = bin(x).count("1")
        if 2 <= k <= n - 2:
            out.add(x)
    return out


def main(ds, rep, methods):
    d = os.path.join(MLDATA, ds, "R%s" % rep)
    names, seqs = rt.read_fasta(os.path.join(d, "true_align.fasta"))
    lens = {n: len(s.replace("-", "")) for n, s in zip(names, seqs)}
    med = statistics.median(lens.values())
    bb = [n for n in names if lens[n] >= 0.5 * med]
    fr = [n for n in names if lens[n] < 0.5 * med]
    index = {n: i for i, n in enumerate(names)}
    bmask = sum(1 << index[n] for n in bb)
    true = splits(dendropy.Tree.get(path=os.path.join(d, "true_tree.tre"), schema="newick",
                                    preserve_underscores=True), index)
    tb = restrict(true, bmask)
    tfr = {f: restrict(true, bmask | 1 << index[f]) for f in fr}
    print("%s R%s: %d backbone, %d fragments" % (ds, rep, len(bb), len(fr)))
    for m in methods:
        p = os.path.join(d, "trees", "true_align.%s.tre" % m)
        if not os.path.exists(p):
            continue
        est = splits(dendropy.Tree.get(path=p, schema="newick", preserve_underscores=True), index)
        eb = restrict(est, bmask)
        bfn = len(tb - eb) / len(tb)
        deltas = []
        for f in fr:
            ef = restrict(est, bmask | 1 << index[f])
            deltas.append(len(tfr[f] - ef) - len(tb - eb))
        print("  %-28s backbone-FN %.1f%%  mean fragment delta-error %.2f edges  (frac fragments misplaced %.0f%%)" % (
            m, 100 * bfn, statistics.mean(deltas), 100 * sum(x > 0 for x in deltas) / len(deltas)))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
