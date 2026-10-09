"""Validate forest_fast against lanl/distphylo's own functions and its shipped output.

Uses distphylo's example alignment (aln_128_500_1.fa, n=128, k=500) and true tree, runs
(a) distphylo's mini_contractor/extender/get_unique_trees (imported unmodified), and
(b) forest_fast.run_forest, for a set of (m, M, tau) grid points, and compares the number of
non-trivial splits and false splits with each other and with distphylo's shipped TSV
(grid_summary_ntips128_1_k500_sorted.tsv).

usage: python validate_forest.py DISTPHYLO_DIR OUT_TSV
"""
import sys, os, time, csv
import numpy as np
import networkx as nx
import dendropy

DP = sys.argv[1]
sys.path.insert(0, DP)
import prune_deep_helper_functions as H  # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import forest_fast as F  # noqa: E402
from common import read_fasta, jc_distance, true_splits  # noqa: E402


def distphylo_forest(D, m, M, tau):
    """Mirror of get_forest_split_check's forest construction (without R / ete3 scoring)."""
    g = H.construct_clustering_graph(D, m)
    comps = list(nx.connected_components(g))
    out = []
    for comp in comps:
        sub = g.subgraph(comp)
        nodes = list(sub.nodes())
        pairs = [(u, v) for u in nodes for v in nodes if u < v]
        final = []
        for lv in pairs:
            mb, _ = H.mini_contractor(sub, D, lv, M, tau)
            final.append(H.extender(sub, mb, lv, D))
        flat = [tuple(frozenset(p) for p in bp) for sl in final for bp in sl]
        uniq = list(set(flat))
        uniq = [tuple(set(p) for p in bp) for bp in uniq]
        sfz = set(H.sort_and_freeze(p) for p in uniq)
        tup = [tuple(map(set, p)) for p in sfz]
        if len(tup) == 0:
            out.append((sorted(comp), None))
            continue
        trees = H.get_unique_trees(tup, len(nodes), 10)
        if len(trees) != 1:
            return None
        out.append((sorted(comp), list(trees)[0]))
    return out


def nontrivial_from_newick(nwk, comp):
    """non-trivial split bitmasks (global idx, normalized away from min leaf) of a newick tree."""
    tns = dendropy.TaxonNamespace([str(i) for i in comp])
    t = dendropy.Tree.get(data=nwk, schema="newick", taxon_namespace=tns, rooting="force-unrooted")
    t.encode_bipartitions()
    full = F.mask_of(comp)
    res = set()
    for e in t.edges():
        if e.head_node.is_leaf() or e.tail_node is None:
            continue
        loc = e.bipartition.leafset_bitmask
        gm = 0
        for tx in t.taxon_namespace:
            if loc & t.taxon_namespace.taxon_bitmask(tx):
                gm |= 1 << int(tx.label)
        k = F.popcount(gm)
        if 2 <= k <= len(comp) - 2:
            if (gm >> min(comp)) & 1:
                gm = full & ~gm
            res.add(gm)
    return res


def false_count(splits, comp, T):
    """splits (global bitmasks within comp) not in the true tree restricted to comp."""
    full = F.mask_of(comp)
    Tr = set()
    for s in T:
        a = s & full
        k = F.popcount(a)
        if 2 <= k <= len(comp) - 2:
            if (a >> int(min(comp))) & 1:
                a = full & ~a
            Tr.add(a)
    return sum(1 for s in splits if s not in Tr)


def main():
    out_tsv = sys.argv[2]
    names, X = read_fasta(os.path.join(DP, "aln_128_500_1.fa"))
    order = np.argsort([int(x) for x in names])  # labels are 0..127
    X = X[order]
    names = [names[i] for i in order]
    D = jc_distance(X)
    tre = open(os.path.join(DP, "true_tree_128.tre")).read().strip()
    tns = dendropy.TaxonNamespace(names)
    tt = dendropy.Tree.get(data=tre, schema="newick", taxon_namespace=tns, rooting="force-unrooted")
    T = true_splits(tt, names)
    pub = {}
    with open(os.path.join(DP, "grid_summary_ntips128_1_k500_sorted.tsv")) as fh:
        lines = fh.read().split("\n\n", 1)[1].strip().splitlines()
        rd = csv.DictReader(lines, delimiter="\t")
        for r in rd:
            key = (round(float(r["m"]), 3), round(float(r["M"]), 3), round(float(r["tau"]), 3))
            pub[key] = (r["sum_split_count"], r["forest_sum_falsesplit_count"], r["num_components"])
    allgrid = [(0.8, 2.3, 0.1), (0.6, 2.0, 0.05), (0.4, 1.0, 0.045), (0.4, 1.0, 0.03), (0.8, 2.0, 0.08)]
    grid = [allgrid[int(i)] for i in sys.argv[3].split(",")] if len(sys.argv) > 3 else allgrid
    rows = []
    for (m, M, tau) in grid:
        t0 = time.time()
        dp = distphylo_forest(D, m, M, tau)
        t_dp = time.time() - t0
        if dp is None:
            dp_n, dp_f, dp_c = "invalid", "", ""
        else:
            sp = [nontrivial_from_newick(nw, c) if nw else set() for c, nw in dp]
            dp_n = sum(len(s) for s in sp)
            dp_f = sum(false_count(s, c, T) for s, (c, _) in zip(sp, dp))
            dp_c = len(dp)
        t0 = time.time()
        r = F.run_forest(D, m, M, [tau])[tau]
        t_ff = time.time() - t0
        if not r["valid"]:
            ff_n, ff_f = "invalid", ""
        else:
            ff_n = sum(len(s) for s in r["splits"])
            ff_f = sum(false_count(s, c, T) for s, c in zip(r["splits"], r["comps"]))
        key = (round(m, 3), round(M, 3), round(tau, 3))
        p = pub.get(key, ("NA", "NA", "NA"))
        rows.append(dict(m=m, M=M, tau=tau, pub_splits=p[0], pub_false=p[1], pub_ncomp=p[2],
                         distphylo_splits=dp_n, distphylo_false=dp_f, distphylo_ncomp=dp_c,
                         fast_splits=ff_n, fast_false=ff_f, fast_ncomp=len(r["comps"]),
                         t_distphylo=round(t_dp, 2), t_fast=round(t_ff, 2)))
        print(rows[-1], flush=True)
    with open(out_tsv, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
