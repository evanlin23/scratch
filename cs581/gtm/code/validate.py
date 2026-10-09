"""Rescore the published trees of Park, Zaharias & Warnow 2021 (Algorithms 14:148,
data doi:10.13012/B2IDB-7008049_V1) against the true trees, FN rate.
Usage: python3 validate.py DATA_ROOT OUT_TSV
"""
import os
import sys
from phylo import read_tree, fn_fp

ROOT = sys.argv[1]
R = os.path.join(ROOT, "d/RNASim1000")
C = os.path.join(ROOT, "cox/Cox1-Het")
CT = os.path.join(ROOT, "cox/Cox1-HET")
M = os.path.join(ROOT, "m1/1000M1_HF")
MT = os.path.join(ROOT, "m1_true")

# published FN (%) from Tables 3-5 (+ Table A2 for TreeMerge-PAUP)
PUB = {
    "RNASim1000": {"FastTree": 14.9, "IQ-TREE": 15.1, "RAxML-NG": 15.1,
                   "CINC/FT": 15.1, "GTM/FT": 14.7, "TreeMerge/FT": 14.7,
                   "CINC/IQ": 14.5, "GTM/IQ": 14.4, "TreeMerge/IQ": 14.4,
                   "TreeMerge-PAUP/FT": 14.7},
    "Cox1-HET": {"FastTree": 23.9, "IQ-TREE": 19.6, "RAxML-NG": 18.2,
                 "CINC/FT": 18.9, "GTM/FT": 18.9, "TreeMerge/FT": 18.9,
                 "CINC/IQ": 19.7, "GTM/IQ": 18.7, "TreeMerge/IQ": 18.7,
                 "TreeMerge-PAUP/FT": 18.9},
    "1000M1-HF": {"FastTree": 50.9, "IQ-TREE": 30.2, "RAxML-NG": 24.9,
                  "CINC/FT": 42.4, "GTM/FT": 42.4, "TreeMerge/FT": 42.5,
                  "CINC/IQ": 28.6, "GTM/IQ": 28.4, "TreeMerge/IQ": 28.5,
                  "TreeMerge-PAUP/FT": 42.5},
}


def rnasim(r):
    cc = f"{R}/CreateConstraintTrees/500"
    return f"{R}/{r}/model/true.tt", {
        "FastTree": f"{cc}/FastTree/{r}/output/fasttree.out",
        "IQ-TREE": f"{cc}/IQTree2/{r}/output/iqtree-full.treefile",
        "RAxML-NG": f"{R}/RAxML-ng/{r}/raxmlng-result.raxml.lastTree.TMP",
        "CINC/FT": f"{R}/Constrained-INC/fasttree_fasttree/500/{r}/node.",
        "GTM/FT": f"{R}/GTM/fasttree_fasttree/500/{r}/branch_length.",
        "TreeMerge/FT": f"{R}/TreeMerge_RAxML-ng/fasttree_fasttree/500/{r}/output/treemerge_node.tree",
        "CINC/IQ": f"{R}/Constrained-INC/IQTree2/500/{r}/node.",
        "GTM/IQ": f"{R}/GTM/IQTree2/500/{r}/branch_length.",
        "TreeMerge/IQ": f"{R}/TreeMerge_RAxML-ng/IQTree2/500/{r}/output/treemerge_node.tree",
        "TreeMerge-PAUP/FT": f"{R}/TreeMerge_PAUP/{r}/output/treemerge.tree",
    }


def cox(r):
    ct = f"{C}/ConstraintTrees/{r}"
    return f"{CT}/{r}/true-tree.tre", {
        "FastTree": f"{ct}/fasttree-centro-500/fasttree.tre",
        "IQ-TREE": f"{ct}/iqtree+GTR+G-centro-500/iqtree-result.treefile",
        "RAxML-NG": f"{C}/RAxML-ng/{r}/raxml+GTR+G.raxml.bestTree",
        "CINC/FT": f"{C}/Constrained-INC/fasttree_fasttree/500/{r}/node_dist.",
        "GTM/FT": f"{C}/GTM/500/fasttree_fasttree/{r}/branch_length.",
        "TreeMerge/FT": f"{C}/TreeMerge_RAxML-ng/fasttree_fasttree/{r}/output/treemerge_node.tree",
        "CINC/IQ": f"{C}/Constrained-INC/IQTree2/500/{r}/node_dist.",
        "GTM/IQ": f"{C}/GTM/500/iqtree/{r}/branch_length.",
        "TreeMerge/IQ": f"{C}/TreeMerge_RAxML-ng/IQTree2/{r}/output/treemerge_node.tree",
        "TreeMerge-PAUP/FT": f"{C}/TreeMerge_PAUP/{r}/output/treemerge.tree",
    }


def m1(r):
    cc = f"{M}/CreateConstraintTrees"
    return f"{MT}/{r}/rose.mt", {
        "FastTree": f"{cc}/FastTree/500/{r}/output/fasttree.out",
        "IQ-TREE": f"{M}/IQTree/{r}/iqtree-result.treefile",
        "RAxML-NG": f"{M}/RAxML-ng/{r}/raxmlng-result.raxml.lastTree.TMP",
        "CINC/FT": f"{M}/Constrained-INC/fasttree_fasttree/500/{r}/node.",
        "GTM/FT": f"{M}/GTM/500/fasttree_fasttree/{r}/branch_length.",
        "TreeMerge/FT": f"{M}/TreeMerge_RAxML-ng/fasttree_fasttree/500/{r}/output/treemerge_node.tree",
        "CINC/IQ": f"{M}/Constrained-INC/IQTree2/500/{r}/node.",
        "GTM/IQ": f"{M}/GTM/500/iqtree/{r}/branch_length.",
        "TreeMerge/IQ": f"{M}/TreeMerge_RAxML-ng/IQTree2/500/{r}/output/treemerge_node.tree",
        "TreeMerge-PAUP/FT": f"{M}/TreeMerge_PAUP/{r}/output/treemerge.tree",
    }


CONDS = [("RNASim1000", rnasim, [str(i) for i in range(1, 6)]),
         ("Cox1-HET", cox, ["R%02d" % i for i in range(1, 11)]),
         ("1000M1-HF", m1, ["R%d" % i for i in range(5)])]

rows = []
for cond, fn, reps in CONDS:
    for r in reps:
        tp, meth = fn(r)
        if not os.path.exists(tp):
            print("missing true tree", cond, r, file=sys.stderr)
            continue
        T = read_tree(tp)
        for m, p in meth.items():
            if not os.path.exists(p):
                print("missing", cond, r, m, p, file=sys.stderr)
                continue
            e = read_tree(p)
            fnr, fpr, nfn, nfp, it, ie = fn_fp(e, T)
            rows.append((cond, r, m, fnr, fpr, nfn, it, ie, len(e.leaves())))

with open(sys.argv[2], "w") as f:
    f.write("condition\treplicate\tmethod\tFN\tFP\tnFN\tinternal_true\tinternal_est\tn_leaves\n")
    for x in rows:
        f.write("%s\t%s\t%s\t%.5f\t%.5f\t%d\t%d\t%d\t%d\n" % x)

# summary
print("%-12s %-18s %8s %8s %5s" % ("condition", "method", "pub FN%", "ours FN%", "nrep"))
for cond, _, _ in CONDS:
    for m in PUB[cond]:
        v = [x[3] for x in rows if x[0] == cond and x[2] == m]
        if v:
            print("%-12s %-18s %8.1f %8.2f %5d" % (cond, m, PUB[cond][m], 100 * sum(v) / len(v), len(v)))
