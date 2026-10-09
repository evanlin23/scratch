"""Rescore Park, Zaharias & Warnow (2021, Algorithms 14:148) published RNASim1000 trees.

    python validate_park.py [PARK_DIR]   (default /opt/data/park2021/RNASim1000)

Data: Illinois Data Bank doi:10.13012/B2IDB-7008049_V1 (RNASim1000.tar.gz + RNASim1000_Analysis.tar.gz).
Prints FN rate per method and replicate, and the mean (compare with the paper's Table 3), and
links each replicate's true alignment/tree into $MLDATA/ParkRNASim1000/R<rep> for runtrees.py.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import treeerr  # noqa: E402

P = sys.argv[1] if len(sys.argv) > 1 else "/opt/data/park2021/RNASim1000"
MLDATA = os.environ.get("MLDATA", "/opt/data/mlcache")
TREES = {
    "FastTree": "CreateConstraintTrees/500/FastTree/{r}/output/fasttree.out",
    "IQ-TREE 2": "CreateConstraintTrees/500/IQTree2/{r}/output/iqtree-full.treefile",
    "RAxML-NG (24h cap, lastTree)": "RAxML-ng/{r}/raxmlng-result.raxml.lastTree.TMP",
    "GTM (FastTree start, 500)": "GTM/fasttree_fasttree/500/{r}/branch_length.",
    "GTM (IQ-TREE start, 500)": "GTM/IQTree2/500/{r}/branch_length.",
    "GTM (IQ-TREE start, 120)": "GTM/IQTree2/120/{r}/branch_length.",
}
res = {}
for r in range(1, 6):
    true = os.path.join(P, str(r), "model", "true.tt")
    d = os.path.join(MLDATA, "ParkRNASim1000", "R%d" % r)
    os.makedirs(d, exist_ok=True)
    for src, dst in (("true.fasta", "true_align.fasta"), ("true.tt", "true_tree.tre")):
        if not os.path.exists(os.path.join(d, dst)):
            os.symlink(os.path.join(P, str(r), "model", src), os.path.join(d, dst))
    for m, pat in TREES.items():
        p = os.path.join(P, pat.format(r=r))
        if os.path.exists(p) and os.path.getsize(p) > 0:
            res.setdefault(m, []).append(treeerr.error(true, p)["fn_rate"])
print("| method (published tree) | FN per replicate 1-5 | mean FN |")
print("|---|---|---|")
for m, v in res.items():
    print("| %s | %s | %.1f%% |" % (m, " ".join("%.1f" % (100 * x) for x in v), 100 * np.mean(v)))


# ---- 1000M1-HF (Table 5): ROSE 1000M1 R0-R4 with half the sequences fragmented ----------
import glob  # noqa: E402

H = os.path.join(os.path.dirname(P), "1000M1_HF")
ROSE = "/opt/data/Datasets/ROSE/1000M1"
HF = {
    "FastTree": "CreateConstraintTrees/FastTree/500/R{r}/output/fasttree.out",
    "IQ-TREE 2": "IQTree/R{r}/iqtree-result.treefile",
    "RAxML-NG (24h cap, lastTree)": "RAxML-ng/R{r}/raxmlng-result.raxml.lastTree.TMP",
    "GTM (FastTree start, 500)": "GTM/500/fasttree_fasttree/R{r}/branch_length.",
    "GTM (IQ-TREE start, 500)": "GTM/500/iqtree/R{r}/branch_length.",
}
if os.path.isdir(H):
    res = {}
    for r in range(5):
        true = os.path.join(ROSE, "R%d" % r, "rose.tt")
        d = os.path.join(MLDATA, "Park1000M1HF", "R%d" % r)
        os.makedirs(d, exist_ok=True)
        # the full fragmentary alignment = union of the (full-width) subset alignments
        if not os.path.exists(os.path.join(d, "true_align.fasta")):
            with open(os.path.join(d, "true_align.fasta"), "w") as out:
                for f in sorted(glob.glob(os.path.join(H, "CreateConstraintTrees/FastTree/500/R%d/output/sequence_partition_*.out" % r))):
                    out.write(open(f).read().rstrip("\n") + "\n")
            os.symlink(true, os.path.join(d, "true_tree.tre"))
        for m, pat in HF.items():
            p = os.path.join(H, pat.format(r=r))
            if os.path.exists(p) and os.path.getsize(p) > 0:
                res.setdefault(m, []).append(treeerr.error(true, p)["fn_rate"])
    print("\n1000M1-HF (Park et al. 2021 Table 5)\n")
    print("| method (published tree) | FN per replicate 0-4 | mean FN |")
    print("|---|---|---|")
    for m, v in res.items():
        print("| %s | %s | %.1f%% |" % (m, " ".join("%.1f" % (100 * x) for x in v), 100 * np.mean(v)))
