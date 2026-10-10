"""Prepare true gene trees (+ true tags) aligned gene-by-gene with the estimated trees.
FastMulRFS: /opt/data/prep/fmrfs/<cond>/<rep>/{true.trees,tt.trees} for the first N genes of genes-with-gt3species.txt
DISCO:      /opt/data/prep/disco/<cond>/<rep>/tt.trees (all genes of g_true.trees)
Usage: python prep.py fmrfs N  |  python prep.py disco COND REPS(01-10)"""
import os, subprocess, sys
H = os.path.dirname(os.path.abspath(__file__))
PY = "/opt/mm/root/envs/gdl/bin/python"
if sys.argv[1] == "fmrfs":
    N = int(sys.argv[2])
    for dl in ("0.0000000001", "0.0000000002", "0.0000000005"):
        for ps in ("10000000", "50000000"):
            for r in range(1, 11):
                rep = "%02d" % r
                d = f"/opt/data/fmrfs/ntaxa-100.dlrate-{dl}.psize-{ps}/{rep}"
                o = f"/opt/data/prep/fmrfs/dl{dl[-1]}e-10_ps{ps[0]}e7/{rep}"
                os.makedirs(o, exist_ok=True)
                ids = [int(x) for x in open(d + "/genes-with-gt3species.txt").read().split()][:N]
                open(o + "/ids.txt", "w").write("\n".join(map(str, ids)) + "\n")
                g = open(d + "/g_trees.trees").readlines()
                open(o + "/true.trees", "w").write("".join(g[i - 1] for i in ids))
                subprocess.run([PY, H + "/truetag.py", d + "/s_tree.trees", d + "/l_trees.trees", d + "/g_trees.trees",
                                o + "/tt.trees", "--genes", o + "/ids.txt"], check=True)
else:
    c = sys.argv[2]
    a, b = (int(x) for x in sys.argv[3].split("-"))
    for r in range(a, b + 1):
        rep = "%02d" % r
        d = f"/opt/data/disco/trees/{c}/{rep}"
        o = f"/opt/data/prep/disco/{c}/{rep}"
        os.makedirs(o, exist_ok=True)
        if not os.path.exists(o + "/tt.trees"):
            subprocess.run([PY, H + "/truetag.py", d + "/s_tree.trees", d + "/l_trees.trees", d + "/g_true.trees",
                            o + "/tt.trees.tmp"], check=True)
            os.rename(o + "/tt.trees.tmp", o + "/tt.trees")
