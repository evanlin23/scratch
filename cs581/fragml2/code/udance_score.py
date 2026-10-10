"""Score uDance (udance/NOTES.md) on M1HF and compare with our arms on uDance's retained leaf set.

    python udance_score.py > ../results/udance.md

uDance drops some taxa (unplaced queries, TreeShrink, stitching), so every other tree is pruned to the
uDance output's leaf set before scoring (FN/FP on the same leaves). CPU = uDance's summed user+sys from
/opt/udance_work/R*/time_v*.txt (R0: summed over the 4 resumed runs, from NOTES.md) + the FastTree
backbone tree it was given.
"""
import glob
import json
import os
import re
import sys

import dendropy
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import treeerr  # noqa: E402

MLDATA = "/opt/data/mlcache"
ARMS = ["base_fasttree", "base_iqfast", "base_raxmlng", "place_ft_0.5_fix_graft", "constr_ft_0.5",
        "place_ft_0.5_fix_rxfast"]


def leaves(p):
    t = dendropy.Tree.get(path=p, schema="newick", preserve_underscores=True)
    return {x.taxon.label for x in t.leaf_node_iter()}


def pruned_err(true, est, keep, tmp):
    t = dendropy.Tree.get(path=est, schema="newick", preserve_underscores=True)
    t.retain_taxa_with_labels(keep)
    t.write(path=tmp, schema="newick", suppress_rooting=True)
    return treeerr.error(true, tmp)


def cpu_of(rep):
    if rep == 0:
        return 1315.0
    s = 0.0
    for f in glob.glob("/opt/udance_work/R%d/time_v*.txt" % rep):
        txt = open(f).read()
        s += sum(float(x) for x in re.findall(r"(?:User|System) time \(seconds\): (\S+)", txt))
    return s


rows = {}
for l in open(os.path.join(HERE, "..", "results", "runs.jsonl")):
    r = json.loads(l)
    if r["dataset"] == "M1HF":
        rows[(r["rep"], r["arm"])] = r
res = {a: [] for a in ["udance"] + ARMS}
print("| rep | leaves kept by uDance | uDance FN | uDance FP | uDance CPU min | " +
      " | ".join("%s FN" % a for a in ARMS) + " |")
print("|---|---|---|---|---|" + "---|" * len(ARMS))
for rep in range(10):
    u = "/opt/udance_work/R%d/output/udance.updates.nwk" % rep
    if not os.path.exists(u):
        continue
    d = os.path.join(MLDATA, "M1HF", "R%d" % rep)
    true = os.path.join(d, "true_tree.tre")
    keep = leaves(u)
    e = treeerr.error(true, u)
    bb = json.load(open(os.path.join(d, "trees", "cache", "true_align.bb_ft_0.5.cost.json")))["cpu"]
    cpu = (cpu_of(rep) + bb) / 60
    res["udance"].append((e["fn_rate"], cpu))
    cells = []
    for a in ARMS:
        if (rep, a) in rows:
            x = pruned_err(true, rows[(rep, a)]["tree"], keep, "/tmp/_ud_%d.tre" % rep)["fn_rate"]
            res[a].append((x, rep, e["fn_rate"]))
            cells.append("%.1f" % (100 * x))
        else:
            cells.append("–")
    print("| R%d | %d | %.1f | %.1f | %.1f | %s |" % (rep, len(keep), 100 * e["fn_rate"], 100 * e["fp_rate"], cpu,
                                                    " | ".join(cells)))
print("\nMean over replicates where both exist (FN on uDance's leaf set; Δ = arm − uDance, points):\n")
print("| arm | n | mean FN | mean uDance FN | Δ |")
print("|---|---|---|---|---|")
for a in ARMS:
    if res[a]:
        print("| %s | %d | %.1f%% | %.1f%% | %+.1f |" % (a, len(res[a]), 100 * np.mean([x[0] for x in res[a]]),
                                                       100 * np.mean([x[2] for x in res[a]]),
                                                       100 * np.mean([x[0] - x[2] for x in res[a]])))
