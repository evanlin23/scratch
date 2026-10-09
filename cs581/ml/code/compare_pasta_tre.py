"""Reproduction check: PASTA's published final tree (pasta.tre, FastTree -gtr -gamma -fastest inside
PASTA) vs our FastTree runs on the published PASTA alignment, per replicate.

    python compare_pasta_tre.py RESULTS.jsonl
"""
import json
import os
import sys

import dendropy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import treeerr  # noqa: E402

MLDATA = os.environ.get("MLDATA", "/opt/data/mlcache")
rows = {(r["dataset"], str(r["rep"]), r["method"]): r for r in map(json.loads, open(sys.argv[1]))
        if r["aln"] == "pasta_align"}
print("| dataset | rep | published pasta.tre FN | our FastTree -fastest FN | our FastTree (default) FN | RF(pasta.tre, ours -fastest) |")
print("|---|---|---|---|---|---|")
for ds in ("1000M2", "1000M3", "1000L1", "RNASim"):
    for rep in ("0", "1", "2"):
        d = os.path.join(MLDATA, ds, "R" + rep)
        pub = treeerr.error(os.path.join(d, "true_tree.tre"), os.path.join(d, "pasta.tre"))["fn_rate"]
        f = rows.get((ds, rep, "fasttree_pasta"))
        g = rows.get((ds, rep, "fasttree"))
        rf = treeerr.error(os.path.join(d, "pasta.tre"), f["tree"])["rf_rate"] if f else None
        print("| %s | %s | %.1f%% | %s | %s | %s |" % (ds, rep, 100 * pub, "%.1f%%" % (100 * f["fn_rate"]) if f else "-",
              "%.1f%%" % (100 * g["fn_rate"]) if g else "-", "%.1f%%" % (100 * rf) if rf is not None else "-"))
