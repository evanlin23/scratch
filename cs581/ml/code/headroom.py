"""Headroom survey: for each (dataset, replicate) with a RAxML-NG tree on the TRUE alignment,
score the true tree and the RAxML-NG tree with the same evaluator (raxml-ng --evaluate, GTR+G)
and report the search gap  lnL(true) - lnL(RAxML-NG)  together with FN errors of all methods.

    python headroom.py OUT.jsonl
Search gap > 0: RAxML-NG missed a better-scoring tree that is also more accurate (search headroom).
Search gap < 0: the ML criterion prefers the (wrong) found tree (no search headroom).
"""
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evaltree  # noqa: E402
import runtrees as rt  # noqa: E402
import treeerr  # noqa: E402

M = rt.MLDATA
PARK_RNA = "/opt/data/park2021/RNASim1000/RAxML-ng/{r}/raxmlng-result.raxml.lastTree.TMP"


def jobs():
    for d in sorted(glob.glob(os.path.join(M, "*", "R*"))):
        ds, rep = d.split("/")[-2], d.split("/")[-1][1:]
        if ds in ("16SMHF", "16SMHFamp", "16S.M", "16S.T"):  # no true tree for real data
            continue
        rx = os.path.join(d, "trees", "true_align.raxmlng.tre")
        if ds == "ParkRNASim1000":
            rx = PARK_RNA.format(r=rep)
        if os.path.exists(rx) and os.path.exists(os.path.join(d, "true_align.fasta")):
            yield ds, rep, d, rx


def main(out):
    done = set()
    if os.path.exists(out):
        done = {(r["dataset"], r["rep"]) for r in map(json.loads, open(out))}
    for ds, rep, d, rx in jobs():
        if (ds, rep) in done:
            continue
        clean = os.path.join(d, "true_align.clean.fasta")
        if not os.path.exists(clean):
            rt.clean_alignment(os.path.join(d, "true_align.fasta"), clean + ".tmp")
            os.replace(clean + ".tmp", clean)
        true = os.path.join(d, "true_tree.tre")
        row = {"dataset": ds, "rep": rep, "lnl_true": evaltree.evaluate(clean, true),
               "lnl_raxmlng": evaltree.evaluate(clean, rx), "fn_raxmlng": treeerr.error(true, rx)["fn_rate"]}
        for m in ("fasttree", "iqtree_fast", "raxmlng_fastmode"):
            p = os.path.join(d, "trees", "true_align.%s.tre" % m)
            if os.path.exists(p):
                row["fn_" + m] = treeerr.error(true, p)["fn_rate"]
                row["lnl_" + m] = evaltree.evaluate(clean, p)
        row["search_gap"] = row["lnl_true"] - row["lnl_raxmlng"]
        with open(out, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
