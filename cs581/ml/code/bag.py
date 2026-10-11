"""Pilot: bootstrap-aggregated ("bagged") ML trees for low-signal alignments.

    python bag.py DATASET REP --method fasttree|iqtree_fast|raxmlng_fastmode [--B 50] [--out OUT.jsonl]

Resample alignment columns with replacement B times, estimate a tree on each replicate with
METHOD, and summarise with the greedy (extended majority-rule, MRE) consensus, which is fully
resolved when possible -- so FN is comparable with single-tree methods. Reports FN/FP/RF of the
consensus vs the true tree, total CPU, and lnL of the consensus under GTR+G (raxml-ng --evaluate).
Model-agnostic: only the substitution model string is data-type specific.
"""
import argparse
import json
import os
import random
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evaltree  # noqa: E402
import runtrees as rt  # noqa: E402
import treeerr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("rep")
    ap.add_argument("--method", default="fasttree")
    ap.add_argument("--B", type=int, default=50)
    ap.add_argument("--out", default=os.path.join(HERE, "..", "results", "bag.jsonl"))
    a = ap.parse_args()
    d = os.path.join(rt.MLDATA, a.dataset, "R%s" % a.rep)
    clean = os.path.join(d, "true_align.clean.fasta")
    if not os.path.exists(clean):
        rt.clean_alignment(os.path.join(d, "true_align.fasta"), clean + ".tmp")
        os.replace(clean + ".tmp", clean)
    names, seqs = rt.read_fasta(clean)
    L = len(seqs[0])
    work = tempfile.mkdtemp(prefix="bag_%s_%s_%s_" % (a.dataset, a.rep, a.method))
    rng = random.Random(1)
    cpu, trees = 0.0, []
    for b in range(a.B):
        cols = [rng.randrange(L) for _ in range(L)]
        ba = os.path.join(work, "b%d.fasta" % b)
        rt.write_fasta(ba, names, ["".join(s[i] for i in cols) for s in seqs])
        bt = os.path.join(work, "b%d.tre" % b)
        wb = os.path.join(work, "w%d" % b)
        os.makedirs(wb)
        rt.estimate(a.method, ba, bt, wb)
        cpu += rt.CPU.seconds
        evaltree.binary(bt, bt + ".bin")  # FastTree may collapse identical sequences into polytomies
        trees.append(open(bt + ".bin").read().strip())
    allt = os.path.join(work, "all.nwk")
    open(allt, "w").write("\n".join(trees) + "\n")
    subprocess.run([rt.RAXMLNG, "--consense", "MRE", "--tree", allt, "--prefix", os.path.join(work, "cons"),
                    "--threads", "1"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    cons = os.path.join(d, "trees", "true_align.bag%d_%s.tre" % (a.B, a.method))
    os.makedirs(os.path.dirname(cons), exist_ok=True)
    os.replace(os.path.join(work, "cons.raxml.consensusTreeMRE"), cons)
    true = os.path.join(d, "true_tree.tre")
    e = treeerr.error(true, cons)
    row = {"dataset": a.dataset, "rep": a.rep, "aln": "true_align", "method": "bag%d_%s" % (a.B, a.method),
           "cpu_seconds": round(cpu, 1), "fn_rate": e["fn_rate"], "fp_rate": e["fp_rate"], "rf_rate": e["rf_rate"],
           "est_internal": e["est_int"], "lnl_eval": evaltree.evaluate(clean, cons), "tree": cons}
    with open(a.out, "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
