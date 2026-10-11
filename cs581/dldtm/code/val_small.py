"""Small-tree validation on an *estimated* alignment: take the 2 largest subsets of a
centroid decomposition (<= 50 taxa) of the TRUE tree, realign their raw sequences with
MAFFT L-INS-i, and estimate trees with each method. (On the true alignment the same
comparison comes for free from the pipeline's subset trees: run_rep.py scores every
subset tree against the induced true tree.)
Usage: python3 val_small.py COND REP OUTROOT METHODS
"""
import json
import os
import subprocess
import sys
import time

import common
from phylo import read_tree, fn_fp
import run_rep_methods as rm

cond, rep, outroot, methods = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4].split(",")
d = f"{outroot}/{cond}/{rep}"
os.makedirs(d, exist_ok=True)
aln, truef = common.paths(cond, rep)
T = read_tree(truef)
seqs = common.read_fasta(aln)
subs = sorted(common.centroid_decomp(T, 50), key=len, reverse=True)[:2]
rows = []
for i, s in enumerate(subs):
    raw = f"{d}/s{i}.raw.fa"
    est = f"{d}/s{i}.mafft.fa"
    if not os.path.exists(est):
        with open(raw, "w") as f:
            for n in s:
                f.write(f">{n}\n{seqs[n].replace('-', '')}\n")
        with open(est + ".tmp", "w") as f:
            subprocess.run(["mafft", "--localpair", "--maxiterate", "1000", "--thread", "1", "--quiet", raw],
                           stdout=f, check=True)
        os.replace(est + ".tmp", est)
        # upper-case for downstream tools
        e = common.read_fasta(est)
        common.write_fasta(e, s, est, strip_gap_cols=True)
    for m in methods:
        out = f"{d}/{m}_s{i}.tre"
        t0 = time.time()
        if not os.path.exists(out):
            rm.sub_method(m)(est, out)
        dt = time.time() - t0
        fn, fp = fn_fp(read_tree(out), T.copy().restrict(s))[:2]
        rows.append(dict(kind="subset_est", cond=cond, rep=rep, subset=i, ntaxa=len(s), method=m, FN=fn, FP=fp,
                         seconds=dt))
        print(cond, rep, i, len(s), m, "%.3f" % fn, flush=True)
with open(f"{d}/scores_est.jsonl", "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")
