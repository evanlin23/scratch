"""One replicate of the DL-subset-tree + GTM pipeline (restartable; every step is cached).

  1. FastTree (GTR) on the full true alignment -> guide tree (also the "FastTree full" baseline)
  2. IQ-TREE 3 (GTR+G, --fast) on the full alignment -> "IQ-TREE full" baseline
  3. centroid decomposition of the guide into subsets of <= MAXSUB taxa
  4. subset trees with each method: FT, IQ (GTR+G --fast), BME (FastME BME+SPR on JC69
     distances, the non-DL distance control), PF (pretrained Phyloformer -> FastME BME+SPR),
     NNJ (NeuralNJ argmax), ...
  5. GTM (github.com/vlasmirnov/GTM, default mode) merges each method's subset trees with the guide
  6. FN/FP of every full tree vs the true tree, and of every subset tree vs the induced true tree

Usage: python3 run_rep.py COND REP OUTROOT MAXSUB METHODS(comma list)
Writes OUTROOT/COND/REP/... and appends JSON lines to OUTROOT/COND/REP/scores.jsonl.
Single-threaded on purpose: four replicates run in parallel on the 4-core machine.
"""
import json
import os
import subprocess
import sys
import time

import common
from phylo import read_tree, fn_fp

cond, rep, outroot, maxsub, methods = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5].split(",")
os.environ.setdefault("OMP_NUM_THREADS", "1")
d = f"{outroot}/{cond}/{rep}"
os.makedirs(d, exist_ok=True)
aln, truef = common.paths(cond, rep)
T = read_tree(truef)
timing_f = f"{d}/timing.json"
timing = json.load(open(timing_f)) if os.path.exists(timing_f) else {}


def timed(key, out, fn):
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return
    t0 = time.time()
    fn(out + ".tmp")
    os.replace(out + ".tmp", out)
    timing[key] = time.time() - t0
    json.dump(timing, open(timing_f, "w"), indent=1)
    print(f"[{cond} {rep}] {key}: {timing[key]:.1f}s", flush=True)


# full alignment in our naming (U->T)
full = f"{d}/full.fa"
if not os.path.exists(full):
    s = common.read_fasta(aln)
    common.write_fasta(s, list(s), full, strip_gap_cols=False)
seqs = common.read_fasta(full)

timed("FT_full", f"{d}/FT_full.tre", lambda o: common.run_fasttree(full, o))
if "IQfull" in methods:
    timed("IQ_full", f"{d}/IQ_full.tre", lambda o: common.run_iqtree(full, o, threads=1, fast=True))
guide = f"{d}/FT_full.tre"

sd = f"{d}/sub{maxsub}"
os.makedirs(sd, exist_ok=True)
subs_f = f"{sd}/subsets.json"
if not os.path.exists(subs_f):
    json.dump(common.centroid_decomp(read_tree(guide), maxsub), open(subs_f, "w"))
subsets = json.load(open(subs_f))
for i, s in enumerate(subsets):
    if not os.path.exists(f"{sd}/s{i}.fa"):
        common.write_fasta(seqs, s, f"{sd}/s{i}.fa")


from run_rep_methods import sub_method  # noqa: E402

sub_methods = [m for m in methods if m != "IQfull"]
for m in sub_methods:
    f = sub_method(m)
    for i in range(len(subsets)):
        timed(f"{m}_s{i}", f"{sd}/{m}_s{i}.tre", lambda o, i=i: f(f"{sd}/s{i}.fa", o))
    trees = [f"{sd}/{m}_s{i}.tre" for i in range(len(subsets))]
    timed(f"GTM_{m}", f"{sd}/GTM_{m}.tre",
          lambda o: subprocess.run([sys.executable, f"{common.GTM}/gtm.py", "-s", guide, "-t", *trees, "-o", o],
                                   check=True, capture_output=True))

# scores
rows = []
fulls = [("FastTree", "FT_full", f"{d}/FT_full.tre"), ("IQ-TREE", "IQ_full", f"{d}/IQ_full.tre")]
fulls += [(f"GTM+{m}", f"GTM_{m}", f"{sd}/GTM_{m}.tre") for m in sub_methods]
for name, key, p in fulls:
    if not os.path.exists(p):
        continue
    fn, fp = fn_fp(read_tree(p), T)[:2]
    t = timing.get(key, 0.0)
    if key.startswith("GTM_"):
        m = key[4:]
        t_sub = sum(v for k, v in timing.items() if k.startswith(f"{m}_s"))
        t = timing.get("FT_full", 0) + t_sub + timing.get(key, 0)
    rows.append(dict(kind="full", cond=cond, rep=rep, maxsub=maxsub, method=name, FN=fn, FP=fp, seconds=t))
for m in sub_methods:
    for i, s in enumerate(subsets):
        if len(s) < 4:
            continue
        fn, fp = fn_fp(read_tree(f"{sd}/{m}_s{i}.tre"), T.copy().restrict(s))[:2]
        rows.append(dict(kind="subset", cond=cond, rep=rep, maxsub=maxsub, method=m, subset=i, ntaxa=len(s),
                         FN=fn, FP=fp, seconds=timing.get(f"{m}_s{i}", 0.0)))
with open(f"{d}/scores_sub{maxsub}.jsonl", "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")
for r in rows:
    if r["kind"] == "full":
        print(cond, rep, maxsub, r["method"], "FN=%.4f" % r["FN"], "t=%.0fs" % r["seconds"])
