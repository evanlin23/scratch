#!/usr/bin/env python3
"""Baseline validation on Nute et al. 2018 data (doi:10.13012/B2IDB-7735354_V1), 'full' (no missing data).

For every condition/replicate with a published ASTRAL log:
  1. rescore the published ASTRAL tree with qtool against the RAxML gene trees used
     (first k of the 1000 for k-gene conditions) and compare with the normalized
     quartet score printed in the published ASTRAL log;
  2. RF (FN rate) of published ASTRAL / ASTRID / MP-EST / SVDquartets trees vs the true tree.
Writes results/validation_published.tsv and prints condition means.
"""
import glob, os, re, subprocess, sys, tempfile
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from treeutil import bipartitions, rf, read_newick_list  # noqa: E402

ROOT = "/opt/data/miss"
QT = "/opt/runs/qtool"
OUT = os.path.join(HERE, "..", "results", "validation_published.tsv")

rows = []
for sd in sorted(glob.glob(f"{ROOT}/25tax-*gen-0bps-*-full-0")):
    cond = os.path.basename(sd)
    ng = int(re.search(r"25tax-(\d+)gen", cond).group(1))
    gd = sd.replace(f"25tax-{ng}gen", "25tax-1000gen")[:-2]  # gene trees live in the 1000-gene dir
    for rd in sorted(glob.glob(f"{sd}/[0-9][0-9]")):
        rep = os.path.basename(rd)
        gfile = f"{gd}/{rep}/raxml-genes.tre"
        log = f"{rd}/astral-raxml-genes.log"
        if not (os.path.exists(gfile) and os.path.exists(log)):
            continue
        m = re.findall(r"Normalized quartet score is: ([0-9.]+)", open(log).read())
        pub = float(m[-1]) if m else float("nan")
        genes = read_newick_list(gfile)[:ng]
        with tempfile.NamedTemporaryFile("w", suffix=".tre", delete=False) as f:
            f.write("\n".join(genes) + "\n")
        out = subprocess.run([QT, "score", "-g", f.name, "-t", f"{rd}/astral-raxml-genes.tre"],
                             capture_output=True, text=True).stdout.split()
        os.unlink(f.name)
        tb = bipartitions(read_newick_list(f"{rd}/true-species.tre")[0])
        r = dict(cond=cond, rep=rep, ngenes=ng, pub_nscore=pub, our_nscore=float(out[1]))
        for meth in ["astral-raxml", "astrid-raxml", "mpest-raxml", "astral-true", "astrid-true", "svdquartets-true"]:
            fn_ = f"{rd}/{meth}-genes.tre"
            if os.path.exists(fn_):
                fn, fp, n = rf(tb, bipartitions(read_newick_list(fn_)[0]))
                r[meth] = fn / n
        rows.append(r)

cols = ["cond", "rep", "ngenes", "pub_nscore", "our_nscore", "astral-raxml", "astrid-raxml", "mpest-raxml",
        "astral-true", "astrid-true", "svdquartets-true"]
with open(OUT, "w") as f:
    f.write("\t".join(cols) + "\n")
    for r in rows:
        f.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
agg = defaultdict(list)
for r in rows:
    agg[r["cond"]].append(r)
print("cond\tn\tmax|pub-our|\t" + "\t".join(cols[5:]))
for c, rs in agg.items():
    dev = max(abs(r["pub_nscore"] - r["our_nscore"]) for r in rs)
    means = [sum(r.get(m, 0) for r in rs) / len(rs) for m in cols[5:]]
    print(f"{c}\t{len(rs)}\t{dev:.2e}\t" + "\t".join(f"{x:.3f}" for x in means))
