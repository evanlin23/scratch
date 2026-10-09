"""Baseline validation on Nute et al. 2018 data (doi:10.13012/B2IDB-7735354_V1).

1. Table S1: average distance (AD), gene-tree error (GTEE), total discord on 1000-gene full data.
2. Rescore published ASTRAL / ASTRID / MP-EST / SVDquartets species trees (FN rate).
3. Re-run our ASTRID re-implementation (ASTRID v1 distance = edges in rooted gene trees,
   FastME) and check topology agreement with the published ASTRID trees.
Output: results/baseline_published.csv, results/baseline_tableS1.csv
"""
import sys, os, csv, itertools
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import stree

ROOT = "/opt/data/nute"
OUT = os.path.join(os.path.dirname(__file__), "..", "results")
HEIGHTS = ["10M", "2M", "500K"]
RATES = ["1E-7", "1E-6"]


def gdir(h, r, pat):
    return f"{ROOT}/25tax-1000gen-0bps-{h}-{r}-{pat}"


def tableS1(args):
    h, r, rep = args
    d = f"{gdir(h, r, 'full')}/{rep:02d}"
    true = open(f"{d}/true-species.tre").read().strip()
    tg = stree.read_trees(f"{d}/true-genes.tre")
    eg = stree.read_trees(f"{d}/raxml-genes.tre")
    ad = np.mean([stree.fn_fp(true, g)[0] for g in tg])
    gte = np.mean([stree.fn_fp(a, b)[0] for a, b in zip(tg, eg)])
    td = np.mean([stree.fn_fp(true, g)[0] for g in eg])
    return dict(height=h, rate=r, rep=rep, AD=ad, GTEE=gte, TD=td)


def published(args):
    h, r, pat, ng, rep = args
    d = f"{ROOT}/25tax-{ng}gen-0bps-{h}-{r}-{pat}/{rep:02d}"
    rows = []
    true = open(f"{d}/true-species.tre").read().strip()
    for m in ["astral", "astrid", "mpest", "svdquartets"]:
        for gt in ["true", "raxml"]:
            f = f"{d}/{m}-{gt}-genes.tre"
            if not os.path.exists(f) or os.path.getsize(f) == 0:
                continue
            est = open(f).read().strip().split("\n")[-1]
            try:
                fn, fp, _ = stree.fn_fp(true, est)
            except Exception as e:
                continue
            rows.append(dict(height=h, rate=r, pattern=pat, ngenes=ng, rep=rep, method=m,
                             genetrees=gt, FN=fn, FP=fp))
    # our ASTRID v1 re-implementation on the same genes
    gd = f"{gdir(h, r, 'full' if pat == 'full-0' else pat)}/{rep:02d}"
    for gt, fname in [("true", "true-genes.tre" if pat == "full-0" else "true-genes-w-missing-data-na.tre"),
                      ("raxml", "raxml-genes.tre")]:
        genes = stree.read_trees(f"{gd}/{fname}", limit=ng)
        taxa = sorted({l for g in genes for l in stree.parse(g).labels(internal=False)})
        M, C = stree.distance_matrix(genes, taxa, count_root=True)
        ours = stree.fastme(M, taxa)
        pubf = f"{d}/astrid-{gt}-genes.tre"
        pub = open(pubf).read().strip().split("\n")[-1] if os.path.exists(pubf) else None
        same = stree.fn_fp(pub, ours)[0] == 0 if pub else None
        rows.append(dict(height=h, rate=r, pattern=pat, ngenes=ng, rep=rep, method="astrid-ours",
                         genetrees=gt, FN=stree.fn_fp(true, ours)[0], FP=stree.fn_fp(true, ours)[1],
                         same_as_published=same))
    return rows


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    with Pool(4) as pool:
        s1 = pool.map(tableS1, [(h, r, rep) for h in HEIGHTS for r in RATES for rep in range(1, 21)])
        with open(f"{OUT}/baseline_tableS1.csv", "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(s1[0]))
            w.writeheader(); w.writerows(s1)
        jobs = [(h, r, pat, ng, rep) for h in HEIGHTS for r in RATES for pat in ["full-0", "rand-30", "rand-60"]
                for ng in [50, 200, 1000] for rep in range(1, 21)]
        rows = [x for rs in pool.map(published, jobs, chunksize=4) for x in rs]
    keys = ["height", "rate", "pattern", "ngenes", "rep", "method", "genetrees", "FN", "FP", "same_as_published"]
    with open(f"{OUT}/baseline_published.csv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader(); w.writerows(rows)
    print("done", len(rows))
