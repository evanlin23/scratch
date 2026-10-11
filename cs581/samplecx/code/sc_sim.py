"""Controlled sample-complexity sweep with true gene trees.

Species trees: caterpillar / balanced, n taxa, EVERY internal branch = f coalescent units.
For each replicate: simulate k_max gene trees under the MSC, then for nested prefixes of
k genes run ASTRID (FastME BME+SPR on mean internode distance), NJst (NJ on the same matrix)
and ASTRAL (ASTER, -r 1 -s 0). Record which true bipartitions are missed.
Output: /opt/runs/sc_sim.jsonl (one row per (job, k, method)); restartable.
Usage: python sc_sim.py [reps_from] [reps_to]
"""
import sys, os, json, time
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import msc, stree

OUT = "/opt/runs/sc_sim.jsonl"
KGRID = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000]
FS = [0.05, 0.1, 0.2, 0.5]
NS = [8, 16, 32, 64]
SHAPES = ["cat", "bal"]
ASTRAL_ARGS = ("-r", "1", "-s", "0")


ACAP = {8: 5000, 16: 5000, 32: 2000, 64: 1000}   # ASTRAL is superlinear in k: cap its grid
DCAP = {8: 5000, 16: 5000, 32: 5000, 64: 2000}   # distance methods


def kgrid(f, n=8):
    return [k for k in KGRID if k <= 60 / f ** 2 and k <= DCAP[n]]


def run(job):
    shape, n, f, rep = job
    sp = msc.caterpillar(n, f) if shape == "cat" else msc.balanced(n, f)
    true = msc.species_newick(sp)
    tb, _ = stree.bipartitions(true)
    tb = sorted(tb, key=lambda s: (len(s), sorted(s)))
    seed = int(1000003 * n + 7919 * round(f * 1000) + 31 * rep + (0 if shape == "cat" else 17))
    ks = kgrid(f, n)
    genes = msc.sim_genes(sp, ks[-1], seed)
    taxa = [f"t{i}" for i in range(n)]
    acc = stree.Accumulator(taxa)
    rows = []
    done = 0
    for k in ks:
        for g in genes[done:k]:
            acc.add(g)
        done = k
        M = acc.matrix()
        ests = {}
        t0 = time.time(); ests["astrid"] = stree.fastme(M, taxa); ta = time.time() - t0
        ests["njst"] = stree.nj(M, taxa)
        if k <= ACAP[n]:
            t0 = time.time(); ests["astral"] = stree.astral(genes[:k], ASTRAL_ARGS); tq = time.time() - t0
        for m, e in ests.items():
            eb, _ = stree.bipartitions(e)
            missed = [i for i, b in enumerate(tb) if b not in eb]
            rows.append(dict(shape=shape, n=n, f=f, rep=rep, k=k, method=m,
                             FN=len(missed) / len(tb), missed=missed,
                             sec=ta if m == "astrid" else tq if m == "astral" else None))
    return rows


if __name__ == "__main__":
    r0, r1 = int(sys.argv[1]), int(sys.argv[2])
    donejobs = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            d = json.loads(line)
            donejobs.add((d["shape"], d["n"], d["f"], d["rep"]))
    jobs = [(s, n, f, r) for r in range(r0, r1) for s in SHAPES for n in NS for f in FS
            if (s, n, f, r) not in donejobs]
    # big jobs first within a rep block for load balance
    print(len(jobs), "jobs", flush=True)
    with Pool(4) as pool, open(OUT, "a") as fh:
        for rows in pool.imap_unordered(run, jobs):
            for r in rows:
                fh.write(json.dumps(r) + "\n")
            fh.flush()
            print(rows[0]["shape"], rows[0]["n"], rows[0]["f"], rows[0]["rep"], flush=True)
