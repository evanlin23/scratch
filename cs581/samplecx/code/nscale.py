"""Focused n-scaling test (Roch 2018: NJst/ASTRID may need m >= linear in n; ASTRAL needs O(log n)).
Caterpillar and balanced trees with every internal branch f = 0.2 CU, n = 8..128, fine k grid,
true gene trees. ASTRAL is run only up to k <= 600 (its k_95 at f = 0.2 is ~200).
Output: /opt/runs/nscale.jsonl (restartable).  Usage: python nscale.py reps nproc
"""
import sys, os, json, time
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
import msc, stree

OUT = "/opt/runs/nscale.jsonl"
F = 0.2
KG = [25, 50, 75, 100, 150, 200, 300, 400, 600, 800, 1200, 1600, 2400]
ACAP = 600


def run(job):
    shape, n, rep = job
    sp = msc.caterpillar(n, F) if shape == "cat" else msc.balanced(n, F)
    true = msc.species_newick(sp)
    tb, _ = stree.bipartitions(true)
    genes = msc.sim_genes(sp, KG[-1], 7777777 + 1000 * n + rep + (0 if shape == "cat" else 500000))
    taxa = [f"t{i}" for i in range(n)]
    acc = stree.Accumulator(taxa)
    rows, done = [], 0
    for k in KG:
        for g in genes[done:k]:
            acc.add(g)
        done = k
        M = acc.matrix()
        ests = {"astrid": stree.fastme(M, taxa), "njst": stree.nj(M, taxa)}
        if k <= ACAP:
            ests["astral"] = stree.astral(genes[:k], ("-r", "1", "-s", "0"))
        for m, e in ests.items():
            eb, _ = stree.bipartitions(e)
            rows.append(dict(shape=shape, n=n, f=F, rep=rep, k=k, method=m, FN=len(tb - eb) / len(tb)))
    return rows


if __name__ == "__main__":
    reps, nproc = int(sys.argv[1]), int(sys.argv[2])
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            d = json.loads(l); done.add((d["shape"], d["n"], d["rep"]))
    jobs = [(s, n, r) for r in range(reps) for s in ["cat", "bal"] for n in [8, 16, 32, 64, 128]
            if (s, n, r) not in done]
    print(len(jobs), flush=True)
    with Pool(nproc) as p, open(OUT, "a") as fh:
        for rows in p.imap_unordered(run, jobs):
            for x in rows: fh.write(json.dumps(x) + "\n")
            fh.flush()
