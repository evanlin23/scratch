"""Exhaustive noise-free identifiability enumeration.

For each topology (all for n<=6, a random sample for n=7,8), each true k in 1..3, a set of true
edge sets S (all, or a sample), and several parameter regimes, compute d exactly and find every
edge set S' with |S'| <= K that reproduces d exactly with a strictly valid solution
(p>0, 0<x<1, ybar>=0). Writes one JSON line per instance. Restartable (one file per topology).
"""
import json, os, sys, itertools, random
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from pdd import Tree, all_topologies, Enumerator

OUT = sys.argv[1] if len(sys.argv) > 1 else "enum_out"
REGIMES = ["generic", "generic_y0", "sym", "sym_y0", "ultra"]


def draw(regime, k, rng, tree, S):
    if regime.startswith("generic"):
        p = rng.dirichlet(np.ones(k))
        while p.min() < 0.05:
            p = rng.dirichlet(np.ones(k))
        x = rng.uniform(0.05, 0.95, k)
        y = 0.0 if regime.endswith("y0") else rng.uniform(0.05, 0.5)
    elif regime.startswith("sym"):
        p = np.ones(k) / k
        x = np.full(k, 0.5)
        y = 0.0 if regime.endswith("y0") else 0.5
    else:  # "ultra": unit lengths, random rational-ish p and x
        p = rng.integers(1, 4, k).astype(float); p /= p.sum()
        x = rng.choice([0.25, 0.5, 0.75], k)
        y = rng.choice([0.0, 0.25, 0.5])
    return p, x, y


def job(args):
    n, ti, edges_ab, max_true, K, seed = args
    fn = os.path.join(OUT, f"n{n}_t{ti}.jsonl")
    if os.path.exists(fn):
        return fn
    rng = np.random.default_rng(seed)
    lens = rng.exponential(1.0, len(edges_ab)) + 0.05
    trees = {"generic": Tree(n, [(a, b, l) for (a, b), l in zip(edges_ab, lens)]),
             "unit": Tree(n, [(a, b, 1.0) for (a, b) in edges_ab])}
    enums = {key: Enumerator(t, K) for key, t in trees.items()}
    lines = []
    for k in (1, 2, 3):
        if k > K:
            continue
        subs = list(itertools.combinations(range(trees["unit"].m), k))
        if len(subs) > max_true:
            idx = rng.choice(len(subs), max_true, replace=False)
            subs = [subs[i] for i in idx]
        for S in subs:
            for regime in REGIMES:
                tkey = "generic" if regime.startswith("generic") else "unit"
                t, E = trees[tkey], enums[tkey]
                p, x, y = draw(regime, k, rng, t, S)
                d = t.mixture(S, p, x, y)
                sols = E.solutions(d, Kmax=min(K, k + 1))
                rec = dict(n=n, topo=ti, k=k, regime=regime, S=list(S), p=p.tolist(), x=list(map(float, x)),
                           y=float(y), true_found=False, true_full_rank=None, alts=[])
                for s in sols:
                    if tuple(s["S"]) == tuple(S):
                        rec["true_found"] = True
                        rec["true_full_rank"] = bool(s["full_rank"])
                    else:
                        z = s["z"]; kk = len(s["S"])
                        rec["alts"].append(dict(S=list(s["S"]), kp=kk, p=z[:kk].round(6).tolist(),
                                                w=z[kk:2 * kk].round(6).tolist(), y=round(float(z[2 * kk]), 6),
                                                full_rank=bool(s["full_rank"]), slack=float(s["slack"])))
                lines.append(rec)
    meta = dict(n=n, topo=ti, edges=[list(e) for e in trees["generic"].edges],
                newick_unit=trees["unit"].newick())
    tmp = fn + ".tmp"
    with open(tmp, "w") as f:
        f.write(json.dumps(dict(meta=meta)) + "\n")
        for r in lines:
            f.write(json.dumps(r) + "\n")
    os.replace(tmp, fn)
    return fn


def main():
    os.makedirs(OUT, exist_ok=True)
    jobs = []
    plan = {4: (None, 10**6, 5), 5: (None, 10**6, 5), 6: (None, 40, 5), 7: (40, 25, 4), 8: (40, 15, 4)}
    for n, (ntop, max_true, K) in plan.items():
        tops = all_topologies(n)
        rnd = random.Random(n)
        idx = list(range(len(tops))) if ntop is None else rnd.sample(range(len(tops)), ntop)
        for ti in idx:
            jobs.append((n, ti, tops[ti], max_true, K, 1000 * n + ti))
    with Pool(4) as pool:
        for i, fn in enumerate(pool.imap_unordered(job, jobs)):
            if i % 20 == 0:
                print(i, len(jobs), fn, flush=True)


if __name__ == "__main__":
    main()
