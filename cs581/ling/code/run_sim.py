"""Paired simulation experiment. Restartable: appends one JSON line per (condition, replicate)
to ../results/sim_runs.jsonl and skips pairs already present.

Usage: python3 run_sim.py [nreps] [nproc] [cond1,cond2,...]
Replicates 0..9 are TRAINING (used only to pick the cap of capped parsimony); 10.. are TEST.
"""
import json
import os
import sys
import time
import zlib
from multiprocessing import Pool

import numpy as np

import methods
import sim
from trees import fn_fp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results', 'sim_runs.jsonl')

BASE = dict(rL=3.0, rP=0.25, p_poly=0.2)
CONDITIONS = {
    # name: overrides of sim.DEFAULTS (on top of BASE)
    'clean':     dict(n_contact=0, hom_L=0.0, hom_M=0.0, hom_P=0.0, p_poly=0.0),
    'moderate':  dict(),                                         # 1 contact, moderate homoplasy, poly 0.2
    'borrow3':   dict(n_contact=3, b_L=0.15, b_P=0.05),          # heavy borrowing
    'poly-high': dict(p_poly=0.5),                               # high polymorphism
    'homoplasy': dict(hom_L=0.25, hom_M=0.5, hom_P=0.6),         # high homoplasy (weights misspecified)
    'lexonly':   dict(nM=0, nP=0, n_contact=2),                  # lexical-only data, 2 contacts
}
CAPS = [2, 3, 4]
WEIGHTS = {'L': 1.0, 'M': 5.0, 'P': 5.0}


def run_one(job):
    cond, rep = job
    seed = zlib.crc32(f'{cond}-{rep}'.encode())
    rng = np.random.default_rng(seed)
    params = dict(BASE); params.update(CONDITIONS[cond])
    t, chars, types, info = sim.simulate(params, rng)
    n = t.n
    rc = sim.resolve_random(chars, rng)
    res = dict(cond=cond, rep=rep, seed=seed, split='train' if rep < 10 else 'test',
               summary=sim.summary(chars, types), contacts=len(info['contacts']),
               borrow=info['stats']['borrow'], fn={}, secs={})

    def rec(name, fn, *a, **k):
        t0 = time.time()
        est = fn(*a, **k)
        if isinstance(est, tuple):
            est = est[0]
        res['secs'][name] = round(time.time() - t0, 3)
        res['fn'][name] = round(fn_fp(t, est)[0], 4)
        return est

    njt = rec('NJ', methods.nj_tree, rc, n)
    srng = np.random.default_rng(seed + 1)
    common = dict(rng=srng, nstarts=8, start_trees=[njt])
    rec('MP', methods.char_search, rc, types, n, 'mp', **common)
    rec('MP-poly', methods.char_search, chars, types, n, 'mp', **common)
    rec('MC', methods.char_search, rc, types, n, 'mc', **common)
    rec('WMC', methods.char_search, rc, types, n, 'mc', weights=WEIGHTS, **common)
    for c in CAPS:
        rec(f'Cap{c}', methods.char_search, rc, types, n, 'cap', cap=c, **common)
        rec(f'WCap{c}', methods.char_search, rc, types, n, 'cap', cap=c, weights=WEIGHTS, **common)
    rec('ML-Mk', methods.iqtree_mk, chars, n, seed=seed % 100000 + 1, rng=np.random.default_rng(seed + 2))
    rec('ML-bin', methods.iqtree_binary, chars, n, seed=seed % 100000 + 1)
    return res


def main():
    nreps = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    conds = sys.argv[3].split(',') if len(sys.argv) > 3 else list(CONDITIONS)
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            d = json.loads(l); done.add((d['cond'], d['rep']))
    jobs = [(c, r) for r in range(nreps) for c in conds if (c, r) not in done]
    print(f'{len(jobs)} jobs', file=sys.stderr)
    with Pool(nproc) as pool:
        for res in pool.imap_unordered(run_one, jobs):
            with open(OUT, 'a') as f:
                f.write(json.dumps(res) + '\n')
            print(res['cond'], res['rep'], res['fn'], file=sys.stderr, flush=True)


if __name__ == '__main__':
    main()
