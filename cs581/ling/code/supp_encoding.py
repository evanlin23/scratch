"""Supplement: is ML-binary's advantage the binary ENCODING or the polymorphism INFORMATION it keeps?

Regenerates the exact TEST replicates of run_sim.py (same seeds; checked by re-scoring NJ) and adds
  MP-bin      : Fitch parsimony on the binary presence/absence encoding (polymorphism kept)
  MP-bin-res  : same, after random resolution of polymorphism (what the multistate methods saw)
  ML-bin-res  : IQ-TREE GTR2+FO+G4 on the binary encoding of the randomly resolved data
Writes ../results/supp_runs.jsonl (restartable).
"""
import json
import os
import sys
import zlib
from multiprocessing import Pool

import numpy as np

import fitch
import methods
import run_sim
import sim
from trees import fn_fp

OUT = os.path.join(run_sim.HERE, '..', 'results', 'supp_runs.jsonl')


def bin_sets(chars, n):
    return [[frozenset([v]) for v in col] for col in methods.binary_matrix(chars, n)]


def run_one(job):
    cond, rep = job
    seed = zlib.crc32(f'{cond}-{rep}'.encode())
    rng = np.random.default_rng(seed)
    params = dict(run_sim.BASE); params.update(run_sim.CONDITIONS[cond])
    t, chars, types, info = sim.simulate(params, rng)
    n = t.n
    rc = sim.resolve_random(chars, rng)
    res = dict(cond=cond, rep=rep, fn={})
    res['fn']['NJ-check'] = round(fn_fp(t, methods.nj_tree(rc, n))[0], 4)
    srng = np.random.default_rng(seed + 7)
    for name, data in [('MP-bin', chars), ('MP-bin-res', rc)]:
        b = bin_sets(data, n)
        est, _ = methods.char_search(b, ['L'] * len(b), n, 'mp', rng=srng, nstarts=8,
                                     start_trees=[methods.nj_tree(rc, n)])
        res['fn'][name] = round(fn_fp(t, est)[0], 4)
    res['fn']['ML-bin-res'] = round(fn_fp(t, methods.iqtree_binary(rc, n, seed=seed % 100000 + 1))[0], 4)
    return res


def main():
    done = set()
    if os.path.exists(OUT):
        done = {(d['cond'], d['rep']) for d in map(json.loads, open(OUT))}
    jobs = [(c, r) for r in range(10, 30) for c in run_sim.CONDITIONS if (c, r) not in done]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as pool:
        for res in pool.imap_unordered(run_one, jobs):
            with open(OUT, 'a') as f:
                f.write(json.dumps(res) + '\n')


if __name__ == '__main__':
    main()
