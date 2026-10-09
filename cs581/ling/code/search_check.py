"""Search adequacy: does the SPR hill-climb find trees at least as good as the TRUE tree under
its own objective? (If the true tree scores better, the search failed.)"""
import numpy as np, sim, fitch, methods, run_sim
out = []
for kind, cap, W in [('mp', None, None), ('mc', None, None), ('cap', 2, None), ('mc', None, run_sim.WEIGHTS)]:
    fails = 0; tot = 0
    for rep in range(12):
        rng = np.random.default_rng(1000 + rep)
        p = dict(run_sim.BASE); p.update(run_sim.CONDITIONS[['moderate', 'borrow3', 'homoplasy'][rep % 3]])
        t, ch, ty, _ = sim.simulate(p, rng)
        rc = sim.resolve_random(ch, rng)
        ls, lb = fitch.encode_sets(rc, 24)
        w = np.ones(len(rc)) if W is None else np.array([W[x] for x in ty])
        obj = fitch.Objective(kind, cap=cap)
        est, sc = methods.char_search(rc, ty, 24, kind, cap=cap, weights=W, rng=rng, nstarts=8,
                                      start_trees=[methods.nj_tree(rc, 24)])
        st, _ = fitch.score_tree(t, ls, w, lb, obj)
        tot += 1; fails += sc > st + 1e-9
    out.append(f'{kind}{"" if cap is None else cap}{" weighted" if W else ""}: search worse than true tree in {fails}/{tot}')
print('\n'.join(out))
open('../results/search_check.txt', 'w').write('\n'.join(out) + '\n')
