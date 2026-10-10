"""Locate which species-tree branch ASTRAL-Pro gets wrong under a tagging-error model
(quartet supports from -C on the true tree). usage: python loc_fail.py SETTING REP NFAM model:par ..."""
import sys, os, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import core
from phylo import parse_newick
setting, rep, nfam = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
d = '/opt/runs/gdlcons/data/%s/rep%d/' % (setting, rep)
spn = open(d + 'species.nwk').read().strip(); print(spn)
sp = parse_newick(spn); species = sorted(sp.label[v] for v in sp.leaves()); spi = {s: i for i, s in enumerate(species)}
lines = [l for l in open(d + 'fams.nwk')][:nfam]
base = [core.true_rt(l, spi) for l in lines]
rng = random.Random(1)
for mp in sys.argv[4:]:
    model, par = mp.split(':'); par = float(par)
    err = [core.apply_error(model, par, rt, tg, spi, rng) for rt, tg in base]
    nw = [core.encode(r, t) for r, t in err]
    est = core.run_apro(nw, [r for r, _ in err], True)
    print(model, par, 'FN', core.rf_error(est, sp)[0], est.newick())
    sc = core.quartet_scores(nw, [r for r, _ in err], True, sp.newick())
    per = {}
    for node, t, lab, s in sc: per.setdefault(node, []).append((t, lab, s))
    for node, v in per.items():
        s1 = [x[2] for x in v if x[0] == 't1'][0]
        mx = max(x[2] for x in v if x[0] != 't1')
        if s1 <= mx * 1.05: print('  CLOSE/FAIL', node, v)
