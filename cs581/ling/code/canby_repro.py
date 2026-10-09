"""MP on the Canby et al. 2024 corrected IE dataset (polymorphic cells 'a/b' kept as state sets).
Canby et al. report two equally optimal (weighted) MP trees, one matching Nakhleh et al. 2005.
We run weighted MP (their weight column: 10000 for the 30 '!' characters, 1 otherwise) from many
starts and list the distinct optimal trees found and their subgroups."""
import csv, os
import numpy as np
import fitch, methods
from ie_repro import CLADES, RWT_NWK, subgroup_hits
from trees import from_nested, parse_newick, to_newick, bipartitions, fn_fp

rows = list(csv.DictReader(open(os.path.join(os.path.dirname(__file__), '..', 'data', 'IE_canby2024.csv'))))
names = list(rows[0])[3:]
n = len(names)
chars = [[frozenset(r[l].split('/')) for l in names] for r in rows]
w = np.array([float(r['weight']) for r in rows])
ls, lb = fitch.encode_sets(chars, n)
rng = np.random.default_rng(7)
start = methods.nj_tree(chars, n)
found = {}
for k in range(40):
    t, sc = fitch.search(ls, lb, w, fitch.Objective('mp'), [start if k == 0 else 'random'], rng)
    key = frozenset(bipartitions(t))
    found.setdefault(key, [sc, t, 0]); found[key][2] += 1
best = min(v[0] for v in found.values())
ref = from_nested(parse_newick(RWT_NWK), names)
out = ['# Weighted MP on Canby et al. 2024 IE dataset (370 chars, polymorphism kept)', '',
       f'reference-tree weighted MP score: {fitch.score_tree(ref, ls, w, lb, fitch.Objective("mp"))[0]:.0f}', '',
       '| score | #starts hitting | RF to ref | ' + ' | '.join(CLADES) + ' |', '|' + '---|' * (3 + len(CLADES))]
for key, (sc, t, cnt) in sorted(found.items(), key=lambda kv: kv[1][0]):
    if sc > best + 2:
        continue
    h = subgroup_hits(t, names)
    out.append(f'| {sc:.0f} | {cnt} | {fn_fp(ref, t)[2]} | ' + ' | '.join('Y' if h[c] else '.' for c in CLADES) + ' |')
    out.append('')
    out.append('`' + to_newick(t, names) + '`')
    out.append('')
open(os.path.join(os.path.dirname(__file__), '..', 'results', 'canby_mp.md'), 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
