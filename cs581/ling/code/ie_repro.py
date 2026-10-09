"""Reproduce published properties of the Ringe-Warnow-Taylor (RWT) IE datasets.

1. Compatibility of the screened / unscreened characters on an RWT-style reference tree
   (deck: '95% of the characters compatible'; all '!' characters required compatible).
2. Run the methods used in Nakhleh et al. 2005 TPS (MP, MC, WMC, NJ, plus Mk and binary ML)
   and check the subgroups the deck says all methods except UPGMA recover.

Usage: python3 ie_repro.py  (writes ../results/ie_repro.md and ie_trees.nwk)
"""
import os
import sys
import time

import numpy as np

import fitch
import methods
from trees import from_nested, parse_newick, to_newick, bipartitions, splits_from_clades, fn_fp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
RES = os.path.join(HERE, '..', 'results')

# Our transcription of the RWT 2002 / Nakhleh-Ringe-Warnow 2005 tree as drawn in the deck.
# The deck figure fixes the major-subgroup backbone; the resolution inside Anatolian, Iranian,
# Germanic and Baltic is our own guess (it does not affect the subgroup checks below).
RWT_NWK = ('((HI,(LU,LY)),((TA,TB),(((LA,(OS,UM)),(OI,WE)),((GO,(ON,(OE,OG))),'
           '(AL,((GK,AR),((VE,(AV,PE)),(OC,(PR,(LI,LT))))))))));')

CLADES = {
    'Anatolian': ['HI', 'LU', 'LY'],
    'Tocharian': ['TA', 'TB'],
    'Indo-Iranian': ['VE', 'AV', 'PE'],
    'Italic': ['LA', 'OS', 'UM'],
    'Celtic': ['OI', 'WE'],
    'Germanic': ['GO', 'ON', 'OE', 'OG'],
    'Balto-Slavic': ['OC', 'PR', 'LI', 'LT'],
    'Baltic': ['PR', 'LI', 'LT'],
    'Iranian': ['AV', 'PE'],
    'Anatolian+Tocharian': ['HI', 'LU', 'LY', 'TA', 'TB'],
    'Greco-Armenian': ['GK', 'AR'],
    'Italo-Celtic': ['LA', 'OS', 'UM', 'OI', 'WE'],
    'Satem core': ['VE', 'AV', 'PE', 'OC', 'PR', 'LI', 'LT'],
}


def load(path):
    lines = [l.split() for l in open(path) if l.strip()]
    names = lines[0]
    chars, types, ids, req, screened_out = [], [], [], [], []
    for row in lines[1:]:
        cid, nm, vals = row[0], row[1], row[2:]
        assert len(vals) == len(names), row
        chars.append([frozenset([int(v)]) for v in vals])
        types.append('P' if nm[0] == 'P' and nm[1:].isdigit() else
                     'M' if nm[0] == 'M' and nm[1].isdigit() else 'L')
        ids.append(nm)
        req.append(cid.startswith('!'))
        screened_out.append(cid.startswith('*'))
    return names, chars, types, ids, np.array(req), np.array(screened_out)


def compat_report(t, chars, types, req, n):
    ls, lb = fitch.encode_sets(chars, n)
    costs = fitch.char_costs(t, ls)
    comp = costs == lb
    out = {'all': (int(comp.sum()), len(chars))}
    for ty in 'LMP':
        m = np.array([x == ty for x in types])
        out[ty] = (int(comp[m].sum()), int(m.sum()))
    out['required(!)'] = (int(comp[req].sum()), int(req.sum()))
    out['extra_steps_total'] = int((costs - lb).sum())
    return out, comp


def subgroup_hits(t, names):
    bp = bipartitions(t)
    return {k: (m in bp) for k, m in zip(CLADES, splits_from_clades(CLADES.values(), names))}


def main():
    rng = np.random.default_rng(2026)
    out = ['# IE reproduction (RWT screened / unscreened datasets)', '']
    trees_out = []
    for label, fn in [('screened', 'IE_screened_RWT.txt'), ('unscreened', 'IE_unscreened_RWT.txt')]:
        names, chars, types, ids, req, sout = load(os.path.join(DATA, fn))
        n = len(names)
        nst = [len(set().union(*c)) for c in chars]
        out.append(f'## {label}: {len(chars)} characters '
                   f'(L={types.count("L")}, M={types.count("M")}, P={types.count("P")}), {n} languages; '
                   f'{int(req.sum())} marked "!" (required compatible), {int(sout.sum())} marked "*" (removed by screening)')
        for ty in 'LMP':
            v = [k for k, t in zip(nst, types) if t == ty]
            out.append(f'- type {ty}: mean #states {np.mean(v):.2f}, median {np.median(v):.0f}, max {max(v)}')
        ref = from_nested(parse_newick(RWT_NWK), names)
        rep, comp = compat_report(ref, chars, types, req, n)
        out.append(f'- **Reference tree** compatibility: {rep}')
        out.append('  incompatible on reference tree: ' + ', '.join(i for i, c in zip(ids, comp) if not c))
        rows = []
        W = {'L': 1.0, 'M': 1.0, 'P': 1.0}
        Wreq = np.where(req, 100.0, 1.0)  # RWT: maximise compatible chars s.t. '!' chars compatible
        start = methods.nj_tree(chars, n)
        runs = {
            'MP': lambda: methods.char_search(chars, types, n, 'mp', rng=rng, nstarts=10, start_trees=[start])[0],
            'MC': lambda: methods.char_search(chars, types, n, 'mc', rng=rng, nstarts=10, start_trees=[start])[0],
            'WMC(!x100)': lambda: _wsearch(chars, n, Wreq, 'mc', rng, start),
            'WMC(M,P x5)': lambda: methods.char_search(chars, types, n, 'mc', weights={'L': 1, 'M': 5, 'P': 5},
                                                       rng=rng, nstarts=10, start_trees=[start])[0],
            'CapPars(2)': lambda: methods.char_search(chars, types, n, 'cap', cap=2, rng=rng, nstarts=10,
                                                      start_trees=[start])[0],
            'NJ': lambda: methods.nj_tree(chars, n),
            'ML-Mk': lambda: methods.iqtree_mk(chars, n, seed=1),
            'ML-binary': lambda: methods.iqtree_binary(chars, n, seed=1),
        }
        out.append('')
        out.append('| method | secs | compat (all/L/M/P/!) | RF to ref | ' + ' | '.join(CLADES) + ' |')
        out.append('|' + '---|' * (4 + len(CLADES)))
        for m, f in [('reference', lambda: ref)] + list(runs.items()):
            t0 = time.time(); t = f(); dt = time.time() - t0
            rep, _ = compat_report(t, chars, types, req, n)
            hits = subgroup_hits(t, names)
            rf = fn_fp(ref, t)[2]
            out.append(f'| {m} | {dt:.1f} | {rep["all"][0]}/{rep["L"][0]}/{rep["M"][0]}/{rep["P"][0]}/{rep["required(!)"][0]} '
                       f'| {rf} | ' + ' | '.join('Y' if hits[k] else '.' for k in CLADES) + ' |')
            trees_out.append(f'[{label} {m}] ' + to_newick(t, names))
            print(label, m, round(dt, 1), rep['all'], file=sys.stderr)
        out.append('')
    open(os.path.join(RES, 'ie_repro.md'), 'w').write('\n'.join(out) + '\n')
    open(os.path.join(RES, 'ie_trees.nwk'), 'w').write('\n'.join(trees_out) + '\n')


def _wsearch(chars, n, w, kind, rng, start):
    ls, lb = fitch.encode_sets(chars, n)
    t, _ = fitch.search(ls, lb, w, fitch.Objective(kind), [start] + ['random'] * 10, rng)
    return t


if __name__ == '__main__':
    main()
