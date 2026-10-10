"""Markdown tables for REPORT.md from the pilot outputs.
Usage: python make_tables.py <results_dir>
Reads <rep>_nested_scores.tsv, <rep>_bscampp_scores.tsv, <rep>_bscampp_times.tsv, <rep>_apples_scores.tsv"""
import sys, glob, os
import pandas as pd
from scipy.stats import wilcoxon

R = sys.argv[1]


def wtl(x, y):
    d = x - y
    w, t, l = (d < 0).sum(), (d == 0).sum(), (d > 0).sum()
    p = wilcoxon(x, y).pvalue if (d != 0).any() else 1.0
    return f'{d.mean():+.3f} | {w}/{t}/{l} | {p:.1g}'


def nested(rep):
    f = f'{R}/{rep}_nested_scores.tsv'
    if not os.path.exists(f):
        return
    d = pd.read_csv(f, sep='\t')
    for qt in ['frag', 'full']:
        x = d[d.qtype == qt]
        p = x.pivot_table(index='query', columns=['mode', 'k'], values='delta')
        modes = [m for m in ['auto', 'fix', 'nopremask', 'noheur', 'rs_on', 'diagcpu', 'pplacer']
                 if m in p.columns.get_level_values(0)]
        ks = sorted(set(p.columns.get_level_values(1)))
        print(f'\n**{rep}, {qt} queries: mean delta error by subtree size** '
              f'(n = queries with all cells present for that mode)\n')
        print('| mode | n | ' + ' | '.join(str(k) for k in ks) + ' |')
        print('|---|---|' + '---|' * len(ks))
        for m in modes:
            cols = [(m, k) for k in ks if (m, k) in p]
            sub = p[cols].dropna()
            vals = [f'{sub[(m, k)].mean():.3f}' if (m, k) in sub else '-' for k in ks]
            print(f'| {m} | {len(sub)} | ' + ' | '.join(vals) + ' |')
        print(f'\nPaired vs stock EPA-ng (auto) at the same size, {qt}: mean diff | W/T/L (W = lower error) | Wilcoxon p\n')
        print('| mode | k | diff | W/T/L | p |')
        print('|---|---|---|---|---|')
        for m in modes:
            if m == 'auto':
                continue
            for k in ks:
                if (m, k) in p and ('auto', k) in p:
                    s = p[[(m, k), ('auto', k)]].dropna()
                    if len(s) == 0:
                        continue
                    print(f'| {m} | {k} | ' + wtl(s[(m, k)], s[('auto', k)]) + ' |')


def bscampp(rep):
    f = f'{R}/{rep}_bscampp_scores.tsv'
    if not os.path.exists(f):
        return
    d = pd.read_csv(f, sep='\t', names=['label', 'query', 'delta', 'lwr'])
    d['variant'] = d.label.str.split('_').str[0]
    d['b'] = d.label.str.split('_').str[1].str[1:].astype(int)
    d['qtype'] = d.label.str.split('_').str[2]
    t = pd.read_csv(f'{R}/{rep}_bscampp_times.tsv', sep='\t', names=['variant', 'b', 'qtype', 'wall', 'rss'])
    t = t.drop_duplicates(['variant', 'b', 'qtype'], keep='last')

    def secs(s):
        parts = [float(x) for x in str(s).split(':')]
        v = 0
        for x in parts:
            v = v * 60 + x
        return v
    t['sec'] = t.wall.map(secs)
    a = None
    f2 = f'{R}/{rep}_apples_scores.tsv'
    if os.path.exists(f2):
        a = pd.read_csv(f2, sep='\t', names=['label', 'query', 'delta', 'lwr'])
    for qt in ['frag', 'full']:
        x = d[d.qtype == qt]
        p = x.pivot_table(index='query', columns=['variant', 'b'], values='delta')
        if ('stock', 2000) not in p:
            continue
        base = p[('stock', 2000)]
        print(f'\n**{rep}, BSCAMPP(e), {qt} queries (n={len(base)})**: baseline = stock EPA-ng, b=2000\n')
        print('| EPA-ng | subtree b | mean delta | diff vs baseline | W/T/L | p | wall (s) | peak RSS (GB) |')
        print('|---|---|---|---|---|---|---|---|')
        for (v, b) in sorted(p.columns, key=lambda c: (c[0] != 'stock', c[1])):
            s = p[[(v, b)]].join(base, rsuffix='_b').dropna()
            tt = t[(t.variant == v) & (t.b == b) & (t.qtype == qt)]
            ws = f'{tt.sec.iloc[0]:.0f}' if len(tt) else '-'
            rs = f'{tt.rss.iloc[0] / 1e6:.1f}' if len(tt) else '-'
            cmp = wtl(s.iloc[:, 0], s.iloc[:, 1]) if (v, b) != ('stock', 2000) else '- | - | -'
            print(f'| {v} | {b} | {p[(v, b)].mean():.3f} | {cmp} | {ws} | {rs} |')
        if a is not None:
            aa = a[a.label == f'apples_{qt}']
            if len(aa):
                print(f'| APPLES-2 (full tree) | - | {aa.delta.mean():.3f} | | | | | |')


for rep in sorted({os.path.basename(f).split('_')[0] for f in glob.glob(f'{R}/R*_*.tsv')}):
    nested(rep)
    bscampp(rep)
