"""Paired comparison of BSCAMPP runs scored by score_multi.py.
Usage: python analyze.py <scores.tsv> <ref_label> [<times.tsv>]
Prints per-run mean delta (over queries placed by every run), and vs ref: mean diff,
W/T/L (W = run has lower delta), Wilcoxon signed-rank p; and also fix vs stock at equal size."""
import sys
import pandas as pd
from scipy.stats import wilcoxon

d = pd.read_csv(sys.argv[1], sep='\t')
ref = sys.argv[2]
p = d.pivot_table(index='query', columns='label', values='delta')
n_all = len(p)
p = p.dropna()
times = {}
if len(sys.argv) > 3:
    t = pd.read_csv(sys.argv[3], sep='\t', header=None, names=['v', 'b', 'q', 'wall', 'rss'])
    for _, r in t.iterrows():
        times[f'{r.v}_b{r.b}_{r.q}'] = (r.wall, r.rss / 1e6)


def cmp(a, b):
    x, y = p[a], p[b]
    diff = x - y
    w, t_, l = (diff < 0).sum(), (diff == 0).sum(), (diff > 0).sum()
    pv = wilcoxon(x, y).pvalue if (diff != 0).any() else 1.0
    return f'{diff.mean():+.3f}', f'{w}/{t_}/{l}', f'{pv:.2g}'


print(f'queries placed by all runs: {len(p)} of {n_all}\n')
print('| run | mean delta | median | % delta=0 | vs ' + ref + ': diff | W/T/L | Wilcoxon p | wall (m:s) | max RSS (GB) |')
print('|---|---|---|---|---|---|---|---|---|')
order = sorted(p.columns, key=lambda c: (int(c.split('_b')[1].split('_')[0]), c))
for c in order:
    s = cmp(c, ref) if c != ref else ('', '', '')
    w, m = times.get(c, ('', float('nan')))
    print(f'| {c} | {p[c].mean():.3f} | {p[c].median():.0f} | {100*(p[c]==0).mean():.0f} | {s[0]} | {s[1]} | {s[2]} | {w} | {m:.1f} |')
print('\nfix vs stock at the same subtree size (W = fix lower):')
for c in order:
    if c.startswith('fix_'):
        sc = 'stock_' + c[4:]
        if sc in p:
            print(f'  {c} vs {sc}: diff/WTL/p = {cmp(c, sc)}')
for v in ('bug2only', 'bug1only'):
    for c in order:
        if c.startswith(v + '_'):
            for o in ('fix_', 'stock_'):
                oc = o + c[len(v) + 1:]
                if oc in p:
                    print(f'  {c} vs {oc}: diff/WTL/p = {cmp(c, oc)}; identical delta on {(p[c]==p[oc]).mean()*100:.1f}% of queries')
