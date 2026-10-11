"""Paired summaries of nested scores: mean delta by mode x k over queries present in all
compared cells, plus paired W/T/L and Wilcoxon signed-rank vs a reference mode.
Usage: python summarize.py scores.tsv qtype ref_mode mode1,mode2,... [k1,k2,...]"""
import sys
import pandas as pd
from scipy.stats import wilcoxon

d = pd.read_csv(sys.argv[1], sep='\t')
qt, ref = sys.argv[2], sys.argv[3]
modes = sys.argv[4].split(',')
d = d[(d.qtype == qt) & d['mode'].isin(modes + [ref])]
if len(sys.argv) > 5:
    d = d[d.k.isin([int(x) for x in sys.argv[5].split(',')])]
p = d.pivot_table(index='query', columns=['mode', 'k'], values='delta')
p = p.dropna()
print(f'qtype={qt}  paired queries={len(p)}')
print(p.mean().unstack('k').round(3).to_string())
print(f'\nvs {ref} at the same k: mean diff (mode - ref), W/T/L (W = mode lower), Wilcoxon p')
for m in modes:
    if m == ref:
        continue
    for k in sorted(set(c[1] for c in p.columns)):
        if (m, k) not in p or (ref, k) not in p:
            continue
        x, y = p[(m, k)], p[(ref, k)]
        diff = x - y
        w, t, l = (diff < 0).sum(), (diff == 0).sum(), (diff > 0).sum()
        pv = wilcoxon(x, y).pvalue if (diff != 0).any() else 1.0
        print(f'{m:>12} k={k:<5} diff={diff.mean():+.3f}  W/T/L={w}/{t}/{l}  p={pv:.2g}')
