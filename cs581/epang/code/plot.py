"""Figure: mean delta error vs subtree size (nested experiment, fragmentary queries).
Usage: python plot.py <nested_scores.tsv> <out.png>"""
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker

d = pd.read_csv(sys.argv[1], sep='\t')
series = [('auto', 'EPA-ng 0.3.8 (stock)', '#2a78d6', 'o'),
          ('fix', 'EPA-ng patched', '#eb6834', 's'),
          ('nopremask', 'EPA-ng --no-pre-mask', '#1baf7a', '^'),
          ('pplacer', 'pplacer', '#eda100', 'D')]
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=False)
for ax, qt in zip(axes, ['frag', 'full']):
    x = d[d.qtype == qt]
    p = x.pivot_table(index='query', columns=['mode', 'k'], values='delta')
    for m, lab, col, mk in series:
        cols = [c for c in p.columns if c[0] == m]
        if not cols:
            continue
        s = p[cols].dropna().mean()
        ks = [c[1] for c in s.index]
        ax.plot(ks, s.values, color=col, marker=mk, ms=7, lw=2, label=lab)
    ax.axvline(2000, color='#999999', lw=1, ls=':')
    ax.text(2050, ax.get_ylim()[1] * 0.97 if ax.get_ylim()[1] else 1, '2,000 tips:\nrate scalers on',
            fontsize=8, color='#555555', va='top')
    ax.set_xscale('log')
    ax.set_xticks([500, 1000, 2000, 3000, 5000, 9000])
    ax.set_xticklabels(['500', '1k', '2k', '3k', '5k', '9k'])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel('placement tree size (leaves)')
    ax.set_ylabel('mean delta error (edges)')
    ax.set_title('fragmentary queries (~154 nt)' if qt == 'frag' else 'full-length queries (stock = patched, lines overlap)')
    ax.grid(axis='y', color='#e5e5e5', lw=0.8)
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
axes[0].legend(frameon=False, fontsize=8)
fig.suptitle('Same queries placed into nested subtrees of growing size (RNASim, 9,000-leaf backbone)', fontsize=10)
fig.tight_layout()
fig.savefig(sys.argv[2], dpi=130)
