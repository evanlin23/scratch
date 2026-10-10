"""(1) KH significance of each successive new-best topology in the default IQ-TREE runs;
(2) speed-up vs lnL-loss frontier plot of nstop-N vs KH rules (per data group)."""
import glob, json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
import replay
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R = '/home/user/scratch/cs581/iqstop/results'
ev = []
for f in sorted(glob.glob(f'{R}/per_dataset/*.json')):
    r = json.load(open(f)); ds = r['ds']
    grp = 'sim' if ds.startswith('rg') else ('emp-AA' if ds.startswith('empAA') else 'emp-DNA')
    names = replay.taxa_of(replay.aln_of(ds))
    rows = [(a, b, c, d, replay.relabel(e, names) if e else None) for a, b, c, d, e in replay.read_trace(f'{replay.W}/{ds}/iqdef1.itertrace')]
    trees = r['trees']; L, _ = replay.evaluate(ds, trees)
    prev = None
    for it, t, sc, best, tr in rows:
        if tr is None: continue
        k = trees.index(tr)
        if prev is not None:
            sig, delta, sd = replay.kh_better(L[:, k], L[:, prev], 0.05)
            ev.append(dict(ds=ds, group=grp, it=it, phase='init(<=20)' if it <= 20 else ('21-100' if it <= 100 else '>100'),
                           delta=delta, sd=sd, sig=bool(sig), t_frac=t / r['wall_default']))
        prev = k
E = pd.DataFrame(ev); E.to_csv(f'{R}/improvements.csv', index=False)
tab = E.groupby(['group', 'phase']).agg(n=('sig', 'size'), frac_KH_sig=('sig', 'mean'), median_gain=('delta', 'median'),
                                       mean_gain=('delta', 'mean')).reset_index()
md = ['# New-best topologies in default IQ-TREE runs: gain vs previous best and one-sided KH (alpha 0.05)\n',
      'Gains are re-evaluated lnL differences (common model fit, branch lengths re-optimized).\n', tab.to_markdown(index=False, floatfmt='.3g')]
open(f'{R}/improvements.md', 'w').write('\n'.join(md) + '\n'); print('\n'.join(md))
# frontier plot
df = pd.read_csv(f'{R}/rows.csv')
fig, axs = plt.subplots(1, 3, figsize=(13, 4))
fam = {'nstop': ['default(nstop100)', 'nstop50', 'nstop20', 'nstop10'], 'KHpat': ['KHpat100', 'KHpat50', 'KHpat20'],
       'KHwin': ['KHwin20', 'KHwin10']}
cols = {'nstop': '#555555', 'KHpat': '#1f77b4', 'KHwin': '#d62728'}
for ax, g in zip(axs, ['sim', 'emp-DNA', 'emp-AA']):
    sub = df[df.group == g]
    if sub.empty: ax.set_title(f'{g}: no data'); continue
    for fm, ms in fam.items():
        xs = [np.exp(np.log(sub[sub.method == m].speedup).mean()) for m in ms]
        ys = [sub[sub.method == m].dlnl.mean() for m in ms]
        ax.plot(xs, ys, 'o-', color=cols[fm], label=fm)
        for m, x, y in zip(ms, xs, ys): ax.annotate(m.replace('default(nstop100)', 'default'), (x, y), fontsize=7)
    if 'iqdef2' in set(sub.method):
        n2 = sub[sub.method == 'iqdef2'].dlnl
        ax.axhspan(-n2.abs().mean(), 0, color='#cccccc', alpha=.4, label='|seed2 - seed1| mean')
    ax.set_xscale('log'); ax.set_xlabel('geo-mean speed-up vs IQ-TREE default'); ax.set_ylabel('mean dlnL vs default tree')
    ax.set_title(f'{g} (n={sub.ds.nunique()})'); ax.legend(fontsize=7)
plt.tight_layout(); plt.savefig(f'{R}/frontier.png', dpi=120)
