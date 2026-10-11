"""RAxML-NG 2.0.3 --fast (KH-mult) vs classic single-start search, paired sequential runs: speed-up and lnL."""
import glob, os, re, numpy as np
from scipy.stats import wilcoxon
W = '/opt/work/runs'; rows = []
for d in sorted(glob.glob(f'{W}/*')):
    a, b = f'{d}/rxfastB.done', f'{d}/rxclassic1.done'
    if not (os.path.exists(a) and os.path.exists(b)): continue
    ta, tb = float(open(a).read().split()[3]), float(open(b).read().split()[3])
    la = float(re.search(r'Final LogLikelihood: (\S+)', open(f'{d}/rxfastB.raxml.log').read()).group(1))
    lb = float(re.search(r'Final LogLikelihood: (\S+)', open(f'{d}/rxclassic1.raxml.log').read()).group(1))
    rows.append((os.path.basename(d), ta, tb, tb / ta, la - lb))
out = ['| dataset | --fast (s) | classic pars{1} (s) | speed-up | lnL(fast) - lnL(classic) |', '|---|---|---|---|---|']
out += [f'| {r[0]} | {r[1]:.1f} | {r[2]:.1f} | {r[3]:.2f} | {r[4]:.2f} |' for r in rows]
if rows:
    s = np.array([r[3] for r in rows]); dl = np.array([r[4] for r in rows])
    emp = [r for r in rows if r[0].startswith('emp')]
    out.append(f'\nAll n={len(rows)}: geo-mean speed-up {np.exp(np.log(s).mean()):.2f}x (DNA 16S only: '
               f'{np.exp(np.mean([np.log(r[3]) for r in emp])):.2f}x, n={len(emp)}); mean dlnL {dl.mean():.2f}, '
               f'W/T/L (tie 0.5) {(dl>0.5).sum()}/{(abs(dl)<=0.5).sum()}/{(dl<-0.5).sum()}, Wilcoxon p={wilcoxon(dl).pvalue:.3g}')
open('/home/user/scratch/cs581/iqstop/results/raxml_repro.md', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
