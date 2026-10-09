"""Summarise ../results/sim_runs.jsonl -> ../results/sim_summary.md

* Mean tree FN rate (binary trees, so FN = FP = RF/2(n-3)) per method x condition, test reps.
* Capped parsimony: cap chosen on TRAINING reps (lowest mean FN pooled over conditions,
  unweighted and weighted variants separately), then evaluated on TEST reps.
* Paired two-sided Wilcoxon signed-rank tests (scipy, zero_method='wilcox'): pairs whose FN
  differs by less than the tie band (1e-9, i.e. exact ties; FN moves in steps of 1/(n-3)=0.048)
  are dropped; W/T/L counts are reported alongside.
"""
import json
import os
from collections import defaultdict

import numpy as np
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
TIE = 1e-9


def load():
    return [json.loads(l) for l in open(os.path.join(RES, 'sim_runs.jsonl'))]


def paired(rows, a, b):
    d = np.array([r['fn'][a] - r['fn'][b] for r in rows])
    w, t, l = int((d < -TIE).sum()), int((np.abs(d) <= TIE).sum()), int((d > TIE).sum())
    nz = d[np.abs(d) > TIE]
    p = wilcoxon(nz).pvalue if len(nz) >= 5 else float('nan')
    return d.mean(), w, t, l, p


def main():
    rows = load()
    conds = list(dict.fromkeys(r['cond'] for r in rows))
    methods = list(rows[0]['fn'])
    train = [r for r in rows if r['split'] == 'train']
    test = [r for r in rows if r['split'] == 'test']
    out = ['# Simulation results', '',
           f'{len(rows)} replicate runs ({len(train)} train, {len(test)} test); n = 24 languages.', '']

    # cap selection on training
    caps = sorted({m for m in methods if m.startswith('Cap')}, key=lambda s: int(s[3:]))
    wcaps = sorted({m for m in methods if m.startswith('WCap')}, key=lambda s: int(s[4:]))
    cand_u = ['MC'] + caps + ['MP']
    cand_w = ['WMC'] + wcaps
    tm = {m: np.mean([r['fn'][m] for r in train]) for m in cand_u + cand_w}
    best_u = min(cand_u, key=lambda m: tm[m]); best_w = min(cand_w, key=lambda m: tm[m])
    out.append('## Cap selection on training replicates (mean FN pooled over conditions)')
    out.append('')
    out.append('| ' + ' | '.join(cand_u + cand_w) + ' |')
    out.append('|' + '---|' * len(cand_u + cand_w))
    out.append('| ' + ' | '.join(f'{tm[m]:.4f}' for m in cand_u + cand_w) + ' |')
    out.append('')
    out.append(f'Selected: unweighted **{best_u}**, weighted **{best_w}** (pilot "new" method = selected cap).')
    out.append('')

    # FN table
    out.append('## Mean FN rate on TEST replicates (lower is better)')
    out.append('')
    out.append('| method | ' + ' | '.join(conds) + ' | all |')
    out.append('|' + '---|' * (len(conds) + 2))
    for m in methods:
        vals = []
        for c in conds:
            v = [r['fn'][m] for r in test if r['cond'] == c]
            vals.append(f'{np.mean(v):.3f}' if v else '-')
        vals.append(f'{np.mean([r["fn"][m] for r in test]):.3f}')
        out.append(f'| {m} | ' + ' | '.join(vals) + ' |')
    out.append('| (reps) | ' + ' | '.join(str(sum(r['cond'] == c for r in test)) for c in conds) + f' | {len(test)} |')
    out.append('')

    # data summaries
    out.append('## Simulated data summaries (all reps; compare with IE: lexical ~14.3 states/char, morph ~9-10, phon 2.3)')
    out.append('')
    keys = list(rows[0]['summary'])
    out.append('| cond | ' + ' | '.join(keys) + ' | borrow events |')
    out.append('|' + '---|' * (len(keys) + 2))
    for c in conds:
        rr = [r for r in rows if r['cond'] == c]
        out.append(f'| {c} | ' + ' | '.join(f'{np.mean([r["summary"][k] for r in rr]):.2f}' for k in keys)
                   + f' | {np.mean([r["borrow"] for r in rr]):.1f} |')
    out.append('')

    # paired tests
    comps = [(best_u, 'MP'), (best_u, 'MC'), (best_u, 'ML-Mk'), (best_u, 'ML-bin'), (best_u, 'NJ'),
             (best_w, 'WMC'), (best_w, best_u), ('MP-poly', 'MP'), ('ML-Mk', 'ML-bin'), ('MP', 'ML-Mk'),
             ('WMC', 'MC'), ('MP', 'NJ')]
    out.append('## Paired Wilcoxon signed-rank tests on TEST reps (A vs B; mean dFN = A - B; W/T/L = A better/tie/worse)')
    out.append('')
    out.append('| A | B | scope | mean dFN | W/T/L | p |')
    out.append('|---|---|---|---|---|---|')
    for a, b in comps:
        for scope in ['all'] + conds:
            rr = test if scope == 'all' else [r for r in test if r['cond'] == scope]
            if not rr:
                continue
            md, w, t, l, p = paired(rr, a, b)
            out.append(f'| {a} | {b} | {scope} | {md:+.4f} | {w}/{t}/{l} | {p:.3g} |')
    out.append('')

    # runtime
    out.append('## Mean runtime per replicate (seconds, one core)')
    out.append('')
    out.append('| ' + ' | '.join(methods) + ' |')
    out.append('|' + '---|' * len(methods))
    out.append('| ' + ' | '.join(f'{np.mean([r["secs"][m] for r in rows]):.2f}' for m in methods) + ' |')
    open(os.path.join(RES, 'sim_summary.md'), 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))


if __name__ == '__main__':
    main()
