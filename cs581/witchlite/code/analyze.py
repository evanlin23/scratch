#!/usr/bin/env python3
"""Tables from results/*.json: pooled SPFN/SPFP, paired per-query difference vs WITCH,
W/T/L (per-query SPFN lower/equal/higher than WITCH), Wilcoxon signed-rank p, runtime vs WITCH.
usage: analyze.py RESULTS_DIR > table.md
"""
import json, sys, glob, os
import numpy as np
from scipy.stats import wilcoxon

ORDER = ['WITCH', 'all/k10', 'all/k3', 'all/k1', 'all/k10_t99', 'all/k10_t95', 'hier/k10', 'hier/k10_t99',
         'hier_es/k10', 'beam/k10', 'beam/k10_t99', 'blastpath/k10', 'blastpath/k10_t99', 'blastpath_sib/k10',
         'blastpath_sib/k10_t99', 'blastpath_sib/k10_t95', 'BLAST']
LABEL = {'WITCH': 'WITCH (default, own run)', 'all/k10': 'all HMMs, k=10 (WITCH re-run via weights)',
         'all/k3': 'all HMMs, k=3', 'all/k1': 'all HMMs, k=1 (UPP-like, adj. bit-score)',
         'all/k10_t99': 'all HMMs, adaptive k (tau=.99)', 'all/k10_t95': 'all HMMs, adaptive k (tau=.95)',
         'hier/k10': 'hier descent (UPP2-style), k=10', 'hier/k10_t99': 'hier descent, adaptive k',
         'hier_es/k10': 'hier EarlyStop (UPP2), k<=10', 'beam/k10': 'beam-2 descent, k=10',
         'beam/k10_t99': 'beam-2 descent, adaptive k', 'blastpath/k10': 'BLAST path, k<=10',
         'blastpath/k10_t99': 'BLAST path, adaptive k', 'blastpath_sib/k10': 'BLAST path+siblings, k=10',
         'blastpath_sib/k10_t99': 'BLAST path+siblings, adaptive k (tau=.99)',
         'blastpath_sib/k10_t95': 'BLAST path+siblings, adaptive k (tau=.95)', 'BLAST': 'BLASTN only (TIPP3-fast)'}

for f in sorted(glob.glob(f'{sys.argv[1]}/*.json')):
    r = json.load(open(f))
    if 'methods' not in r or 'WITCH' not in r['methods']:
        continue
    M = r['methods']
    wt = r['witch_default_time']
    base = np.array(M['WITCH']['spfn_q'])
    print(f"\n### {os.path.basename(f)[:-5]}\n")
    print(f"WITCH default wall {wt['wall']:.0f} s = decomposition {wt['decomposition']:.0f} s + "
          f"all-vs-all hmmsearch {wt['search']:.0f} s + alignment/merge/other {wt['rest']:.0f} s "
          f"(search = {100 * wt['search'] / wt['wall']:.0f}% of wall)\n")
    print('| method | SPFN % | SPFP % | dSPFN vs WITCH (pp) | W/T/L vs WITCH | Wilcoxon p | HMMs scored/query | mean k | '
          'time (s) | % of WITCH time | % excl. decomposition |')
    print('|---|---|---|---|---|---|---|---|---|---|---|')
    for m in ORDER:
        if m not in M:
            continue
        d = M[m]
        q = np.array(d['spfn_q'])
        diff = q - base
        w, t, l = int((diff < -1e-12).sum()), int((abs(diff) <= 1e-12).sum()), int((diff > 1e-12).sum())
        p = wilcoxon(q, base).pvalue if (w + l) > 0 else 1.0
        tm = d['time']; tn = d['time_nodecomp']
        sp = d.get('scores_per_query', '')
        sp = f'{sp:.1f}' if sp != '' else ('all' if m == 'WITCH' else '-')
        mk = d.get('mean_k', '')
        mk = f'{mk:.1f}' if mk != '' else ('10' if m == 'WITCH' else '-')
        print(f"| {LABEL[m]} | {100 * d['spfn']:.2f} | {100 * d['spfp']:.2f} | {100 * (d['spfn'] - M['WITCH']['spfn']):+.2f} | "
              f"{w}/{t}/{l} | {p:.2g} | {sp} | {mk} | {tm:.0f} | {100 * tm / wt['wall']:.0f}% | "
              f"{100 * tn / (wt['wall'] - wt['decomposition']):.0f}% |")
