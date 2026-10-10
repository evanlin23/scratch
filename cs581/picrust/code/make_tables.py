#!/usr/bin/env python3
"""Collect results/*_placement.json and *_compare.json into markdown tables (stdout)."""
import json, os, glob
R = '/home/user/scratch/cs581/picrust/results'
def L(p):
    return json.load(open(p)) if os.path.exists(p) else None
rows = []
print('| dataset | ASVs | aln cols | first col (median) | logL changed | edge changed | median / max edge dist | LWR≥0.99 stock→fix | mean LWR of moved, stock→fix |')
print('|---|---|---|---|---|---|---|---|---|')
for ds in ['ocean', 'mammal', 'indian', 'hmp', 'loo_v4', 'loo_v3_130']:
    p = L(f'{R}/{ds}_placement.json')
    if not p: continue
    n = p['n']
    print(f"| {ds} | {n} | {p.get('aln_cols','')} | {p.get('first_col_median','')} | {p['logl_changed']} ({100*p['logl_changed']/n:.1f}%) | "
          f"{p['edge_changed']} ({100*p['edge_changed']/n:.1f}%) | {p['edge_dist_median']} / {p['edge_dist_max']} | "
          f"{p['lwr99_stock']}→{p['lwr99_fix']} | " + (f"{p['lwr_changed_stock']:.2f}→{p['lwr_changed_fix']:.2f}" if p['lwr_changed_stock'] is not None else '–') + ' |')
print()
print('| dataset | ASV NSTI changed (up/down) | mean ASV NSTI stock→fix | closest ref changed | ASV KO vectors changed | KO presence changes / sample (median, max) | per-sample Spearman stock vs fix KO (median, min) | pathway (median, min) | KO rel-abun shift ½L1 (median, max) |')
print('|---|---|---|---|---|---|---|---|---|')
for ds in ['ocean', 'mammal', 'indian', 'hmp']:
    c = L(f'{R}/{ds}_compare.json')
    if not c: continue
    print(f"| {ds} | {c['asv_nsti_changed']} ({c['asv_nsti_fix_higher']}/{c['asv_nsti_fix_lower']}) | {c['asv_nsti_mean_stock']:.3f}→{c['asv_nsti_mean_fix']:.3f} | "
          f"{c['asv_closest_ref_changed']} | {c['asv_ko_vector_changed']} | {c['KO_presence_changes_per_sample_median']:.0f}, {c['KO_presence_changes_per_sample_max']} | "
          f"{c['KO_stock_vs_fix_spearman_median']:.3f}, {c['KO_stock_vs_fix_spearman_min']:.3f} | {c['PATH_stock_vs_fix_spearman_median']:.3f}, {c['PATH_stock_vs_fix_spearman_min']:.3f} | "
          f"{c['KO_relabun_L1half_median']:.3f}, {c['KO_relabun_L1half_max']:.3f} |")
print()
print('| dataset | samples | KO Spearman vs MGS stock→fix | fix better/worse | Wilcoxon p | KO precision stock→fix | KO recall stock→fix | pathway Spearman stock→fix | better/worse | p |')
print('|---|---|---|---|---|---|---|---|---|---|')
for ds in ['ocean', 'mammal', 'indian', 'hmp']:
    c = L(f'{R}/{ds}_compare.json')
    if not c: continue
    k = c['acc_KO_spearman_fix_better_worse_tie']; q = c['acc_PATH_spearman_fix_better_worse_tie']
    print(f"| {ds} | {c['acc_KO_n_samples']} | {c['acc_KO_spearman_stock_mean']:.4f}→{c['acc_KO_spearman_fix_mean']:.4f} | {k[0]}/{k[1]} | {c['acc_KO_wilcoxon_p']:.2g} | "
          f"{c['acc_KO_precision_stock_fix'][0]:.3f}→{c['acc_KO_precision_stock_fix'][1]:.3f} | {c['acc_KO_recall_stock_fix'][0]:.3f}→{c['acc_KO_recall_stock_fix'][1]:.3f} | "
          f"{c['acc_PATH_spearman_stock_mean']:.4f}→{c['acc_PATH_spearman_fix_mean']:.4f} | {q[0]}/{q[1]} | {c['acc_PATH_wilcoxon_p']:.2g} |")
