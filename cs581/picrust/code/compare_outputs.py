#!/usr/bin/env python3
"""Compare PICRUSt2 outputs from stock vs patched EPA-ng for one dataset, and score both against the
paired shotgun metagenomes (HUMAnN2 tables from the PICRUSt2 manuscript repo).
Usage: compare_outputs.py <dataset>   (expects /opt/work/<ds>/pipe_stock and pipe_fix)
Writes a JSON summary to stdout and per-sample tables to results/."""
import sys, json, gzip, numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon

ds = sys.argv[1]
W = f'/opt/work/{ds}'
MS = '/home/user/picrust/ms/data'
RES = '/home/user/scratch/cs581/picrust/results'
A, B = f'{W}/pipe_stock', f'{W}/pipe_fix'
out = {'dataset': ds}

def rd(p):
    return pd.read_csv(p, sep='\t', index_col=0, comment=None)

# --- per-ASV NSTI and domain ---
na, nb = rd(f'{A}/combined_marker_predicted_and_nsti.tsv.gz'), rd(f'{B}/combined_marker_predicted_and_nsti.tsv.gz')
common = na.index.intersection(nb.index)
out['asv_n'] = int(len(common))
d = (nb.loc[common, 'metadata_NSTI'] - na.loc[common, 'metadata_NSTI'])
out['asv_nsti_changed'] = int((d.abs() > 1e-9).sum())
out['asv_nsti_mean_stock'] = float(na.loc[common, 'metadata_NSTI'].mean())
out['asv_nsti_mean_fix'] = float(nb.loc[common, 'metadata_NSTI'].mean())
out['asv_nsti_absdiff_mean_changed'] = float(d[d.abs() > 1e-9].abs().mean()) if out['asv_nsti_changed'] else 0.0
out['asv_nsti_fix_lower'] = int((d < -1e-9).sum()); out['asv_nsti_fix_higher'] = int((d > 1e-9).sum())
out['asv_closest_ref_changed'] = int((na.loc[common, 'closest_reference_genome'] != nb.loc[common, 'closest_reference_genome']).sum())
out['asv_domain_changed'] = int((na.loc[common, 'best_domain'] != nb.loc[common, 'best_domain']).sum())
for t in (2.0,):
    out[f'asv_over_nsti{t}_stock'] = int((na['metadata_NSTI'] > t).sum())
    out[f'asv_over_nsti{t}_fix'] = int((nb['metadata_NSTI'] > t).sum())

# --- per-ASV KO predictions (bacteria-domain predictions, as used downstream) ---
ka, kb = rd(f'{A}/bac_KO_predicted.tsv.gz'), rd(f'{B}/bac_KO_predicted.tsv.gz')
ci = ka.index.intersection(kb.index); cc = ka.columns.intersection(kb.columns)
ka, kb = ka.loc[ci, cc], kb.loc[ci, cc]
diffrow = (ka.values != kb.values).any(axis=1)
out['asv_ko_vector_changed'] = int(diffrow.sum())
pa, pb = ka.values > 0, kb.values > 0
flips = (pa != pb).sum(axis=1)
out['asv_ko_presence_flips_mean_over_changed'] = float(flips[diffrow].mean()) if diffrow.any() else 0.0
out['asv_ko_presence_flips_total'] = int(flips.sum())

# --- sample-level KO / EC / pathway ---
def sample_cmp(tag, fa, fb):
    x, y = rd(fa), rd(fb)
    idx = x.index.union(y.index)
    x, y = x.reindex(idx).fillna(0), y.reindex(idx).fillna(0)
    y = y[x.columns]
    rho = [spearmanr(x[s], y[s])[0] for s in x.columns]
    pres = [int(((x[s] > 0) != (y[s] > 0)).sum()) for s in x.columns]
    rel = (x / x.sum()).fillna(0); rely = (y / y.sum()).fillna(0)
    bc = [0.5 * np.abs(rel[s] - rely[s]).sum() for s in x.columns]
    out[f'{tag}_n_features'] = int(len(idx))
    out[f'{tag}_stock_vs_fix_spearman_min'] = float(np.min(rho))
    out[f'{tag}_stock_vs_fix_spearman_median'] = float(np.median(rho))
    out[f'{tag}_presence_changes_per_sample_median'] = float(np.median(pres))
    out[f'{tag}_presence_changes_per_sample_max'] = int(np.max(pres))
    out[f'{tag}_relabun_L1half_median'] = float(np.median(bc))
    out[f'{tag}_relabun_L1half_max'] = float(np.max(bc))
    return x, y

ko_a, ko_b = sample_cmp('KO', f'{A}/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz', f'{B}/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz')
sample_cmp('EC', f'{A}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz', f'{B}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz')
pw_a, pw_b = sample_cmp('PATH', f'{A}/pathways_out/path_abun_unstrat.tsv.gz', f'{B}/pathways_out/path_abun_unstrat.tsv.gz')
wa, wb = rd(f'{A}/KO_metagenome_out/weighted_nsti.tsv.gz'), rd(f'{B}/KO_metagenome_out/weighted_nsti.tsv.gz')
out['weighted_nsti_median_stock'] = float(wa.iloc[:, 0].median()); out['weighted_nsti_median_fix'] = float(wb.iloc[:, 0].median())
out['weighted_nsti_maxabsdiff'] = float((wa.iloc[:, 0] - wb.loc[wa.index].iloc[:, 0]).abs().max())

# --- accuracy vs metagenomes (manuscript protocol, simplified): per-sample Spearman over the KOs
# that both PICRUSt2 (current ko.txt.gz) and HUMAnN2 can call; missing = 0 ---
def accuracy(tag, pa_, pb_, mgs_file, possible):
    strip = lambda df: df.rename(index=lambda i: i[3:] if i.startswith('ko:') else i)
    pa_, pb_ = strip(pa_), strip(pb_)
    m = rd(mgs_file)
    m = m.loc[[i for i in m.index if i not in ('UNMAPPED', 'UNGROUPED', 'UNINTEGRATED')]]
    m = m.loc[[i for i in m.index if '|' not in i]]
    feats = sorted(set(possible) & set(m.index)) if possible is not None else sorted(set(m.index) | set(pa_.index))
    samples = [s for s in pa_.columns if s in m.columns]
    M = m.reindex(feats).fillna(0)[samples]
    X = pa_.reindex(feats).fillna(0)[samples]; Y = pb_.reindex(feats).fillna(0)[samples]
    ra = np.array([spearmanr(X[s], M[s])[0] for s in samples]); rb = np.array([spearmanr(Y[s], M[s])[0] for s in samples])
    # presence/absence precision & recall per sample
    def pr(P):
        p, r = [], []
        for s in samples:
            pp, mm = P[s] > 0, M[s] > 0
            tp = (pp & mm).sum(); p.append(tp / max(pp.sum(), 1)); r.append(tp / max(mm.sum(), 1))
        return np.array(p), np.array(r)
    pa1, ra1 = pr(X); pb1, rb1 = pr(Y)
    out[f'acc_{tag}_n_samples'] = len(samples); out[f'acc_{tag}_n_features'] = len(feats)
    out[f'acc_{tag}_spearman_stock_mean'] = float(ra.mean()); out[f'acc_{tag}_spearman_fix_mean'] = float(rb.mean())
    dd = rb - ra
    out[f'acc_{tag}_spearman_diff_mean'] = float(dd.mean())
    out[f'acc_{tag}_spearman_fix_better_worse_tie'] = [int((dd > 1e-12).sum()), int((dd < -1e-12).sum()), int((np.abs(dd) <= 1e-12).sum())]
    try:
        out[f'acc_{tag}_wilcoxon_p'] = float(wilcoxon(rb, ra).pvalue) if (np.abs(dd) > 1e-12).any() else 1.0
    except Exception:
        out[f'acc_{tag}_wilcoxon_p'] = None
    out[f'acc_{tag}_precision_stock_fix'] = [float(pa1.mean()), float(pb1.mean())]
    out[f'acc_{tag}_recall_stock_fix'] = [float(ra1.mean()), float(rb1.mean())]
    pd.DataFrame({'sample': samples, 'rho_stock': ra, 'rho_fix': rb, 'prec_stock': pa1, 'prec_fix': pb1,
                  'rec_stock': ra1, 'rec_fix': rb1}).to_csv(f'{RES}/{ds}_{tag}_accuracy_per_sample.tsv', sep='\t', index=False)

poss_ko = set(pd.read_csv('/opt/mm/root/envs/picrust2/lib/python3.12/site-packages/picrust2/default_files/bacteria/ko.txt.gz',
                          sep='\t', nrows=1, index_col=0).columns)
poss_ko = {c[3:] if c.startswith('ko:') else c for c in poss_ko}
poss_ko &= set(open(f'{MS}/16S_validation/possible_ko/humann2_ko.txt').read().split())
accuracy('KO', ko_a, ko_b, f'{MS}/mgs_validation/{ds}/humann2_ko_unstrat.tsv', poss_ko)
accuracy('PATH', pw_a, pw_b, f'{MS}/mgs_validation/{ds}/humann2_pathabun_unstrat.tsv',
         set(open(f'{MS}/16S_validation/possible_path/picrust2_path.txt').read().split()))
print(json.dumps(out))
