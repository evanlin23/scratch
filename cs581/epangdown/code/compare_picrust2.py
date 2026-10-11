"""Compare two PICRUSt2 runs (stock vs fixed EPA-ng): placements, NSTI, predicted ASV gene
content and predicted metagenome. Usage: python compare_picrust2.py <run_a> <run_b>"""
import sys, os, json, glob
import numpy as np, pandas as pd
from scipy.stats import spearmanr
a, b = sys.argv[1], sys.argv[2]


def best(jp):
    J = json.load(open(jp))
    f = J['fields']; ie, il = f.index('edge_num'), f.index('like_weight_ratio')
    out = {}
    for p in J['placements']:
        top = max(p['p'], key=lambda r: r[il])
        for n in p.get('n', []) + [x[0] for x in p.get('nm', [])]:
            out[n] = (top[ie], top[il])
    return out, J['tree']


for jp in sorted(glob.glob(f'{a}/intermediate/place_seqs*/*/epa_out/epa_result.jplace') +
                 glob.glob(f'{a}/intermediate/*/epa_out/epa_result.jplace')):
    rel = os.path.relpath(jp, a)
    pa, ta = best(jp); pb, tb = best(os.path.join(b, rel))
    same_tree = ta == tb
    common = set(pa) & set(pb)
    ch = [q for q in common if pa[q][0] != pb[q][0]]
    print(f'{rel}: same tree string {same_tree}; queries {len(common)}; best edge changed for {len(ch)} '
          f'({100*len(ch)/max(1,len(common)):.1f}%)')
for name in ('marker_predicted_and_nsti.tsv.gz',):
    fa = glob.glob(f'{a}/**/{name}', recursive=True)
    for f in fa:
        rel = os.path.relpath(f, a)
        A = pd.read_csv(f, sep='\t', index_col=0); B = pd.read_csv(os.path.join(b, rel), sep='\t', index_col=0)
        A, B = A.align(B, join='inner', axis=0)
        print(f'{rel}: NSTI changed for {(abs(A.metadata_NSTI - B.metadata_NSTI) > 1e-9).sum()} of {len(A)} ASVs; '
              f'mean NSTI {A.metadata_NSTI.mean():.4f} vs {B.metadata_NSTI.mean():.4f}')
for name in ('KO_predicted.tsv.gz', 'EC_predicted.tsv.gz'):
    for f in glob.glob(f'{a}/**/{name}', recursive=True):
        rel = os.path.relpath(f, a)
        A = pd.read_csv(f, sep='\t', index_col=0); B = pd.read_csv(os.path.join(b, rel), sep='\t', index_col=0)
        A, B = A.align(B, join='inner')
        rows = (A != B).any(axis=1)
        print(f'{rel}: predicted gene-family profile changed for {rows.sum()} of {len(A)} ASVs; '
              f'cells changed {(A != B).values.sum()} of {A.size}')
for name in ('KO_metagenome_out/pred_metagenome_unstrat.tsv.gz', 'EC_metagenome_out/pred_metagenome_unstrat.tsv.gz',
             'pathways_out/path_abun_unstrat.tsv.gz'):
    fa, fb = os.path.join(a, name), os.path.join(b, name)
    if not os.path.exists(fa):
        continue
    A = pd.read_csv(fa, sep='\t', index_col=0); B = pd.read_csv(fb, sep='\t', index_col=0)
    A, B = A.align(B, join='outer', fill_value=0)
    rel = (abs(A - B).sum() / A.sum())  # per-sample L1 relative difference
    rho = [spearmanr(A[c], B[c]).statistic for c in A.columns]
    print(f'{name}: {len(A)} features x {A.shape[1]} samples; per-sample relative L1 diff '
          f'median {rel.median():.4f} max {rel.max():.4f}; Spearman min {min(rho):.5f}; '
          f'features with any change {((A - B).abs() > 1e-9).any(axis=1).sum()}')
