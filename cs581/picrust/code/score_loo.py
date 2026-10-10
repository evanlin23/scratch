#!/usr/bin/env python3
"""Score leave-out placements and KO predictions for stock vs patched EPA-ng.
For each held-out genome g:
  * true edge = edge of the pruned reference tree where g attached in the full tree (the edge above
    the clade formed by g's sibling leaves, after removing other held-out leaves)
  * placement error = number of edges between best-LWR edge and true edge (0 = exact)
  * KO accuracy = precision / recall / F1 of predicted KO presence (copy >0) vs g's true KO table row,
    and Spearman of copy numbers over the union of KOs present in either.
Usage: score_loo.py setdir"""
import sys, json, re, gzip, collections
import numpy as np, pandas as pd
from scipy.stats import spearmanr, wilcoxon
S = sys.argv[1]
B = '/opt/mm/root/envs/picrust2/lib/python3.12/site-packages/picrust2/default_files/bacteria'
sys.setrecursionlimit(10**6)

def parse(s):
    """minimal newick parser -> list of nodes {par, kids, label, edge}"""
    nodes = [{'par': None, 'kids': [], 'label': '', 'edge': None}]
    cur = 0; i = 0; n = len(s)
    while i < n:
        c = s[i]
        if c == '(':
            nodes.append({'par': cur, 'kids': [], 'label': '', 'edge': None}); nodes[cur]['kids'].append(len(nodes) - 1)
            cur = len(nodes) - 1; i += 1
        elif c == ',':
            p = nodes[cur]['par']
            nodes.append({'par': p, 'kids': [], 'label': '', 'edge': None}); nodes[p]['kids'].append(len(nodes) - 1)
            cur = len(nodes) - 1; i += 1
        elif c == ')':
            cur = nodes[cur]['par']; i += 1
        elif c == ';':
            break
        else:
            j = i
            while j < n and s[j] not in '(),;': j += 1
            m = re.match(r'([^:{\[]*)(?::[^{\[]*)?(?:\[[^\]]*\])?(?:\{(\d+)\})?', s[i:j])
            nodes[cur]['label'] = m.group(1)
            if m.group(2) is not None: nodes[cur]['edge'] = int(m.group(2))
            i = j
    return nodes

def leafsets(nodes):
    ls = [None] * len(nodes)
    order = []; st = [0]
    while st:
        x = st.pop(); order.append(x); st.extend(nodes[x]['kids'])
    for x in reversed(order):
        ls[x] = frozenset([nodes[x]['label']]) if not nodes[x]['kids'] else frozenset().union(*(ls[k] for k in nodes[x]['kids']))
    return ls

held = [l.strip() for l in open(f'{S}/heldout.txt') if l.strip()]
hs = set(held)
full = parse(open(f'{B}/bac_ref/bac_ref.tre').read().strip())
fls = leafsets(full)
lab2node = {nd['label']: i for i, nd in enumerate(full) if not nd['kids']}
# sibling clade of each held-out leaf in the full tree, minus held-out leaves; climb if empty
true_clade = {}
for g in held:
    x = lab2node[g]
    while True:
        p = full[x]['par']
        sib = frozenset().union(*(fls[k] for k in full[p]['kids'] if k != x)) - hs
        if sib: break
        x = p
    true_clade[g] = sib

res = {}
for v in ('stock', 'fix'):
    jp = json.load(open(f'{S}/int_{v}/epa_out/epa_result.jplace'))
    T = parse(jp['tree'].strip()); tls = leafsets(T)
    allleaves = tls[0]
    edge2node = {nd['edge']: i for i, nd in enumerate(T) if nd['edge'] is not None}
    # map clade (or its complement, unrooted) -> node
    clade2node = {}
    for i, l in enumerate(tls):
        clade2node[l] = i; clade2node[allleaves - l] = i
    depth = [0] * len(T)
    for i in range(1, len(T)): depth[i] = depth[T[i]['par']] + 1  # parents always precede kids
    def dist(a, b):
        d = 0
        while depth[a] > depth[b]: a = T[a]['par']; d += 1
        while depth[b] > depth[a]: b = T[b]['par']; d += 1
        while a != b: a = T[a]['par']; b = T[b]['par']; d += 2
        return d
    f = jp['fields']; ie, iw = f.index('edge_num'), f.index('like_weight_ratio')
    best = {}
    for p in jp['placements']:
        b = max(p['p'], key=lambda r: r[iw])
        for nme in (p.get('n') or [x[0] for x in p['nm']]):
            best[nme] = (b[ie], b[iw], len(p['p']))
    out = {}
    for g in held:
        if g not in best: continue
        tn = clade2node.get(true_clade[g])
        e, lwr, npl = best[g]
        out[g] = dict(edge=e, lwr=lwr, npl=npl, err=dist(edge2node[e], tn) if tn is not None else None)
    res[v] = out

# KO predictions
true = pd.read_csv(f'{B}/ko.txt.gz', sep='\t', index_col=0)
true = true.loc[[g for g in held if g in true.index]]
rows = []
for v in ('stock', 'fix'):
    pr = pd.read_csv(f'{S}/ko_{v}.tsv.gz', sep='\t', index_col=0)
    nsti = pr['metadata_NSTI']; pr = pr.drop(columns=['metadata_NSTI'])
    for g in held:
        if g not in pr.index or g not in res[v]: continue
        t = true.loc[g]; p = pr.loc[g].reindex(true.columns).fillna(0)
        tp = ((t > 0) & (p > 0)).sum(); P = (p > 0).sum(); Tn = (t > 0).sum()
        prec, rec = tp / max(P, 1), tp / max(Tn, 1)
        u = (t > 0) | (p > 0)
        rows.append(dict(genome=g, version=v, err=res[v][g]['err'], lwr=res[v][g]['lwr'], npl=res[v][g]['npl'],
                         edge=res[v][g]['edge'], nsti=nsti[g], precision=prec, recall=rec,
                         f1=2 * prec * rec / max(prec + rec, 1e-12), rho=spearmanr(t[u], p[u])[0]))
df = pd.DataFrame(rows)
df.to_csv(f'{S}/loo_scores.tsv', sep='\t', index=False)
w = df.pivot(index='genome', columns='version')
ch = w[('edge', 'stock')] != w[('edge', 'fix')]
summ = {'n': int(len(w)), 'placement_changed': int(ch.sum())}
for m in ('err', 'lwr', 'nsti', 'precision', 'recall', 'f1', 'rho'):
    a, b = w[(m, 'stock')].astype(float), w[(m, 'fix')].astype(float)
    summ[m] = {'stock_mean': float(a.mean()), 'fix_mean': float(b.mean()),
               'changed_stock_mean': float(a[ch].mean()) if ch.any() else None,
               'changed_fix_mean': float(b[ch].mean()) if ch.any() else None}
    d = (b - a)[ch]
    if m in ('err', 'nsti', 'f1', 'rho', 'precision', 'recall') and (d.abs() > 1e-12).sum() > 0:
        summ[m]['fix_better_worse_tie_on_changed'] = ([int((d < -1e-12).sum()), int((d > 1e-12).sum())] if m in ('err', 'nsti') else
                                                       [int((d > 1e-12).sum()), int((d < -1e-12).sum())]) + [int((d.abs() <= 1e-12).sum())]
        summ[m]['wilcoxon_p'] = float(wilcoxon(b[ch], a[ch]).pvalue)
summ['exact_edge_stock_fix'] = [int((w[('err', 'stock')] == 0).sum()), int((w[('err', 'fix')] == 0).sum())]
summ['lwr_ge_0.99_stock_fix'] = [int((w[('lwr', 'stock')] >= 0.99).sum()), int((w[('lwr', 'fix')] >= 0.99).sum())]
print(json.dumps(summ, indent=1))
