#!/usr/bin/env python3
"""Offline selection study: replay HMM-selection rules on the cached all-HMM bit-scores
(lite_all/allscores.json) without re-running hmmsearch. Bit-scores are deterministic, so a rule's
selected set and its WITCH weights are exactly what an online run would produce.

Reports per rule: HMM scores/query, recall of WITCH's top-1 HMM, and the fraction of WITCH's
(full-ensemble) top-10 weight mass covered by the rule's selected HMMs.
Writes weights_<rule>.txt (k=10) into OUTDIR for rules passed with --write.
usage: sim.py INST OUTDIR [--write rule1,rule2]
"""
import sys, os, json, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lite import load_ensemble, adj, weights_from_scores

ap = argparse.ArgumentParser()
ap.add_argument('inst'); ap.add_argument('outdir'); ap.add_argument('--write', default='')
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
H = f'{a.inst}/witch_default/tree_decomp/root'
nodes, roots = load_ensemble(H)
S = {q: {int(i): b for i, b in d.items()} for q, d in json.load(open(f'{a.inst}/lite_all/allscores.json')).items()}
bd = f'{a.inst}/blast_sens2' if os.path.exists(f'{a.inst}/blast_sens2/tophits5.json') else f'{a.inst}/blast_sens'
top5 = json.load(open(f'{bd}/tophits5.json'))
Q = list(S)

depth = {}
def dep(i):
    if i not in depth:
        p = nodes[i]['parent']; depth[i] = 0 if p is None else dep(p) + 1
    return depth[i]
leaf_of = {}
for i, n in nodes.items():
    for nm in n['names']:
        if nm not in leaf_of or dep(i) > dep(leaf_of[nm]):
            leaf_of[nm] = i

def path(i):
    out = []
    while i is not None:
        out.append(i); i = nodes[i]['parent']
    return out

def A(q, i):
    return adj(S[q][i], nodes[i]['size'])

def beam(q, b):
    sel = set(roots); fr = sorted(roots, key=lambda r: -A(q, r))[:b]
    while True:
        ch = [c for f in fr for c in nodes[f]['children']]
        if not ch:
            return sel
        sel.update(ch)
        fr = sorted(ch, key=lambda c: -A(q, c))[:b]

def blastpaths(q, m, sib):
    hits = top5.get(q, [])[:m]
    if not hits:
        return None
    sel = set()
    for h, _ in hits:
        for i in path(leaf_of[h[0] if isinstance(h, list) else h]):
            sel.add(i)
            if sib and nodes[i]['parent'] is not None:
                sel.update(nodes[nodes[i]['parent']]['children'])
    return sel

def hybrid(q, m, b, minbits):
    """BLAST path(s)+siblings when the top hit is strong, else beam-b descent"""
    hits = top5.get(q, [])
    if hits and hits[0][1] >= minbits:
        return blastpaths(q, m, True)
    return beam(q, b)

RULES = {
    'beam1': lambda q: beam(q, 1), 'beam2': lambda q: beam(q, 2), 'beam3': lambda q: beam(q, 3),
    'beam4': lambda q: beam(q, 4), 'beam6': lambda q: beam(q, 6),
    'bp1': lambda q: blastpaths(q, 1, False), 'bp1sib': lambda q: blastpaths(q, 1, True),
    'bp3': lambda q: blastpaths(q, 3, False), 'bp3sib': lambda q: blastpaths(q, 3, True),
    'bp5sib': lambda q: blastpaths(q, 5, True),
    'beam2+bp1sib': lambda q: beam(q, 2) | (blastpaths(q, 1, True) or set()),
    'beam3+bp3sib': lambda q: beam(q, 3) | (blastpaths(q, 3, True) or set()),
    'hyb_b3_100': lambda q: hybrid(q, 3, 3, 100), 'hyb_b4_150': lambda q: hybrid(q, 3, 4, 150),
}
allset = set(nodes)
full = {q: weights_from_scores(S[q], nodes, 10, 1.0) for q in Q}
print(f'{os.path.basename(a.inst)}: {len(nodes)} HMMs, {len(Q)} queries')
print(f"{'rule':16s} scores/q  top1-recall  top10-weight-covered")
res = {}
for name, f in RULES.items():
    ns, t1, cov = [], [], []
    sels = {}
    for q in Q:
        sel = f(q)
        if sel is None:
            sel = allset
        sels[q] = sel
        ns.append(len(sel))
        t1.append(full[q][0][0] in sel)
        cov.append(sum(w for i, w in full[q] if i in sel) / sum(w for _, w in full[q]))
    res[name] = dict(scores_per_query=float(np.mean(ns)), top1=float(np.mean(t1)), cov=float(np.mean(cov)))
    print(f"{name:16s} {np.mean(ns):8.1f}  {np.mean(t1):11.3f}  {np.mean(cov):10.3f}")
    if name in a.write.split(','):
        with open(f'{a.outdir}/weights_{name}.txt', 'w') as fo:
            for q in Q:
                w = weights_from_scores({i: S[q][i] for i in sels[q]}, nodes, 10, 1.0)
                fo.write('{}:{}\n'.format(q, tuple(w)))
json.dump(res, open(f'{a.outdir}/sim.json', 'w'), indent=1)
