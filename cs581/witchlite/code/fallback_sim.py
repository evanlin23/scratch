#!/usr/bin/env python3
"""Offline test of a confidence-triggered fallback for beam-3 (score all HMMs if the best bit-score found is low).
usage: fallback_sim.py INST"""
import sys, json, os, numpy as np
sys.argv=[sys.argv[0]]+sys.argv[1:]
sys.path.insert(0,'/home/user/scratch/cs581/witchlite/code')
from lite import load_ensemble, adj, weights_from_scores
inst=sys.argv[1]
nodes,roots=load_ensemble(f'{inst}/witch_default/tree_decomp/root')
S={q:{int(i):b for i,b in d.items()} for q,d in json.load(open(f'{inst}/lite_all/allscores.json')).items()}
def beam(q,b):
    A=lambda i: adj(S[q][i],nodes[i]['size'])
    sel=set(roots); fr=sorted(roots,key=lambda r:-A(r))[:b]
    while True:
        ch=[c for f in fr for c in nodes[f]['children']]
        if not ch: return sel
        sel.update(ch); fr=sorted(ch,key=lambda c:-A(c))[:b]
Q=list(S); N=len(nodes)
full={q:weights_from_scores(S[q],nodes,10,1.0) for q in Q}
rows=[]
for q in Q:
    sel=beam(q,3)
    cov=sum(w for i,w in full[q] if i in sel)
    best=max(S[q][i] for i in sel)
    # length-normalised best bitscore (bits per residue)
    rows.append((q,cov,best))
cov=np.array([r[1] for r in rows]); best=np.array([r[2] for r in rows])
print(os.path.basename(inst), 'beam3 mean cov %.3f'%cov.mean(), 'queries cov<0.5: %d'%(cov<0.5).sum())
for thr in [20,40,60,80,100]:
    fb=best<thr
    c2=np.where(fb,1.0,cov)
    print(f'  fallback if best bits<{thr}: frac fallback {fb.mean():.3f}  cov {c2.mean():.3f}  low-cov left {(c2<0.5).sum()}  scores/q {(np.where(fb,N,36)).mean():.1f}')
