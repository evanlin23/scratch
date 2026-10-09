"""Run GTM-Blend (constrained parsimony SPR) from the published GTM tree.
Usage: python3 run_blend.py COND REP GUIDE OUTDIR [radius]
Writes OUTDIR/COND/GUIDE/REP/{blend.tre,blend.json}."""
import json
import os
import sys
import time
from datasets import inputs, read_fasta
from phylo import read_tree, fn_fp, is_induced
from blend import BlendSearch, compress

cond, rep, guide, outdir = sys.argv[1:5]
radius = int(sys.argv[5]) if len(sys.argv) > 5 else 6
x = inputs(cond, rep, guide)
T = read_tree(x["true"])
G = read_tree(x["published_gtm"])
subs = [read_tree(p) for p in x["subsets"]]
seqs = read_fasta(x["aln"])
names = sorted(G.label.values())
states, w = compress(seqs, names)
t0 = time.time()
bs = BlendSearch(G, names, states, w, [[t.label[v] for v in t.leaves()] for t in subs], radius=radius)
p0 = bs.length
fn0 = fn_fp(G, T)[0]
traj = []


def log(m, L, d):
    if m % 10 == 0:
        print(f"  move {m} pars {L} ({time.time() - t0:.0f}s)", flush=True)
    traj.append(L)


moves = bs.run_parsimony(log=log)
nwk = bs.to_newick()
od = f"{outdir}/{cond}/{guide}/{rep}"
os.makedirs(od, exist_ok=True)
with open(f"{od}/blend.tre", "w") as f:
    f.write(nwk + "\n")
B = read_tree(f"{od}/blend.tre")
ok = all(is_induced(B, t) for t in subs)
res = dict(cond=cond, rep=rep, guide=guide, radius=radius, moves=moves, pars_gtm=p0, pars_blend=bs.length,
           fn_gtm=fn0, fn_blend=fn_fp(B, T)[0], constraints_ok=ok, seconds=time.time() - t0,
           sites=int(w.sum()), patterns=len(w))
json.dump(res, open(f"{od}/blend.json", "w"), indent=1)
print(json.dumps(res))
