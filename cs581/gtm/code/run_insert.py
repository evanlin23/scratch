"""Constrained-insertion blending merger + optional constrained-SPR parsimony polish.
Usage: python3 run_insert.py COND REP GUIDE MODE OUTDIR [polish_radius]
MODE: pars | guide | oracle. Writes OUTDIR/COND/GUIDE/REP/ins_MODE[_spr].{tre,json}"""
import json, os, sys, time
from datasets import inputs, read_fasta
from phylo import read_tree, fn_fp, is_induced, parse_newick
from blend import BlendSearch, compress
import insert

cond, rep, guide, mode, outdir = sys.argv[1:6]
rad = int(sys.argv[6]) if len(sys.argv) > 6 else 0
x = inputs(cond, rep, guide)
T = read_tree(x["true"]); G = read_tree(x["published_gtm"]); GU = read_tree(x["guide"])
subs = [read_tree(p) for p in x["subsets"]]
seqs = read_fasta(x["aln"])
names = sorted(G.label.values())
states, w = compress(seqs, names)
ungapped = {n: sum(c not in "-?N" for c in seqs[n]) for n in names}
t0 = time.time()
ref = T if mode == "oracle" else GU
g = insert.run(names, states, w, subs, mode, ref_tree=ref, order_key=lambda n: -ungapped[n])
t_ins = time.time() - t0
I = parse_newick(g.to_newick())
od = f"{outdir}/{cond}/{guide}/{rep}"; os.makedirs(od, exist_ok=True)
res = dict(cond=cond, rep=rep, guide=guide, mode=mode, fn_gtm=fn_fp(G, T)[0], fn_ins=fn_fp(I, T)[0],
           ok_ins=all(is_induced(I, t) for t in subs), sec_ins=t_ins)
bs = BlendSearch(G, names, states, w, [[t.label[v] for v in t.leaves()] for t in subs], radius=1)
res["pars_gtm"] = bs.length
bi = BlendSearch(I, names, states, w, [[t.label[v] for v in t.leaves()] for t in subs], radius=max(rad, 1))
res["pars_ins"] = bi.length
open(f"{od}/ins_{mode}.tre", "w").write(g.to_newick() + "\n")
if rad:
    t1 = time.time()
    res["spr_moves"] = bi.run_parsimony()
    P = parse_newick(bi.to_newick())
    res.update(fn_spr=fn_fp(P, T)[0], pars_spr=bi.length, ok_spr=all(is_induced(P, t) for t in subs),
               sec_spr=time.time() - t1)
    open(f"{od}/ins_{mode}_spr.tre", "w").write(bi.to_newick() + "\n")
json.dump(res, open(f"{od}/ins_{mode}{'_spr' if rad else ''}.json", "w"), indent=1)
print(json.dumps(res), flush=True)
