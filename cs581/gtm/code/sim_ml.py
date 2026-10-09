"""ML-scored constrained SPR (mlspr.py) from GTM on one simulated replicate.
Usage: python3 sim_ml.py SIMDIR STMODE [radius]   (SIMDIR from sim.py; writes SIMDIR/sub_STMODE/mlspr.json)"""
import glob, json, os, sys, time
from phylo import read_tree, fn_fp, parse_newick, is_induced
from datasets import read_fasta
import mlspr
d, st = sys.argv[1], sys.argv[2]
rad = int(sys.argv[3]) if len(sys.argv) > 3 else 4
sd = f"{d}/sub_{st}"
if os.path.exists(f"{sd}/mlspr.json"):
    print(open(f"{sd}/mlspr.json").read()); sys.exit()
T = read_tree(f"{d}/true.tre"); seqs = read_fasta([f"{d}/aln.fa"]); names = sorted(seqs)
subs = [read_tree(p) for p in sorted(glob.glob(f"{sd}/s*.tre"))]
G = read_tree(f"{sd}/gtm.tre")
t0 = time.time()
wd = f"{sd}/mlspr_work"
bs, h = mlspr.run(G, names, seqs, [[t.label[v] for v in t.leaves()] for t in subs], wd, radius=rad,
                  model_tree=f"{d}/guide.tre", log=lambda s: None)
B = parse_newick(bs.to_newick())
open(f"{sd}/mlspr.tre", "w").write(bs.to_newick() + "\n")
mtxt = open(f"{wd}/model.raxml.bestModel").read().split(",")[0]
tl = mlspr.raxml_eval([open(f"{d}/true.tre").read().strip()], f"{wd}/full.fa", mtxt, f"{wd}/true")[0][0]
q = json.load(open(f"{sd}/quick.json")) if os.path.exists(f"{sd}/quick.json") else {}
res = dict(q, dir=d, subtrees=st, radius=rad, fn_gtm=fn_fp(G, T)[0], fn_mlspr=fn_fp(B, T)[0], moves=len(h) - 1,
           logL_gtm=h[0], logL_mlspr=h[-1], logL_true=tl, constraints_ok=all(is_induced(B, s) for s in subs),
           sec=time.time() - t0)
json.dump(res, open(f"{sd}/mlspr.json", "w"), indent=1)
print(json.dumps(res), flush=True)
