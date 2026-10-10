"""ML-scored constrained SPR (GTM-Blend-ML) from the published GTM tree.
Usage: python3 run_mlspr_pub.py COND REP GUIDE OUTDIR [radius]"""
import json, os, sys, time
from datasets import inputs, read_fasta
from phylo import read_tree, fn_fp, parse_newick, is_induced
import mlspr
cond, rep, guide, outdir = sys.argv[1:5]
rad = int(sys.argv[5]) if len(sys.argv) > 5 else 4
od = f"{outdir}/{cond}/{guide}/{rep}"
os.makedirs(od, exist_ok=True)
if os.path.exists(f"{od}/mlspr.json"):
    print(open(f"{od}/mlspr.json").read()); sys.exit()
x = inputs(cond, rep, guide)
T = read_tree(x["true"]); G = read_tree(x["published_gtm"])
subs = [read_tree(p) for p in x["subsets"]]
seqs = read_fasta(x["aln"]); names = sorted(G.label.values())
t0 = time.time()
logf = open(f"{od}/mlspr.log", "w")
bs, h = mlspr.run(G, names, seqs, [[t.label[v] for v in t.leaves()] for t in subs], f"{od}/work", radius=rad,
                  model_tree=x["guide"], max_rounds=int(os.environ.get("MLSPR_MAX_ROUNDS", "30")), log=lambda s: (logf.write(s + "\n"), logf.flush()))
B = parse_newick(bs.to_newick())
open(f"{od}/mlspr.tre", "w").write(bs.to_newick() + "\n")
res = dict(cond=cond, rep=rep, guide=guide, radius=rad, max_rounds=int(os.environ.get("MLSPR_MAX_ROUNDS", "30")), fn_gtm=fn_fp(G, T)[0], fn_mlspr=fn_fp(B, T)[0],
           moves=len(h) - 1, logL_gtm=h[0], logL_mlspr=h[-1], constraints_ok=all(is_induced(B, s) for s in subs),
           sec=time.time() - t0)
json.dump(res, open(f"{od}/mlspr.json", "w"), indent=1)
print(json.dumps(res), flush=True)
