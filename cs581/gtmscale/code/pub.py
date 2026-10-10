"""GTM-Blend-FT on the published GTM-pipeline inputs (Park, Zaharias & Warnow 2021):
published guide + IQ-TREE subset trees + published GTM tree -> Blend-FT and the
unconstrained FastTree polish control.
Usage: python3 pub.py OUTROOT COND REP GUIDE   (GUIDE in FT, IQ)
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../gtm/code"))
from datasets import inputs, read_fasta  # noqa: E402
from phylo import read_tree, fn_fp, is_induced  # noqa: E402
from blend import resolve  # noqa: E402
from pipe import constraint_aln, timed, FT, FTOPT  # noqa: E402

out, cond, rep, guide = sys.argv[1:5]
if os.path.exists(f"{out}/skip_{cond}"):  # budget: drop a condition without restarting the queue
    sys.exit(0)
w = f"{out}/{cond}/{guide}/{rep}"
os.makedirs(w, exist_ok=True)
x = inputs(cond, rep, guide)
T = read_tree(x["true"])
res = dict(cond=cond, rep=rep, guide=guide)
seqs = read_fasta(x["aln"])
names = sorted(seqs)
if not os.path.exists(f"{w}/aln.fa"):
    with open(f"{w}/aln.fa", "w") as f:
        for n in names:
            f.write(">%s\n%s\n" % (n, seqs[n]))
subs = [read_tree(p) for p in x["subsets"]]
G = read_tree(x["published_gtm"])
open(f"{w}/gtm_bin.tre", "w").write(resolve(G).to_newick() + "\n")
res["ncol"] = constraint_aln(subs, names, f"{w}/constraints.fa")
res["fn_guide"] = fn_fp(read_tree(x["guide"]), T)[0]
res["fn_gtm"] = fn_fp(G, T)[0]
for arm, extra in (("blendft", ["-constraints", f"{w}/constraints.fa"]), ("polishft", [])):
    p = f"{w}/{arm}.tre"
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        with open(p, "w") as f:
            wall, mem, rc = timed([FT, *FTOPT, *extra, "-intree", f"{w}/gtm_bin.tre", f"{w}/aln.fa"], stdout=f)
        json.dump(dict(wall=wall, mem=mem, rc=rc), open(p + ".time", "w"))
    E = read_tree(p)
    res["fn_" + arm] = fn_fp(E, T)[0]
    res["induced_" + arm] = sum(is_induced(E, s) for s in subs)
    res.update({f"{k}_{arm}": v for k, v in json.load(open(p + ".time")).items()})
res["k"] = len(subs)
json.dump(res, open(f"{w}/result.json", "w"), indent=1)
os.remove(f"{w}/aln.fa")
print(json.dumps(res), flush=True)
