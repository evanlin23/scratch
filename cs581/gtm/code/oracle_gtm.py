"""Best possible UNBLENDED merger: GTM with the true tree as guide (GTM returns the
unblended merger closest in RF to its guide). Compare to the any-DTM floor (ceiling.py).
Usage: python3 oracle_gtm.py GTM_SRC OUTDIR OUT_TSV"""
import os, subprocess, sys
from datasets import REPS, inputs
from phylo import read_tree, fn_fp
gtm, outdir, out = sys.argv[1:4]
rows = []
for cond, reps in REPS.items():
    for rep in reps:
        for g in ("FT", "IQ"):
            x = inputs(cond, rep, g)
            o = f"{outdir}/{cond}/{g}/{rep}/gtm_trueguide.tre"
            os.makedirs(os.path.dirname(o), exist_ok=True)
            subprocess.run([sys.executable, f"{gtm}/gtm.py", "-s", x["true"], "-t", *x["subsets"], "-o", o],
                           check=True, capture_output=True)
            T = read_tree(x["true"])
            rows.append((cond, rep, g, fn_fp(read_tree(o), T)[0]))
with open(out, "w") as f:
    f.write("condition\treplicate\tguide\tFN_GTM_true_guide\n")
    for r in rows:
        f.write("%s\t%s\t%s\t%.5f\n" % r)
import collections
a = collections.defaultdict(list)
for r in rows: a[(r[0], r[2])].append(r[3])
for k, v in a.items(): print(k, "%.2f%%" % (100 * sum(v) / len(v)))
