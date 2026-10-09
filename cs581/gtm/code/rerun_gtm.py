"""Rerun GTM (github.com/vlasmirnov/GTM) on the published guide + subset trees and
compare to (a) the published GTM tree (RF distance) and (b) the true tree (FN).
Usage: python3 rerun_gtm.py GTM_SRC OUTDIR OUT_TSV
"""
import os
import subprocess
import sys
import time
from datasets import REPS, inputs
from phylo import read_tree, fn_fp

gtm, outdir, out_tsv = sys.argv[1:4]
rows = []
for cond, reps in REPS.items():
    for rep in reps:
        for g in ("FT", "IQ"):
            x = inputs(cond, rep, g)
            T = read_tree(x["true"])
            pub = read_tree(x["published_gtm"])
            for mode in ("convex", "old"):
                o = f"{outdir}/{cond}/{g}/{rep}/gtm_{mode}.tre"
                os.makedirs(os.path.dirname(o), exist_ok=True)
                t0 = time.time()
                subprocess.run([sys.executable, f"{gtm}/gtm.py", "-s", x["guide"], "-t", *x["subsets"],
                                "-o", o, "-m", mode], check=True, capture_output=True)
                dt = time.time() - t0
                e = read_tree(o)
                fn = fn_fp(e, T)[0]
                rf_pub = fn_fp(e, pub)[2]  # #bipartitions of published GTM tree missing in ours
                pub_fn = fn_fp(pub, T)[0]
                rows.append((cond, rep, g, mode, fn, pub_fn, rf_pub, dt))
                print(*rows[-1], sep="\t", flush=True)

with open(out_tsv, "w") as f:
    f.write("condition\treplicate\tguide\tmode\tFN_rerun\tFN_published\tRF_to_published\tseconds\n")
    for r in rows:
        f.write("%s\t%s\t%s\t%s\t%.5f\t%.5f\t%d\t%.2f\n" % r)
