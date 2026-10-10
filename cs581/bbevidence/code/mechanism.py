"""Does the change in evidence quality predict the change in MAGUS's error? Joins results and diag
rows (all variants vs the control on the same replicate) and reports correlations.

    python3 mechanism.py RESULTS_DIR   -> markdown on stdout
"""

import collections
import glob
import json
import os
import sys

import numpy as np
from scipy.stats import spearmanr

d = sys.argv[1]
res, diag = collections.defaultdict(dict), collections.defaultdict(dict)
for kind, store in (("results", res), ("diag", diag)):
    for f in glob.glob(os.path.join(d, "*.{}.jsonl".format(kind))):
        for line in open(f):
            r = json.loads(line)
            store[r["rep"]][r["variant"]] = r

rows = []
for rep in res:
    if "linsi" not in res[rep] or "linsi" not in diag[rep]:
        continue
    c, cd = res[rep]["linsi"], diag[rep]["linsi"]
    for v, r in res[rep].items():
        if v == "linsi" or v not in diag[rep] or "|" in v:  # masked variants: diag rows predate the fix
            continue
        g = diag[rep][v]
        rows.append({
            "rep": rep, "variant": v,
            "dErr": 100 * (r["avgErr"] - c["avgErr"]), "dFN": 100 * (r["SPFN"] - c["SPFN"]),
            "dFP": 100 * (r["SPFP"] - c["SPFP"]),
            "dPrec": 100 * (g["unit_prec"] - cd["unit_prec"]),
            "dTrueW": 100 * (g["w_true_edges"] / g.get("nbb_norm", 1) - cd["w_true_edges"]) / cd["w_true_edges"],
            "dFalseW": 100 * (g["w_false_edges"] - cd["w_false_edges"]) / cd["w_false_edges"],
            "dKeepFalse": 100 * (g["bb_keeps_subset_false"] - cd["bb_keeps_subset_false"]),
            "dFalseSurv": 100 * (g.get("false_w_survival", np.nan) - cd.get("false_w_survival", np.nan)),
            "nuc": not rep.startswith("BBA"),
        })

print("rows: {} (variant x replicate pairs, unmasked variants)".format(len(rows)))
print()
print("| predictor (variant - control) | Spearman with ΔSPFP | with ΔSPFN | with Δerror |")
print("|---|---|---|---|")
for p in ("dPrec", "dFalseW", "dTrueW", "dKeepFalse", "dFalseSurv"):
    xs = [r[p] for r in rows if not np.isnan(r[p])]
    out = []
    for y in ("dFP", "dFN", "dErr"):
        ys = [r[y] for r in rows if not np.isnan(r[p])]
        rho, pv = spearmanr(xs, ys)
        out.append("{:+.2f} (p={:.1g})".format(rho, pv))
    print("| {} | {} |".format(p, " | ".join(out)))

# per evidence-amount normalised: evidence units per backbone differ; report relative change in
# total true and false evidence weight as well (already relative above)
print()
print("| rep | variant | Δprec (pts) | Δfalse weight % | Δtrue weight % | ΔSPFP | ΔSPFN | Δerr |")
print("|---|---|---|---|---|---|---|---|")
for r in sorted(rows, key=lambda r: (r["rep"], r["dErr"])):
    if r["variant"] in ("clustalo", "linsi&clustalo", "linsi+clustalo", "linsi&fftns2", "fftns2", "famsa",
                        "clustalo&fftns2"):
        print("| {} | {} | {:+.2f} | {:+.1f} | {:+.1f} | {:+.2f} | {:+.2f} | {:+.2f} |".format(
            r["rep"].replace("_R0", ""), r["variant"], r["dPrec"], r["dFalseW"], r["dTrueW"], r["dFP"],
            r["dFN"], r["dErr"]))
