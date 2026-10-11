"""Union of published and rerun (soft-MAGUS) alignments on replicates present in both sets.

    python analyze_union.py RESULTS_DIR   (reads pub_table.csv written by analyze_pub.py,
                                           /opt/runs/soft/out/results.jsonl, trees_soft.jsonl)
Within-replicate regressions of FastTree FN on SPFN, SPFP, compression with up to 9 alignments
per replicate (TRUE, 5 published methods, MAGUS rerun, MAGUS(Slow) rerun, slow-soft-m3).
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.formula.api as smf  # noqa: E402
from scipy import stats  # noqa: E402

out = sys.argv[1]
pub = pd.read_csv(os.path.join(out, "pub_table.csv"))
pub = pub[["rid", "cond", "method", "fn", "SPFN", "SPFP", "Compression", "TC"]]
s = pd.DataFrame([json.loads(l) for l in open("/opt/runs/soft/out/results.jsonl")])
s = s[s.avgErr.notna()]
t = pd.DataFrame([json.loads(l) for l in open("/opt/runs/alncrit/trees_soft.jsonl")])
t = t[t.method == "ft"]
t[["dataset", "variant"]] = t.key.str.split("/", expand=True)
m = t.merge(s, on=["dataset", "variant"])
names = {"default": "MAGUS-rerun", "slow": "MAGUS(Slow)-rerun", "slow-soft-m3": "slow-soft-m3"}
soft = pd.DataFrame({"rid": m.dataset.str.replace("_", "/"), "cond": m.dataset.str.split("_").str[0],
                     "method": m.variant.map(names), "fn": 100 * m.fn_rate, "SPFN": 100 * m.SPFN,
                     "SPFP": 100 * m.SPFP, "Compression": m.Compression, "TC": m.TC})
both = set(pub.rid) & set(soft.rid)
u = pd.concat([pub[pub.rid.isin(both)], soft[soft.rid.isin(both)]])
u["log_comp"] = np.log(u.Compression)
u["TC_err"] = 100 * (1 - u.TC)
u.to_csv(os.path.join(out, "union_table.csv"), index=False)
L = ["Replicates in both sets: %d (%s); alignments per replicate: %s\n" % (
    len(both), ", ".join(sorted(both)), ", ".join(sorted(u.method.unique())))]
for label, d0 in (("incl. TRUE", u), ("estimated only", u[u.method != "TRUE"])):
    d = d0.copy()
    for c in ("fn", "SPFN", "SPFP", "log_comp", "TC_err"):
        d[c] = d0[c] - d0.groupby("rid")[c].transform("mean")
    g = pd.factorize(d.rid)[0]
    L.append("#### Union, %s (n=%d)\n" % (label, len(d)))
    L.append("| model | coefs | cluster-robust p | within R² |")
    L.append("|---|---|---|---|")
    for f in ("SPFN", "SPFP", "log_comp", "TC_err", "SPFN + SPFP", "SPFN + SPFP + log_comp"):
        fit = smf.ols("fn ~ " + f, data=d).fit(cov_type="cluster", cov_kwds={"groups": g})
        ts = f.split(" + ")
        L.append("| %s | %s | %s | %.3f |" % (f, ", ".join("%.3f" % fit.params[x] for x in ts),
                                           ", ".join("%.2g" % fit.pvalues[x] for x in ts), fit.rsquared))
    r = stats.pearsonr(d.SPFN, d.SPFP)[0]
    L.append("\nwithin-replicate correlation of SPFN and SPFP: r = %.3f; SPFN vs log compression r = %.3f\n"
             % (r, stats.pearsonr(d.SPFN, d.log_comp)[0]))
    if label == "estimated only":
        fig, ax = plt.subplots(1, 2, figsize=(10, 4))
        for i, c in enumerate(("SPFP", "log_comp")):
            for meth, gg in d.groupby(d0.method):
                ax[i].scatter(gg[c], gg.fn, s=14, label=meth if i == 0 else None)
            ax[i].set_xlabel("%s (demeaned within replicate)" % c)
            ax[i].set_title("r = %.2f" % stats.pearsonr(d[c], d.fn)[0])
        ax[0].set_ylabel("FastTree FN % (demeaned)")
        fig.legend(loc="center right", fontsize=7, frameon=False)
        fig.tight_layout(rect=(0, 0, 0.82, 1))
        fig.savefig(os.path.join(out, "union_scatter.png"), dpi=130)
open(os.path.join(out, "union.md"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
