"""Paired tree-error comparison of MAGUS, MAGUS(Slow) and slow-soft-m3 alignments.

    python analyze_soft.py SOFT_RESULTS.jsonl TREES.jsonl OUT.md [--tie 0.1]

SOFT_RESULTS = /opt/runs/soft/out/results.jsonl (FastSP scores from gcmx.experiment);
TREES = runtrees.py output with keys NAME/variant (FastTree and/or IQ-TREE --fast).
"""
import argparse
import json

import pandas as pd
from scipy import stats

VAR = {"default": "MAGUS", "slow": "MAGUS(Slow)", "slow-soft-m3": "slow-soft-m3", "true": "TRUE"}


def wtl(d, tie):
    return "%d/%d/%d" % ((d < -tie).sum(), (d.abs() <= tie).sum(), (d > tie).sum())


def main():
    p = argparse.ArgumentParser()
    p.add_argument("soft")
    p.add_argument("trees")
    p.add_argument("out")
    p.add_argument("--tie", type=float, default=0.1)
    a = p.parse_args()
    s = pd.DataFrame([json.loads(l) for l in open(a.soft)])
    s = s[s.get("avgErr").notna()][["dataset", "variant", "SPFN", "SPFP", "avgErr", "Compression", "TC", "seconds"]]
    t = pd.DataFrame([json.loads(l) for l in open(a.trees)])
    t[["dataset", "variant"]] = t.key.str.split("/", expand=True)
    t["fn"] = 100 * t.fn_rate
    L = []
    P = L.append
    P("Tie band for tree FN: |diff| <= %.2f points. Negative diff = variant better than reference.\n" % a.tie)
    for method, g in t.groupby("method"):
        piv = g.pivot_table(index="dataset", columns="variant", values="fn")
        piv = piv.dropna(subset=[v for v in ("default", "slow", "slow-soft-m3") if v in piv])
        P("### Tree FN %% (%s), %d replicates\n" % ({"ft": "FastTree GTR+G", "iqfast": "IQ-TREE --fast GTR+G4"}.get(method, method), len(piv)))
        P("Mean FN: " + ", ".join("%s %.2f" % (VAR.get(v, v), piv[v].mean()) for v in piv.columns) + "\n")
        P("| comparison | n | mean diff (FN pts) | W/T/L | Wilcoxon p (two-sided) |")
        P("|---|---|---|---|---|")
        for x, y in (("slow-soft-m3", "default"), ("slow-soft-m3", "slow"), ("slow", "default"),
                     ("default", "true"), ("slow-soft-m3", "true")):
            if x in piv and y in piv:
                d = (piv[x] - piv[y]).dropna()
                pv = stats.wilcoxon(d).pvalue if (d != 0).any() else 1.0
                P("| %s − %s | %d | %+.2f | %s | %.3g |" % (VAR[x], VAR[y], len(d), d.mean(), wtl(d, a.tie), pv))
        P("")
        P("| replicate | " + " | ".join(VAR.get(v, v) for v in piv.columns) + " |")
        P("|---" * (len(piv.columns) + 1) + "|")
        for r, row in piv.iterrows():
            P("| %s | " % r + " | ".join("%.2f" % v for v in row.values) + " |")
        P("")
    sp = s.pivot_table(index="dataset", columns="variant", values=["avgErr", "SPFN", "SPFP", "Compression"])
    P("### Alignment criteria on the same replicates (means)\n")
    P("| variant | n | SPFN % | SPFP % | (SPFN+SPFP)/2 % | compression (est/true length) | merge seconds |")
    P("|---|---|---|---|---|---|---|")
    for v in ("default", "slow", "slow-soft-m3"):
        g = s[s.variant == v]
        P("| %s | %d | %.2f | %.2f | %.2f | %.3f | %.0f |" % (VAR[v], len(g), 100 * g.SPFN.mean(), 100 * g.SPFP.mean(),
                                                         100 * g.avgErr.mean(), g.Compression.mean(), g.seconds.mean()))
    d = 100 * (sp["avgErr"]["slow-soft-m3"] - sp["avgErr"]["default"]).dropna()
    P("\nslow-soft-m3 − MAGUS alignment error: %+.2f pts, W/T/L %s (tie 0.05), Wilcoxon p = %.3g\n"
      % (d.mean(), wtl(d, 0.05), stats.wilcoxon(d).pvalue))
    # criteria vs tree FN within replicate on this set (TRUE has SPFN = SPFP = 0, compression 1)
    import numpy as np
    import statsmodels.formula.api as smf
    sc = s.copy()
    sc["SPFN"] *= 100
    sc["SPFP"] *= 100
    sc["TC_err"] = 100 * (1 - sc.TC)
    sc["log_comp"] = np.log(sc.Compression)
    for method, g in t.groupby("method"):
        m = g.merge(sc, on=["dataset", "variant"], how="left")
        tr = m.variant == "true"
        m.loc[tr, ["SPFN", "SPFP", "TC_err", "log_comp"]] = [0.0, 0.0, 0.0, 0.0]
        m = m.dropna(subset=["SPFN", "fn"])
        k = m.groupby("dataset").variant.transform("count")
        m = m[k == k.max()]
        for label, d0 in (("incl. TRUE", m), ("estimated only", m[m.variant != "true"])):
            d = d0.copy()
            for c in ("fn", "SPFN", "SPFP", "TC_err", "log_comp"):
                d[c] = d0[c] - d0.groupby("dataset")[c].transform("mean")
            grp = pd.factorize(d.dataset)[0]
            P("#### %s: within-replicate regressions of tree FN (%s, %d replicates, n=%d)\n" % (method, label, d.dataset.nunique(), len(d)))
            P("| model | coefs | cluster-robust p | within R² |")
            P("|---|---|---|---|")
            for f in ("SPFN", "SPFP", "log_comp", "TC_err", "SPFN + SPFP", "SPFN + SPFP + log_comp"):
                fit = smf.ols("fn ~ " + f, data=d).fit(cov_type="cluster", cov_kwds={"groups": grp})
                ts = f.split(" + ")
                P("| %s | %s | %s | %.3f |" % (f, ", ".join("%.3f" % fit.params[x] for x in ts),
                                            ", ".join("%.2g" % fit.pvalues[x] for x in ts), fit.rsquared))
            P("")
    open(a.out, "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
