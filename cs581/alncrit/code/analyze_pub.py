"""Alignment criteria vs FastTree tree error on the MAGUS paper's published alignments.

    python analyze_pub.py TREES.jsonl SCORES.jsonl ALNSTATS.jsonl OUTDIR

SCORES = cs581/experiments/validation/published_scores.jsonl (FastSP on the same files).
Writes OUTDIR/pub_table.csv, OUTDIR/criteria.md and plots.
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

NAMES = {"true_align": "TRUE", "gcm": "MAGUS", "gcm_slow": "MAGUS(Slow)", "pasta_align": "PASTA(3)",
         "pasta_1_align": "PASTA(1)", "pasta_3_gcm_align": "PASTA(3)+GCM"}
CRIT = ["SPFN", "SPFP", "avgErr", "TC_err", "comp_dev", "log_comp", "indel_ratio"]


def load(trees, scores, alnstats):
    rows = []
    for r in map(json.loads, open(trees)):
        cond, rep, m = r["key"].split("/")
        rows.append({"cond": cond, "rep": rep, "file": m + ".txt", "method": NAMES.get(m, m),
                     "fn": 100 * r["fn_rate"], "fp": 100 * r["fp_rate"], "tree_s": r["seconds"], "aln": r["aln"]})
    t = pd.DataFrame(rows)
    s = pd.DataFrame([json.loads(l) for l in open(scores)])
    s = s[["dataset", "rep", "file", "SPFN", "SPFP", "TC", "Compression"]].rename(columns={"dataset": "cond"})
    df = t.merge(s, on=["cond", "rep", "file"], how="left")
    tr = df["method"] == "TRUE"
    df.loc[tr, ["SPFN", "SPFP", "TC", "Compression"]] = [0.0, 0.0, 1.0, 1.0]
    a = pd.DataFrame([json.loads(l) for l in open(alnstats)])
    df = df.merge(a[["aln", "indels", "gapfrac"]], on="aln", how="left")
    true_ind = df[tr].set_index(["cond", "rep"])["indels"]
    df["indel_ratio"] = [i / true_ind.get((c, r), np.nan) for c, r, i in zip(df.cond, df.rep, df.indels)]
    df["SPFN"] *= 100
    df["SPFP"] *= 100
    df["avgErr"] = (df.SPFN + df.SPFP) / 2
    df["TC_err"] = 100 * (1 - df.TC)
    df["log_comp"] = np.log(df.Compression)
    df["comp_dev"] = 100 * (1 - df.Compression)  # % shorter than true
    df["rid"] = df.cond + "/" + df.rep
    # keep replicates with all methods present
    n = df.groupby("rid")["method"].transform("count")
    return df[n == n.max()].copy()


def demean(df, cols):
    out = df.copy()
    for c in cols:
        out[c] = df[c] - df.groupby("rid")[c].transform("mean")
    return out


def partial(d, y, x, z):
    """Partial correlation of y and x controlling for z (list), on already-demeaned data."""
    X = np.column_stack([np.ones(len(d))] + [d[c] for c in z])
    ry = d[y] - X @ np.linalg.lstsq(X, d[y], rcond=None)[0]
    rx = d[x] - X @ np.linalg.lstsq(X, d[x], rcond=None)[0]
    return stats.pearsonr(rx, ry)


def main():
    trees, scores, alnst, outdir = sys.argv[1:5]
    os.makedirs(outdir, exist_ok=True)
    df = load(trees, scores, alnst)
    df.to_csv(os.path.join(outdir, "pub_table.csv"), index=False)
    L = []
    P = L.append
    nrep = df.rid.nunique()
    P("Replicates with all %d alignments: %d (%s)\n" % (df.method.nunique(), nrep,
      ", ".join("%s:%d" % (c, k) for c, k in df.groupby("cond").rid.nunique().items())))

    # 1. per-condition means
    P("### Mean per condition (FN = FastTree missing-branch rate %, SP errors %, compression = est/true length)\n")
    order = [m for m in NAMES.values() if m in set(df.method)]
    for col, fmt in (("fn", "%.2f"), ("SPFN", "%.2f"), ("SPFP", "%.2f"), ("Compression", "%.3f")):
        tab = df.pivot_table(index="method", columns="cond", values=col, aggfunc="mean").reindex(order)
        tab["all"] = df.groupby("method")[col].mean()
        P("**%s**\n" % col)
        P("| method | " + " | ".join(tab.columns) + " |")
        P("|---" * (len(tab.columns) + 1) + "|")
        for m, row in tab.iterrows():
            P("| %s | " % m + " | ".join(fmt % v for v in row.values) + " |")
        P("")

    # 2. paired: each method vs TRUE and vs MAGUS
    P("### Paired tree error vs MAGUS (FN points; negative = better tree than MAGUS; tie band |d| <= 0.1)\n")
    P("| method | n | mean diff | W/T/L | Wilcoxon p |")
    P("|---|---|---|---|---|")
    piv = df.pivot_table(index="rid", columns="method", values="fn")
    for m in order:
        if m == "MAGUS":
            continue
        d = piv[m] - piv["MAGUS"]
        w, t, l = (d < -0.1).sum(), (d.abs() <= 0.1).sum(), (d > 0.1).sum()
        p = stats.wilcoxon(d).pvalue if (d != 0).any() else 1.0
        P("| %s | %d | %+.2f | %d/%d/%d | %.3g |" % (m, len(d), d.mean(), w, t, l, p))
    P("")

    # 3. within-replicate correlations (rep fixed effects)
    est = df[df.method != "TRUE"]
    for label, data in (("all alignments incl. TRUE", df), ("estimated alignments only", est)):
        d = demean(data, ["fn"] + CRIT)
        P("### Within-replicate association with tree FN (%s; values demeaned within replicate, n=%d)\n"
          % (label, len(d)))
        P("| criterion | Pearson r | p | Spearman rho | p |")
        P("|---|---|---|---|---|")
        for c in CRIT:
            ok = d[c].notna()
            r, p = stats.pearsonr(d[c][ok], d.fn[ok])
            rho, p2 = stats.spearmanr(d[c][ok], d.fn[ok])
            P("| %s | %.3f | %.2g | %.3f | %.2g |" % (c, r, p, rho, p2))
        r1, p1 = partial(d, "fn", "SPFN", ["SPFP"])
        r2, p2 = partial(d, "fn", "SPFP", ["SPFN"])
        r3, p3 = partial(d, "fn", "log_comp", ["SPFN", "SPFP"])
        P("\nPartial r(FN, SPFN | SPFP) = %.3f (p=%.2g); partial r(FN, SPFP | SPFN) = %.3f (p=%.2g); "
          "partial r(FN, log compression | SPFN, SPFP) = %.3f (p=%.2g)\n" % (r1, p1, r2, p2, r3, p3))
        # fixed-effects regressions, cluster-robust by replicate
        P("| model (FN ~ ... + replicate FE) | coef (FN pts per criterion pt) | cluster-robust p | within R² |")
        P("|---|---|---|---|")
        groups = pd.factorize(data.rid)[0]
        for f in ("SPFN", "SPFP", "avgErr", "TC_err", "log_comp", "SPFN + SPFP", "SPFN + SPFP + log_comp",
                  "SPFN + SPFP + TC_err"):
            fit = smf.ols("fn ~ %s" % f, data=d).fit(cov_type="cluster", cov_kwds={"groups": groups})
            terms = f.split(" + ")
            P("| %s | %s | %s | %.3f |" % (f, ", ".join("%.3f" % fit.params[t] for t in terms),
                                        ", ".join("%.2g" % fit.pvalues[t] for t in terms), fit.rsquared))
        P("")

    # 4. per-condition within-replicate correlations
    P("### Per condition: within-replicate Pearson r of FN with each criterion (estimated alignments only)\n")
    P("| condition | n | " + " | ".join(CRIT) + " |")
    P("|---" * (len(CRIT) + 2) + "|")
    for c, g in est.groupby("cond"):
        d = demean(g, ["fn"] + CRIT)
        vals = []
        for k in CRIT:
            ok = d[k].notna() & (d[k].std() > 0)
            vals.append("%.2f" % stats.pearsonr(d[k][ok], d.fn[ok])[0] if ok.sum() > 3 else "–")
        P("| %s | %d | " % (c, len(d)) + " | ".join(vals) + " |")
    P("")

    # 5. between-condition / raw (no FE) association, the usual way such plots are made
    P("### Without replicate fixed effects (raw pooled values, estimated alignments)\n")
    for c in CRIT:
        ok = est[c].notna()
        rho, p = stats.spearmanr(est[c][ok], est.fn[ok])
        P("- %s: Spearman rho = %.3f (p=%.2g)" % (c, rho, p))
    gm = est.groupby(["cond", "method"])[["fn"] + CRIT].mean()
    P("\nCondition×method means (n=%d): Spearman rho(FN, avgErr) = %.3f, rho(FN, SPFN) = %.3f, rho(FN, SPFP) = %.3f\n"
      % (len(gm), stats.spearmanr(gm.fn, gm.avgErr)[0], stats.spearmanr(gm.fn, gm.SPFN)[0],
         stats.spearmanr(gm.fn, gm.SPFP)[0]))
    excess = (df.pivot_table(index="rid", columns="method", values="fn").sub(
        df[df.method == "TRUE"].set_index("rid").fn, axis=0)).drop(columns="TRUE")
    P("Mean FN excess over the TRUE-alignment tree, all conditions: " +
      ", ".join("%s %+.2f" % (m, excess[m].mean()) for m in order if m != "TRUE") + "\n")

    open(os.path.join(outdir, "criteria.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))

    # plots
    d = demean(est, ["fn"] + CRIT)
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    conds = sorted(d.cond.unique())
    cmap = plt.get_cmap("tab20")
    for i, c in enumerate(("SPFN", "SPFP", "log_comp")):
        for j, k in enumerate(conds):
            g = d[d.cond == k]
            ax[i].scatter(g[c], g.fn, s=12, color=cmap(j % 20), label=k if i == 0 else None, alpha=0.8)
        r = stats.pearsonr(d[c], d.fn)[0]
        ax[i].set_xlabel("%s (demeaned within replicate)" % c)
        ax[i].set_title("r = %.2f" % r)
        ax[i].axhline(0, color="#999", lw=0.5)
        ax[i].axvline(0, color="#999", lw=0.5)
    ax[0].set_ylabel("tree FN % (demeaned within replicate)")
    fig.legend(loc="center right", fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, 0, 0.92, 1))
    fig.savefig(os.path.join(outdir, "within_rep_scatter.png"), dpi=130)

    fig, ax = plt.subplots(figsize=(11, 4))
    tab = df.pivot_table(index="cond", columns="method", values="fn", aggfunc="mean")[order]
    tab.plot.bar(ax=ax, width=0.8)
    ax.set_ylabel("FastTree FN %")
    ax.set_xlabel("")
    ax.legend(fontsize=7, ncol=6)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "fn_by_method.png"), dpi=130)


if __name__ == "__main__":
    main()
