"""Fill the sample-complexity placeholders of REPORT.md from results/*.csv (run after the analyses).
Reads REPORT.template.md, writes REPORT.md."""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
REP = os.path.join(HERE, "..", "REPORT.md")
TPL = os.path.join(HERE, "..", "REPORT.template.md")


def md(df):
    cols = list(df.columns)
    out = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(x) for x in r.values) + " |")
    return "\n".join(out)


def main():
    s = open(TPL).read()
    k = pd.read_csv(f"{RES}/sc_sim_k95.csv")
    for shape, tag in [("cat", "__SCTABLE_CAT__"), ("bal", "__SCTABLE_BAL__")]:
        x = k[k["shape"] == shape].pivot_table(index=["n", "f"], columns="method", values="k95_tree", aggfunc="first").reset_index()
        x = x[["n", "f", "astral", "astrid", "njst"]]
        s = s.replace(tag, md(x))
    s = s.replace("__SCREPS__", f"{int(k.reps.min())}–{int(k.reps.max())}")
    # c = k95 * f^2 averaged in log space over f values with a finite k95 > 10
    for shape in ["cat", "bal"]:
        for m in ["astral", "astrid"]:
            vals = []
            for n in [8, 16, 32, 64]:
                y = k[(k["shape"] == shape) & (k.method == m) & (k.n == n)].dropna(subset=["k95_tree_v"])
                y = y[y.k95_tree_v > 10]
                vals.append(f"{np.exp(np.mean(np.log(y.k95_tree_v * y.f ** 2))):.1f} ({len(y)} f)" if len(y) else "–")
            s = s.replace(f"| __C_{shape}_{m}__ |", "| " + " | ".join(vals) + " |")
    ns = pd.read_csv(f"{RES}/nscale_k95.csv")
    t = ns.copy()
    t["k95 [90% CI]"] = [f"{r.k95} [{r.ci90_lo:.0f}–{r.ci90_hi:.0f}]" if np.isfinite(r.ci90_lo) else str(r.k95)
                         for r in t.itertuples()]
    t = t.pivot_table(index=["shape", "n"], columns="method", values="k95 [90% CI]", aggfunc="first").reset_index()
    s = s.replace("__NSTABLE__", md(t[["shape", "n", "astral", "astrid", "njst"]]))
    s = s.replace("__NSREPS__", f"{int(ns.reps.min())}–{int(ns.reps.max())}")
    fits = pd.read_csv(f"{RES}/nscale_fits.csv")
    s = s.replace("__NSEXP__", "; ".join(f"{r.shape} {r.method} n^{r.exponent_in_n}" for r in fits.itertuples()))
    open(REP, "w").write(s)
    left = [w for w in s.split() if w.startswith("__")]
    print("unfilled:", left)


if __name__ == "__main__":
    main()
