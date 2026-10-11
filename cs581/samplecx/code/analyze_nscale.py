"""Analyse nscale.jsonl: k_95 (P(species tree exactly correct) >= 0.95) vs n at f = 0.2, and the
mean FN at fixed k. Fits log k95 = a + b log n (power law) per shape/method.
Outputs results/nscale_k95.csv, results/nscale.png"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analyze_sc import k_at, COL

RES = os.path.join(os.path.dirname(__file__), "..", "results")


def main():
    d = pd.DataFrame([json.loads(l) for l in open("/opt/runs/nscale.jsonl")])
    d["ok"] = d.FN == 0
    g = d.groupby(["shape", "n", "method", "k"]).agg(p=("ok", "mean"), FN=("FN", "mean"), reps=("rep", "nunique")).reset_index()
    out = []
    for (s, n, m), x in g.groupby(["shape", "n", "method"]):
        x = x.sort_values("k")
        lab, v = k_at(x.k, x.p)
        # bootstrap CI over replicates
        sub = d[(d["shape"] == s) & (d.n == n) & (d.method == m)]
        P = sub.pivot_table(index="rep", columns="k", values="ok").dropna()
        rng = np.random.default_rng(0)
        bs = []
        for _ in range(300):
            idx = rng.integers(0, len(P), len(P))
            bs.append(k_at(P.columns, P.iloc[idx].mean(axis=0).values)[1])
        bs = np.array(bs, float)
        lo, hi = (np.nanpercentile(bs, [5, 95]) if np.isfinite(bs).mean() > 0.5 else (np.nan, np.nan))
        out.append(dict(shape=s, n=n, method=m, reps=len(P), k95=lab, k95_v=v, ci90_lo=lo, ci90_hi=hi,
                        frac_boot_reached=np.isfinite(bs).mean()))
    k = pd.DataFrame(out)
    k.to_csv(f"{RES}/nscale_k95.csv", index=False)
    fits = []
    for (s, m), x in k.dropna(subset=["k95_v"]).groupby(["shape", "method"]):
        x = x[x.k95_v > 25]
        if len(x) >= 3:
            b, a = np.polyfit(np.log(x.n), np.log(x.k95_v), 1)
            fits.append(dict(shape=s, method=m, exponent_in_n=round(b, 2), npts=len(x)))
    fits = pd.DataFrame(fits)
    fits.to_csv(f"{RES}/nscale_fits.csv", index=False)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for i, s in enumerate(["cat", "bal"]):
        ax = axs[i]
        for m in ["astral", "astrid", "njst"]:
            x = k[(k["shape"] == s) & (k.method == m)].dropna(subset=["k95_v"]).sort_values("n")
            if not len(x):
                continue
            ax.errorbar(x.n, x.k95_v, yerr=[np.nan_to_num((x.k95_v - x.ci90_lo).clip(lower=0)), np.nan_to_num((x.ci90_hi - x.k95_v).clip(lower=0))], color=COL[m], marker="o",
                        capsize=3, label=m)
        ax.set_xscale("log", base=2); ax.set_yscale("log")
        ax.set_xlabel("n taxa"); ax.set_ylabel("k_95 (P(tree exactly correct) = 0.95)")
        ax.set_title(f"{'caterpillar' if s == 'cat' else 'balanced'}, all internal branches f = 0.2 CU, true gene trees", fontsize=9)
        ax.legend()
    fig.tight_layout(); fig.savefig(f"{RES}/nscale.png", dpi=130)
    pd.set_option("display.width", 200)
    print(k.round(1).to_string()); print(fits.to_string())
    print(g[g.k.isin([100, 200, 400, 800])].pivot_table(index=["shape", "n", "k"], columns="method", values="FN").round(4).to_string())


if __name__ == "__main__":
    main()
