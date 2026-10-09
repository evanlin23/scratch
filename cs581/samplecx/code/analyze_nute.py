"""Analyse nute_curves.jsonl (Nute et al. complete data, 25+1 taxa, true / RAxML gene trees).

Outputs (results/):
  nute_fn_vs_k.png / .csv : mean species-tree FN vs genes, per ILS level, method, gene-tree type
  nute_branch_bins.csv    : P(branch recovered) by CU-length bin and k
  nute_k95.csv            : k_95(f) from a logistic fit logit P = b0 + b1 log k + b2 log f
                            (fit on branches with f < 1 CU), plus slope d log k95 / d log f = -b2/b1
  nute_k95.png
"""
import os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import minimize

RES = os.path.join(os.path.dirname(__file__), "..", "results")
COL = {"astral": "#1f77b4", "astrid": "#d62728", "njst": "#ff7f0e"}
ILS = {"10M": "low ILS (10M)", "2M": "high ILS (2M)", "500K": "very high ILS (500K)"}
BINS = [0, 0.01, 0.03, 0.1, 0.3, 1, np.inf]


def logistic_fit(k, f, y):
    X = np.column_stack([np.ones_like(k), np.log(k), np.log(f)])

    def nll(b):
        z = X @ b
        return np.sum(np.logaddexp(0, z) - y * z)

    def grad(b):
        z = X @ b
        return X.T @ (1 / (1 + np.exp(-z)) - y)
    r = minimize(nll, np.zeros(3), jac=grad, method="BFGS")
    # bootstrap-free SE from inverse Hessian approx
    return r.x, r.hess_inv


def main():
    rows = [json.loads(l) for l in open("/opt/runs/nute_curves.jsonl")]
    d = pd.DataFrame([{k: v for k, v in r.items() if k != "branches"} for r in rows])
    fnk = d.groupby(["height", "genetrees", "method", "k"]).agg(FN=("FN", "mean"), se=("FN", "sem"),
                                                                reps=("rep", "count")).reset_index()
    fnk.to_csv(f"{RES}/nute_fn_vs_k.csv", index=False)
    fig, axs = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    for i, gt in enumerate(["true", "raxml"]):
        for j, h in enumerate(["10M", "2M", "500K"]):
            ax = axs[i][j]
            for m in ["astral", "astrid", "njst"]:
                x = fnk[(fnk.height == h) & (fnk.genetrees == gt) & (fnk.method == m)].sort_values("k")
                ax.errorbar(x.k, x.FN, yerr=x.se, color=COL[m], label=m, marker="o", ms=3, capsize=2)
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_title(f"{ILS[h]}, {'true' if gt == 'true' else 'RAxML'} gene trees", fontsize=9)
            if i == 1: ax.set_xlabel("genes k")
            if j == 0: ax.set_ylabel("species-tree FN rate")
    axs[0][0].legend()
    fig.tight_layout(); fig.savefig(f"{RES}/nute_fn_vs_k.png", dpi=130)

    # branch-level table
    br = []
    for r in rows:
        for key, (f, ok) in r["branches"].items():
            br.append((r["height"], r["rate"], r["rep"], r["genetrees"], r["method"], r["k"], f, ok))
    b = pd.DataFrame(br, columns=["height", "rate", "rep", "genetrees", "method", "k", "f", "ok"])
    b["fbin"] = pd.cut(b.f, BINS, right=False)
    tab = b.groupby(["genetrees", "method", "fbin", "k"], observed=True).ok.agg(["mean", "count"]).reset_index()
    tab.to_csv(f"{RES}/nute_branch_bins.csv", index=False)

    out = []
    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    for i, gt in enumerate(["true", "raxml"]):
        ax = axs[i]
        for m in ["astral", "astrid", "njst"]:
            s = b[(b.genetrees == gt) & (b.method == m) & (b.f < 1) & (b.f > 0)]
            beta, H = logistic_fit(s.k.values.astype(float), s.f.values, s.ok.values.astype(float))
            slope = -beta[2] / beta[1]
            fs = np.array([0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5])
            k95 = np.exp((np.log(0.95 / 0.05) - beta[0] - beta[2] * np.log(fs)) / beta[1])
            for f, kk in zip(fs, k95):
                out.append(dict(genetrees=gt, method=m, f=f, k95=kk, slope=slope, b0=beta[0], b1=beta[1], b2=beta[2],
                                nobs=len(s)))
            ax.plot(fs, k95, "-o", color=COL[m], ms=3, label=f"{m} (slope {slope:.2f})")
            # empirical: per bin, smallest k with P >= 0.95
            for (lo, hi) in [(0.01, 0.03), (0.03, 0.1), (0.1, 0.3), (0.3, 1)]:
                t = s[(s.f >= lo) & (s.f < hi)].groupby("k").ok.mean()
                hit = t[t >= 0.95]
                if len(hit):
                    ax.plot(np.sqrt(lo * hi), hit.index.min(), "x", color=COL[m], ms=7)
        ax.plot(fs, 10 * fs ** -2.0, "k:", lw=0.8, label="10 f^-2")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("branch length f (CU)")
        ax.set_ylabel("k_95 (genes)"); ax.set_title(f"Nute et al. data, {gt} gene trees (lines: logistic fit; x: empirical bins)", fontsize=8)
        ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f"{RES}/nute_k95.png", dpi=130)
    k95 = pd.DataFrame(out); k95.to_csv(f"{RES}/nute_k95.csv", index=False)
    pd.set_option("display.width", 200)
    print(d.groupby(["genetrees", "method"]).rep.count())
    print(fnk.pivot_table(index=["height", "genetrees", "k"], columns="method", values="FN").round(3).to_string())
    print(k95.pivot_table(index=["genetrees", "f"], columns="method", values="k95").round(0).to_string())
    print(k95.groupby(["genetrees", "method"]).slope.first().round(2))
    print(tab[tab.fbin.astype(str).str.startswith("[0.01")].pivot_table(index=["genetrees", "k"], columns="method", values="mean").round(3))


if __name__ == "__main__":
    main()
