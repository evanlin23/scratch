"""Analyse sc_sim.jsonl: P(recover) vs k, k_95(f) per method/shape/n, slope of log k95 vs log f.

k_95 is estimated per (shape, n, f, method) from the empirical success curve by log-linear
interpolation in k between the grid points that bracket 0.95 (">" = not reached on the grid).
Two targets: (a) per-branch: mean over internal branches of P(branch recovered);
(b) whole tree: P(FN = 0).
"""
import sys, os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = os.path.join(os.path.dirname(__file__), "..", "results")
COL = {"astral": "#1f77b4", "astrid": "#d62728", "njst": "#ff7f0e"}


def load(path="/opt/runs/sc_sim.jsonl"):
    rows = [json.loads(l) for l in open(path)]
    d = pd.DataFrame(rows)
    d["ok"] = d.FN == 0
    d["branch_rec"] = 1 - d.FN
    return d


def k_at(ks, ps, target=0.95):
    ks, ps = np.asarray(ks, float), np.asarray(ps, float)
    if ps[0] >= target:
        return f"<={int(ks[0])}", ks[0]
    for i in range(1, len(ks)):
        if ps[i] >= target:
            lk = np.log(ks[i - 1]) + (target - ps[i - 1]) / (ps[i] - ps[i - 1] + 1e-12) * (np.log(ks[i]) - np.log(ks[i - 1]))
            v = float(np.exp(lk))
            return f"{v:.0f}", v
    return f">{int(ks[-1])}", np.nan


def main():
    d = load()
    nrep = d.groupby(["shape", "n", "f"]).rep.nunique()
    g = d.groupby(["shape", "n", "f", "method", "k"]).agg(p_tree=("ok", "mean"), p_branch=("branch_rec", "mean"),
                                                         reps=("rep", "nunique"), sec=("sec", "mean")).reset_index()
    g.to_csv(f"{RES}/sc_sim_curves.csv", index=False)
    out = []
    for (s, n, f, m), x in g.groupby(["shape", "n", "f", "method"]):
        x = x.sort_values("k")
        lab_t, v_t = k_at(x.k, x.p_tree)
        lab_b, v_b = k_at(x.k, x.p_branch)
        out.append(dict(shape=s, n=n, f=f, method=m, reps=int(x.reps.min()), k95_tree=lab_t, k95_tree_v=v_t,
                        k95_branch=lab_b, k95_branch_v=v_b))
    k95 = pd.DataFrame(out)
    k95.to_csv(f"{RES}/sc_sim_k95.csv", index=False)
    # slopes: log k95 ~ a + b log f, per shape, n, method (only where finite at >= 3 f values)
    sl = []
    for (s, n, m), x in k95.groupby(["shape", "n", "method"]):
        for tgt in ["k95_tree_v", "k95_branch_v"]:
            y = x.dropna(subset=[tgt])
            y = y[y[tgt] > 10]
            if len(y) >= 3:
                b, a = np.polyfit(np.log(y.f), np.log(y[tgt]), 1)
                sl.append(dict(shape=s, n=n, method=m, target=tgt.split("_")[1], slope=round(b, 2),
                               c_f2=round(float(np.exp(np.mean(np.log(y[tgt]) + 2 * np.log(y.f)))), 1), npts=len(y)))
    pd.DataFrame(sl).to_csv(f"{RES}/sc_sim_slopes.csv", index=False)

    # Figure 1: P(tree correct) vs k, one panel per (shape, f), lines = methods, line style = n
    shapes, fs, ns = sorted(d["shape"].unique()), sorted(d.f.unique()), sorted(d.n.unique())
    fig, axs = plt.subplots(len(shapes), len(fs), figsize=(4 * len(fs), 3.4 * len(shapes)), squeeze=False, sharey=True)
    ls = {8: ":", 16: "-.", 32: "--", 64: "-"}
    for i, s in enumerate(shapes):
        for j, f in enumerate(fs):
            ax = axs[i][j]
            for m in ["astral", "astrid", "njst"]:
                for n in ns:
                    x = g[(g["shape"] == s) & (g.f == f) & (g.method == m) & (g.n == n)].sort_values("k")
                    if len(x):
                        ax.plot(x.k, x.p_tree, ls[n], color=COL[m], lw=1.4,
                                label=f"{m} n={n}" if (i == 0 and j == 0) else None)
            ax.axhline(0.95, color="grey", lw=0.6)
            ax.set_xscale("log"); ax.set_title(f"{s}, every internal branch f={f} CU", fontsize=9)
            ax.set_xlabel("genes k")
            if j == 0:
                ax.set_ylabel("P(species tree exactly correct)")
    fig.legend(loc="lower center", ncol=6, fontsize=7)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(f"{RES}/sc_curves_tree.png", dpi=130)

    # Figure 2: k95 (whole tree) vs n for each f, methods
    fig, axs = plt.subplots(1, len(shapes), figsize=(5.5 * len(shapes), 4), squeeze=False)
    for i, s in enumerate(shapes):
        ax = axs[0][i]
        for m in ["astral", "astrid", "njst"]:
            for f, mk in zip(fs, "osD^v"):
                x = k95[(k95["shape"] == s) & (k95.method == m) & (k95.f == f)].sort_values("n")
                ax.plot(x.n, x.k95_tree_v, "-" + mk, color=COL[m], label=f"{m} f={f}", ms=4)
        ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.set_xlabel("n taxa")
        ax.set_ylabel("k_95 (genes for P(tree correct) = 0.95)"); ax.set_title(s)
    axs[0][0].legend(fontsize=6, ncol=2)
    fig.tight_layout(); fig.savefig(f"{RES}/sc_k95_vs_n.png", dpi=130)

    # Figure 3: k95 per branch vs f (log-log) with f^-2 reference
    fig, axs = plt.subplots(1, len(shapes), figsize=(5.5 * len(shapes), 4), squeeze=False)
    for i, s in enumerate(shapes):
        ax = axs[0][i]
        for m in ["astral", "astrid", "njst"]:
            for n in ns:
                x = k95[(k95["shape"] == s) & (k95.method == m) & (k95.n == n)].sort_values("f")
                ax.plot(x.f, x.k95_tree_v, ls[n], color=COL[m], marker="o", ms=3, label=f"{m} n={n}")
        ff = np.array(fs); ax.plot(ff, 20 * ff ** -2, color="k", lw=0.8, label="20 f^-2")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("f (CU)"); ax.set_ylabel("k_95 (whole tree)")
        ax.set_title(s)
    axs[0][0].legend(fontsize=6, ncol=2)
    fig.tight_layout(); fig.savefig(f"{RES}/sc_k95_vs_f.png", dpi=130)
    print(nrep.to_string())
    print(k95.pivot_table(index=["shape", "n", "f"], columns="method", values="k95_tree", aggfunc="first").to_string())
    print(pd.DataFrame(sl).to_string())


if __name__ == "__main__":
    main()
