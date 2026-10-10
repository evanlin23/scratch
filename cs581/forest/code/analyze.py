"""Tables, paired tests and plots for the Forest+DTM pilot.
usage: python analyze.py RESULTS_DIR   (reads main.jsonl [+ large.jsonl], writes summary.md + PNGs)"""
import sys, os, json
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R = sys.argv[1]
REG_NAME = {"U:0.005:0.05": "short U[0.005,0.05] (Kim et al.)", "U:0.05:0.1": "long U[0.05,0.1] (Kim et al.)",
            "U:0.1:0.4": "deep U[0.1,0.4]", "UH:2.0:1.0": "ultrametric h=2, lognormal(1) rates, K2P"}
METHODS = ["NJ", "BIONJ", "FastME", "FastTree", "Forest", "FGTM_NJ_sub", "FGTM_NJ_ind", "FGTM_FastME_sub",
           "FGTM_FastME_ind", "CompGTM_NJ", "DecGTM_NJ"]
LABEL = {"FGTM_NJ_sub": "Forest+GTM(NJ)", "FGTM_NJ_ind": "Forest+GTM(NJ, induced)",
         "FGTM_FastME_sub": "Forest+GTM(FastME)", "FGTM_FastME_ind": "Forest+GTM(FastME, induced)",
         "CompGTM_NJ": "Forest comps only+GTM(NJ)", "DecGTM_NJ": "Centroid dec.+GTM(NJ)"}
COMPARE = [("FGTM_NJ_sub", "NJ"), ("FGTM_FastME_sub", "FastME"), ("FGTM_NJ_sub", "FastME"),
           ("FGTM_NJ_sub", "CompGTM_NJ"), ("FGTM_NJ_sub", "DecGTM_NJ"), ("FGTM_NJ_sub", "FastTree")]


def load(fn):
    rows = [json.loads(l) for l in open(fn) if l.strip()]
    err = [r for r in rows if "error" in r]
    rows = [r for r in rows if "error" not in r]
    df = pd.DataFrame(rows)
    if "key" in df:
        df = df.drop_duplicates("key", keep="first")  # a resumed run can re-do in-flight replicates
    return df, err


def paired(df, a, b, n):
    d = (df[f"FN_{a}"] - df[f"FN_{b}"]).dropna()
    band = 0.5 / (n - 3)  # less than one split
    w = int((d < -band).sum()); l = int((d > band).sum()); t = len(d) - w - l
    p = 1.0 if np.allclose(d, 0) else wilcoxon(d, zero_method="wilcox", alternative="two-sided").pvalue
    return d.mean(), w, t, l, p, len(d)


def fmt(x, nd=3):
    return "" if pd.isna(x) else f"{x:.{nd}f}"


def main():
    out = []
    df, err = load(os.path.join(R, "main.jsonl"))
    out.append(f"Replicates: {len(df)} ok, {len(err)} errors.\n")
    for reg in REG_NAME:
        sub = df[df.regime == reg]
        if sub.empty:
            continue
        n = int(sub.n.iloc[0])
        out.append(f"\n### {REG_NAME[reg]}, n={n}\n")
        out.append("Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).\n")
        cols = ["k", "reps", "sat"] + [LABEL.get(m, m) for m in METHODS] + ["Forest #comp", "Forest false splits"]
        out.append("| " + " | ".join(cols) + " |")
        out.append("|" + "---|" * len(cols))
        for k, g in sub.groupby("k"):
            vals = [str(k), str(len(g)), fmt(g.frac_saturated.mean(), 2)]
            for m in METHODS:
                v = fmt(g[f"FN_{m}"].mean()) if f"FN_{m}" in g and g[f"FN_{m}"].notna().any() else "—"
                if m == "Forest":
                    v += f" [{fmt(g.FP_Forest.mean())}]"
                vals.append(v)
            vals += [fmt(g.forest_ncomp.mean(), 1), fmt(g.forest_false.mean(), 2)]
            out.append("| " + " | ".join(vals) + " |")
        out.append("\nPaired two-sided Wilcoxon on FN rate (A − B; negative = A better). "
                   "W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).\n")
        out.append("| k | A | B | mean diff | W/T/L | p |")
        out.append("|---|---|---|---|---|---|")
        for k, g in sub.groupby("k"):
            for a, b in COMPARE:
                if g[f"FN_{b}"].isna().all():
                    continue
                md, w, t, l, p, nn = paired(g, a, b, n)
                out.append(f"| {k} | {LABEL.get(a, a)} | {LABEL.get(b, b)} | {md:+.4f} | {w}/{t}/{l} | {p:.2g} |")
        # reproduction-style metrics (Kim et al.): P(no false split), P(fully resolved & correct)
        out.append("\nKim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).\n")
        out.append("| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |")
        out.append("|---|---|---|---|---|---|")
        for k, g in sub.groupby("k"):
            out.append(f"| {k} | {(g.FN_NJ == 0).mean():.2f} | {(g.forest_false == 0).mean():.2f} | "
                       f"{(g.FN_NJ == 0).mean():.2f} | {((g.forest_false == 0) & (g.forest_splits == n - 3)).mean():.2f} | "
                       f"{g.forest_splits.mean():.1f} / {n - 3} |")
    # runtime
    out.append("\n### Runtime (mean seconds per replicate, single core)\n")
    out.append("| regime | n | k | NJ | FastME | FastTree | Forest (grid search) | Forest+GTM(NJ) total | Centroid dec.+GTM(NJ) |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    for (reg, n, k), g in df.groupby(["regime", "n", "k"]):
        out.append(f"| {reg} | {n} | {k} | {g.t_NJ.mean():.2f} | {g.t_FastME.mean():.2f} | "
                   f"{fmt(g.t_FastTree.mean(), 1) if 't_FastTree' in g else ''} | {g.t_Forest.mean():.1f} | "
                   f"{g.t_FGTM_NJ_sub.mean():.1f} | {g.t_DecGTM_NJ.mean():.2f} |")
    lf = os.path.join(R, "large.jsonl")
    if os.path.exists(lf):
        dl, errl = load(lf)
        if not dl.empty:
            out.append(f"\n### Large trees (n=500; reduced grid; {len(dl)} replicates, {len(errl)} errors)\n")
            ms = ["NJ", "FastME", "FastTree", "Forest", "FGTM_NJ_sub", "FGTM_FastME_sub", "DecGTM_NJ"]
            out.append("| regime | k | reps | " + " | ".join(LABEL.get(m, m) for m in ms) + " | Forest #comp | t_Forest (s) |")
            out.append("|---|---|---|" + "---|" * (len(ms) + 2))
            for (reg, k), g in dl.groupby(["regime", "k"]):
                out.append(f"| {reg} | {k} | {len(g)} | " + " | ".join(fmt(g[f'FN_{m}'].mean()) for m in ms) +
                           f" | {g.forest_ncomp.mean():.0f} | {g.t_Forest.mean():.0f} |")
            out.append("\nPaired tests (n=500, all k pooled):\n")
            out.append("| A | B | mean diff | W/T/L | p | N |")
            out.append("|---|---|---|---|---|---|")
            for a, b in COMPARE:
                md, w, t, l, p, nn = paired(dl, a, b, 500)
                out.append(f"| {LABEL.get(a, a)} | {LABEL.get(b, b)} | {md:+.4f} | {w}/{t}/{l} | {p:.2g} | {nn} |")
    ff = os.path.join(R, "followup.jsonl")
    if os.path.exists(ff):
        du, erru = load(ff)
        if not du.empty:
            out.append(f"\n### Follow-up: saturation handling and FastME-guided controls ({len(du)} replicates, {len(erru)} errors)\n")
            cs = ["NJ_cap2", "NJ_cap1.2", "NJ_cap5", "NJ_pcap", "FastME_cap2", "FastME_cap1.2", "FastME_cap5",
                  "FastME_pcap", "FGTM_FastME", "CompGTM_FastME", "DecGTM_FastME_25", "DecGTM_FastME_50",
                  "FGTM_FastME_pcap"]
            out.append("| regime | k | reps | " + " | ".join(cs) + " |")
            out.append("|---|---|---|" + "---|" * len(cs))
            for (reg, k), g in du.groupby(["regime", "k"]):
                out.append(f"| {reg} | {k} | {len(g)} | " + " | ".join(fmt(g[f'FN_{c}'].mean()) for c in cs) + " |")
            out.append("\n| regime | k | A | B | mean diff | W/T/L | p |")
            out.append("|---|---|---|---|---|---|---|")
            for (reg, k), g in du.groupby(["regime", "k"]):
                for a, b in (("FGTM_FastME_pcap", "FastME_pcap"), ("FGTM_FastME", "FastME_cap2"),
                             ("FGTM_FastME", "CompGTM_FastME"), ("FGTM_FastME", "DecGTM_FastME_25"),
                             ("FastME_pcap", "FastME_cap2")):
                    md, w, t, l, pv, nn = paired(g, a, b, 100)
                    out.append(f"| {reg} | {k} | {a} | {b} | {md:+.4f} | {w}/{t}/{l} | {pv:.2g} |")
    open(os.path.join(R, "summary.md"), "w").write("\n".join(out) + "\n")
    # plots
    series = [("NJ", "#2a78d6", "o"), ("FastME", "#eb6834", "s"), ("FastTree", "#1baf7a", "^"),
              ("FGTM_NJ_sub", "#eda100", "D"), ("DecGTM_NJ", "#e87ba4", "v")]
    regs = [r for r in REG_NAME if r in set(df.regime)]
    fig, axes = plt.subplots(1, len(regs), figsize=(4.2 * len(regs), 3.8), sharey=False)
    for ax, reg in zip(np.atleast_1d(axes), regs):
        g = df[df.regime == reg].groupby("k")
        for m, c, mk in series:
            s = g[f"FN_{m}"].mean().dropna()
            ax.plot(s.index, s.values, color=c, marker=mk, ms=6, lw=2, label=LABEL.get(m, m))
        s = g["FN_Forest"].mean()
        ax.plot(s.index, s.values, color="#7a7a74", ls="--", lw=1.5, marker="x", ms=6, label="Forest (FN = unresolved)")
        ax.set_xscale("log")
        ax.set_title(REG_NAME[reg], fontsize=9)
        ax.set_xlabel("sequence length k")
        ax.grid(alpha=0.25, lw=0.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
    np.atleast_1d(axes)[0].set_ylabel("mean FN rate")
    h, l = np.atleast_1d(axes)[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=6, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(os.path.join(R, "fn_vs_k.png"), dpi=130)
    if os.path.exists(ff) and not du.empty:
        fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8))
        ser = [("FastME_cap2", "#2a78d6", "o", "-", "FastME, cap2 (main)"),
               ("FGTM_FastME", "#2a78d6", "o", "--", "Forest+GTM(FastME), cap2"),
               ("FastME_pcap", "#eb6834", "s", "-", "FastME, p clipped"),
               ("FGTM_FastME_pcap", "#eb6834", "s", "--", "Forest+GTM(FastME), p clipped")]
        for ax, reg in zip(axes, ["U:0.1:0.4", "UH:2.0:1.0"]):
            g = du[du.regime == reg].groupby("k")
            for c, col, mk, ls, lab in ser:
                v = g[f"FN_{c}"].mean()
                ax.plot(v.index, v.values, color=col, marker=mk, ls=ls, lw=2, ms=6, label=lab)
            ax.set_xscale("log"); ax.set_title(REG_NAME[reg], fontsize=9); ax.set_xlabel("sequence length k")
            ax.grid(alpha=0.25, lw=0.5)
            for sp in ("top", "right"):
                ax.spines[sp].set_visible(False)
        axes[0].set_ylabel("mean FN rate")
        h, l = axes[0].get_legend_handles_labels()
        fig.legend(h, l, loc="lower center", ncol=2, fontsize=8, frameon=False)
        fig.tight_layout(rect=(0, 0.14, 1, 1))
        fig.savefig(os.path.join(R, "saturation_followup.png"), dpi=130)
    print("\n".join(out))


if __name__ == "__main__":
    main()
