"""Aggregate /opt/runs/gdlcons/curves/*.jsonl into results/curves_summary.md,
results/curves_runs.jsonl (copy) and error-vs-families plots."""
import collections
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
rows = []
for f in sorted(glob.glob("/opt/runs/gdlcons/curves/*.jsonl")):
    rows += [json.loads(l) for l in open(f)]
with open(os.path.join(RES, "curves_runs.jsonl"), "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")

agg = collections.defaultdict(list)
errs = collections.Counter()
for r in rows:
    key = (r["setting"], r["method"], r["nfam"])
    if "FN" in r:
        agg[key].append(r["FN"] / r["nI"])
    else:
        errs[(r["setting"], r["method"])] += 1

settings = sorted({k[0] for k in agg})
ATLAS = ["astral-pro", "astrid-multi", "astrid-disco", "astral-disco", "fastmulrfs", "stag"]
out = ["# Error vs number of gene families (mean RF/FN rate over replicates; true gene trees)", "",
       "Each cell: mean FN rate (n replicates). FN rate = missing true bipartitions / (n-3).", ""]
for s in settings:
    ks = sorted({k[2] for k in agg if k[0] == s})
    meths = sorted({k[1] for k in agg if k[0] == s}, key=lambda m: (m not in ATLAS, ATLAS.index(m) if m in ATLAS else 0, m))
    out.append("## %s" % s)
    out.append("")
    out.append("| method | " + " | ".join(str(k) for k in ks) + " |")
    out.append("|---|" + "---|" * len(ks))
    for m in meths:
        cells = []
        for k in ks:
            v = agg.get((s, m, k))
            cells.append("%.3f (%d)" % (sum(v) / len(v), len(v)) if v else "")
        out.append("| %s | %s |" % (m, " | ".join(cells)))
    out.append("")
if errs:
    out.append("Runs that failed (method error or timeout): " + ", ".join("%s/%s: %d" % (a, b, c) for (a, b), c in sorted(errs.items())))
open(os.path.join(RES, "curves_summary.md"), "w").write("\n".join(out) + "\n")

# plots: one panel per setting; atlas methods (top row) and ASTRAL-Pro error models (separate figure)
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MARK = ["o", "s", "^", "D", "v", "P", "X", "*"]


def panel_fig(methods_of, fname, title):
    n = len(settings)
    cols = min(4, n)
    rows_ = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows_, cols, figsize=(3.4 * cols, 2.8 * rows_), squeeze=False, sharey=True)
    for i, s in enumerate(settings):
        ax = axes[i // cols][i % cols]
        ms = methods_of(s)
        for j, m in enumerate(ms):
            ks = sorted(k[2] for k in agg if k[0] == s and k[1] == m)
            if not ks:
                continue
            ys = [sum(agg[(s, m, k)]) / len(agg[(s, m, k)]) for k in ks]
            ax.plot(ks, ys, color=COLORS[j % 8], marker=MARK[j % 8], ms=4, lw=1.5, label=m)
        ax.set_xscale("log")
        ax.set_title(s, fontsize=9)
        ax.grid(alpha=0.25, lw=0.5)
        ax.tick_params(labelsize=7)
        if i % cols == 0:
            ax.set_ylabel("species-tree FN rate", fontsize=8)
        ax.set_xlabel("gene families", fontsize=8)
    for i in range(n, rows_ * cols):
        axes[i // cols][i % cols].axis("off")
    hs, ls = [], []
    for ax in axes.flat:
        h, l = ax.get_legend_handles_labels()
        for a, b in zip(h, l):
            if b not in ls:
                hs.append(a)
                ls.append(b)
    fig.legend(hs, ls, loc="lower center", ncol=min(4, len(ls)), fontsize=7, frameon=False)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout(rect=(0, 0.08 + 0.03 * (len(ls) // 5), 1, 0.95))
    fig.savefig(os.path.join(RES, fname), dpi=130)


panel_fig(lambda s: [m for m in ATLAS], "atlas_curves.png", "Method atlas: error vs number of families")
EM = ["apro-fixed:true:0", "apro-fixed:ovl:0", "apro-fixed:rovl:0.3", "apro-fixed:rovl:1",
      "apro-fixed:flip:0.15", "apro-fixed:flip:0.3", "apro-fixed:d2s:0.5", "apro-fixed:d2s:1"]
settings_all = settings
settings = [s for s in settings_all if s.startswith("gdl")]
if settings:
    panel_fig(lambda s: EM, "apro_error_curves.png", "ASTRAL-Pro under rooting/tagging error models (GDL)")
print("\n".join(out))
