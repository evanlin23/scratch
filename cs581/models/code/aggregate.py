"""Collect results.json files into results/results.csv, a per-cell summary, and plots.

    python aggregate.py SIMDIR OUTDIR
"""
import csv
import glob
import json
import os
import statistics
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from characterize import aln_stats, read_fasta  # noqa: E402

METHODS = ["famsa", "mafft-auto", "mafft-linsi", "pasta", "magus", "true"]
AXES = {
    "clock (log-normal sigma)": ["clock0", "base", "clock0.7", "clock1.2", "clock2.0"],
    "indel length distribution": ["base", "pow1.7", "pow1.5"],
    "tree shape": ["base", "balanced", "caterpillar"],
}
LABEL = {"clock0": "0", "base": "base", "clock0.7": "0.7", "clock1.2": "1.2", "clock2.0": "2.0",
         "pow1.7": "Zipf 1.7", "pow1.5": "Zipf 1.5", "balanced": "balanced", "caterpillar": "caterpillar"}
COLORS = {"famsa": "#8c6bb1", "mafft-auto": "#999999", "mafft-linsi": "#e6550d", "pasta": "#3182bd",
          "magus": "#31a354", "true": "#000000"}


def main():
    simdir, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    rows = []
    for rj in sorted(glob.glob(os.path.join(simdir, "*_n*_r*", "results.json"))):
        d = os.path.dirname(rj)
        cell, n, rep = os.path.basename(d).rsplit("_", 2)
        params = json.load(open(os.path.join(d, "params.json")))
        stats_path = os.path.join(d, "stats.json")
        if not os.path.exists(stats_path):
            json.dump(aln_stats(read_fasta(os.path.join(d, "true.fasta")), pairs=500), open(stats_path, "w"))
        st = json.load(open(stats_path))
        for m, r in json.load(open(rj)).items():
            if "error" in r:
                continue
            rows.append({"cell": cell, "n": int(n[1:]), "rep": int(rep[1:]), "method": m,
                         "rtt_cv": round(params["rtt_cv"], 4), "mean_p": round(st["mean_p"], 4),
                         "gap_frac": round(st["gap_frac"], 4), "true_len": st["aln_len"],
                         **{k: r.get(k) for k in ("SPFN", "SPFP", "avgErr", "treeFN", "wall", "LenEst")}})
    keys = list(rows[0])
    with open(os.path.join(out, "results.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, keys)
        w.writeheader()
        w.writerows(rows)

    def cellmean(cell, m, k, n=None):
        v = [r[k] for r in rows if r["cell"] == cell and r["method"] == m and r[k] is not None and (n is None or r["n"] == n)]
        return (statistics.mean(v), statistics.stdev(v) / len(v) ** 0.5 if len(v) > 1 else 0, len(v)) if v else None

    # paired MAGUS-vs-PASTA gain per replicate
    def gain(cell, k):
        g = []
        for rep in {r["rep"] for r in rows if r["cell"] == cell}:
            p = [r for r in rows if r["cell"] == cell and r["rep"] == rep]
            mg = [r[k] for r in p if r["method"] == "magus"]
            pa = [r[k] for r in p if r["method"] == "pasta"]
            if mg and pa and pa[0]:
                g.append((pa[0] - mg[0], (pa[0] - mg[0]) / pa[0]))
        return g

    lines = ["| cell | reps | rtt CV | mean p | gap frac | " + " | ".join("{} SP err".format(m) for m in METHODS[:-1])
             + " | MAGUS-PASTA abs (rel) | " + " | ".join("{} tree FN".format(m) for m in METHODS) + " | ranking (SP err) |",
             "|" + "---|" * (6 + 2 * len(METHODS)) + "--|"]
    cells = [c for c in LABEL if any(r["cell"] == c for r in rows)]
    summary = {}
    for c in cells:
        cr = [r for r in rows if r["cell"] == c and r["method"] == "true"]
        reps = len(cr)
        sp = {m: cellmean(c, m, "avgErr") for m in METHODS[:-1]}
        tf = {m: cellmean(c, m, "treeFN") for m in METHODS}
        g = gain(c, "avgErr")
        gt = gain(c, "treeFN")
        rank = " < ".join(m for m, v in sorted(((m, v) for m, v in sp.items() if v), key=lambda x: x[1][0]))
        summary[c] = {"sp": sp, "tree": tf, "gain_abs": [x[0] for x in g], "gain_rel": [x[1] for x in g],
                      "tree_gain_abs": [x[0] for x in gt], "rank": rank}
        fmt = lambda v, s=100: "{:.1f}±{:.1f}".format(v[0] * s, v[1] * s) if v else "–"
        gtxt = "{:+.1f} ({:+.0f}%)".format(100 * statistics.mean(x[0] for x in g), 100 * statistics.mean(x[1] for x in g)) if g else "–"
        lines.append("| {} | {} | {:.2f} | {:.2f} | {:.2f} | {} | {} | {} | {} |".format(
            c, reps, statistics.mean(r["rtt_cv"] for r in cr), statistics.mean(r["mean_p"] for r in cr),
            statistics.mean(r["gap_frac"] for r in cr), " | ".join(fmt(sp[m]) for m in METHODS[:-1]), gtxt,
            " | ".join(fmt(tf[m]) for m in METHODS), rank))
    open(os.path.join(out, "summary.md"), "w").write("\n".join(lines) + "\n")
    json.dump(summary, open(os.path.join(out, "summary.json"), "w"), indent=1)

    # plots: one row per axis, columns = SP error, tree FN, MAGUS-PASTA gain
    fig, axs = plt.subplots(len(AXES), 3, figsize=(14, 3.6 * len(AXES)))
    for i, (axis, cs) in enumerate(AXES.items()):
        cs = [c for c in cs if c in cells]
        x = list(range(len(cs)))
        for j, (k, ms) in enumerate((("avgErr", METHODS[:-1]), ("treeFN", METHODS))):
            ax = axs[i][j]
            for m in ms:
                v = [cellmean(c, m, k) for c in cs]
                if not any(v):
                    continue
                xs = [xx for xx, vv in zip(x, v) if vv]
                ax.errorbar(xs, [100 * vv[0] for vv in v if vv], yerr=[100 * vv[1] for vv in v if vv], marker="o",
                            label=m if m != "true" else "true alignment", color=COLORS[m], capsize=3,
                            linestyle="--" if m == "true" else "-")
            ax.set_xticks(x, [LABEL[c] for c in cs])
            ax.set_xlabel(axis)
            ax.set_ylabel(("SP error (SPFN+SPFP)/2, %" if k == "avgErr" else "FastTree tree FN, %"))
            ax.grid(alpha=0.3)
            if i == 0:
                ax.legend(fontsize=8)
        ax = axs[i][2]
        for k, lab, col in (("gain_rel", "SP error", "#31a354"),):
            vals = [summary[c][k] for c in cs]
            ax.bar(x, [100 * statistics.mean(v) if v else 0 for v in vals], color=col, alpha=0.6, label="relative SP-error reduction")
            for xx, v in zip(x, vals):
                ax.scatter([xx] * len(v), [100 * y for y in v], color="k", s=10, zorder=3)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(x, [LABEL[c] for c in cs])
        ax.set_xlabel(axis)
        ax.set_ylabel("MAGUS gain over PASTA, % of PASTA error")
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "error_vs_factor.png"), dpi=110)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
