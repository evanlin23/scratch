"""Paired comparison of every variant against MAGUS (same replicate, same machine) + Pareto plot.

    python analyze.py   ->  results/summary.md, results/pareto.png
"""

import glob
import json
import os
from collections import defaultdict

import numpy as np
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
PUB = os.path.join(HERE, "..", "..", "experiments", "validation", "published_scores.jsonl")
TIE = 0.05


def load():
    rows = {}
    for path in sorted(glob.glob(os.path.join(RES, "*.jsonl"))):
        for line in open(path):
            r = json.loads(line)
            rows.setdefault(r["dataset"], {}).update(r)  # later rows add variants
    return rows


def published():
    """(dataset, rep) -> method -> (err %, seconds), from the authors' published alignments/logs."""
    pub = defaultdict(dict)
    for line in open(PUB):
        r = json.loads(line)
        if r.get("published_seconds"):
            key = r["rep"].replace("RV100_", "") + "_R0" if r["dataset"] == "balibase" else "{}_{}".format(r["dataset"], r["rep"])
            pub[key][r["method"]] = (100 * r["avgErr"], r["published_seconds"])
    return pub


def main():
    rows = load()
    variants = []
    for r in rows.values():
        for k, v in r.items():
            if isinstance(v, dict) and "avgErr" in v and k not in variants:
                variants.append(k)
    lines = ["# Fast MAGUS pilot: paired comparison with MAGUS (paper flags), same machine", "",
             "Error = (SPFN+SPFP)/2 in %. Delta < 0: lower error than MAGUS. W/T/L: replicates where the "
             "variant is better / within {} pts / worse. Runtime ratio = MAGUS wall / variant wall "
             "(> 1: faster), mean over replicates.".format(TIE), "",
             "## Per replicate (error %, wall s)", ""]
    ds = sorted(rows)
    lines.append("| variant | " + " | ".join(ds) + " |")
    lines.append("|---|" + "---|" * len(ds))
    for v in variants:
        cells = []
        for d in ds:
            x = rows[d].get(v)
            cells.append("{:.2f} ({:.0f})".format(100 * x["avgErr"], x["wall"]) if x and "avgErr" in x else "-")
        lines.append("| {} | {} |".format(v, " | ".join(cells)))
    lines += ["", "## Paired summary", "",
              "| variant | n | mean err | mean delta (pts) | W/T/L | Wilcoxon p | mean wall s | speedup vs MAGUS (mean of ratios) | CPU ratio |",
              "|---|---|---|---|---|---|---|---|---|"]
    points = {}
    for v in variants:
        pairs = [(rows[d]["magus"], rows[d][v]) for d in ds
                 if "avgErr" in rows[d].get("magus", {}) and "avgErr" in rows[d].get(v, {})]
        if not pairs:
            continue
        delta = np.array([100 * (b["avgErr"] - a["avgErr"]) for a, b in pairs])
        speed = np.array([a["wall"] / b["wall"] for a, b in pairs])
        cpu = np.array([a["cpu"] / b["cpu"] for a, b in pairs])
        w, t, l = (delta < -TIE).sum(), (abs(delta) <= TIE).sum(), (delta > TIE).sum()
        p = wilcoxon(delta).pvalue if len(delta) >= 2 and np.any(delta != 0) else float("nan")
        lines.append("| {} | {} | {:.2f} | {:+.2f} | {}/{}/{} | {} | {:.0f} | {:.2f}x | {:.2f}x |".format(
            v, len(pairs), np.mean([100 * b["avgErr"] for _, b in pairs]), delta.mean(), w, t, l,
            "{:.3g}".format(p) if p == p else "-", np.mean([b["wall"] for _, b in pairs]), speed.mean(), cpu.mean()))
        points[v] = (1 / speed.mean(), delta.mean(), len(pairs))

    # the same comparison against the authors' published MAGUS(Fast) alignment of each replicate
    # (a different random backbone draw: shows how much of a delta is backbone-draw noise)
    pub = published()
    lines += ["", "## Against the published MAGUS(Fast) alignment of the same replicate", "",
              "| variant | n | mean delta vs published MAGUS(Fast) (pts) | W/T/L |", "|---|---|---|---|"]
    for v in variants:
        dl = [100 * rows[d][v]["avgErr"] - pub[d]["MAGUS(Fast)"][0] for d in ds
              if "avgErr" in rows[d].get(v, {}) and "MAGUS(Fast)" in pub.get(d, {})]
        if dl:
            dl = np.array(dl)
            lines.append("| {} | {} | {:+.2f} | {}/{}/{} |".format(v, len(dl), dl.mean(), (dl < -TIE).sum(),
                                                               (abs(dl) <= TIE).sum(), (dl > TIE).sum()))

    # stage breakdown of MAGUS runs
    lines += ["", "## MAGUS (paper flags) stage breakdown on this machine", "",
              "| replicate | wall s | guide tree s (%) | subsets+backbones until graph built s (%) | cluster+trace s (%) |",
              "|---|---|---|---|---|"]
    for d in ds:
        m = rows[d].get("magus", {})
        st = m.get("stages", {})
        if "total" not in st:
            continue
        tot = st["total"]
        bb = st.get("graph", 0)
        rest = tot - st.get("decomp", 0) - bb
        lines.append("| {} | {:.0f} | {:.0f} ({:.0f}%) | {:.0f} ({:.0f}%) | {:.0f} ({:.0f}%) |".format(
            d, m["wall"], st["decomp"], 100 * st["decomp"] / tot, bb, 100 * bb / tot, rest, 100 * rest / tot))

    # PASTA reference (published logs, same replicates): ratio of runtimes, difference of errors
    ref = []
    for d in ds:
        p = pub.get(d, {})
        if "PASTA(3)" in p and "MAGUS(Fast)" in p:
            ref.append((p["PASTA(3)"][1] / p["MAGUS(Fast)"][1], p["PASTA(3)"][0] - p["MAGUS(Fast)"][0]))
    if ref:
        rx, ry = np.mean([r[0] for r in ref]), np.mean([r[1] for r in ref])
        lines += ["", "PASTA(3) (published alignments and logs, same {} replicates): {:.2f}x MAGUS(Fast)'s runtime, "
                  "{:+.2f} pts error vs MAGUS(Fast).".format(len(ref), rx, ry)]
    open(os.path.join(RES, "summary.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.5, 5))
    xs = {v: p[0] for v, p in points.items()}
    front = []
    for v, (x, y, n) in sorted(((v, p) for v, p in points.items() if p[2] == max(q[2] for q in points.values())),
                               key=lambda kv: kv[1][0]):
        if not front or y < front[-1][1]:
            front.append((x, y))
    ax.plot([f[0] for f in front], [f[1] for f in front], "-", color="#999", lw=1, zorder=1, label="Pareto front (variants run on all replicates)")
    for v, (x, y, n) in points.items():
        if y > 8:
            continue  # alignment-free guide trees: far off the scale, listed in the tables
        c = "#d62728" if v == "magus" else ("#1f77b4" if "+ss" in v else "#7f7f7f")
        ax.scatter(x, y, color=c, zorder=2)
        ax.annotate("{} (n={})".format(v, n), (x, y), textcoords="offset points", xytext=(4, 4), fontsize=8)
    if ref:
        ax.scatter(rx, ry, marker="s", color="#2ca02c", zorder=2)
        ax.annotate("PASTA(3), published", (rx, ry), textcoords="offset points", xytext=(4, 4), fontsize=8)
    ax.axhline(0, color="#ccc", lw=0.8)
    ax.axvline(1, color="#ccc", lw=0.8)
    ax.axvline(0.5, color="#ccc", lw=0.8, ls="--")
    ax.set_ylim(-1.5, 7)
    ax.set_xscale("log")
    ax.set_xlabel("wall-clock relative to MAGUS (log; dashed = 2x faster)")
    ax.set_ylabel("error minus MAGUS error (pts)")
    ax.set_title("Accuracy vs runtime, paired by replicate (n = {})".format(max(p[2] for p in points.values())))
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "pareto.png"), dpi=130)


if __name__ == "__main__":
    main()
