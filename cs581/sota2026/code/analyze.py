"""Tables + Pareto plots for the sota2026 benchmark.

    python3 analyze.py   -> results/full_table.md, results/subset_table.md, results/pareto_*.png

Reference points for MAGUS(Fast) / PASTA(3): FastSP on the authors' published alignments of the
same replicate (cs581/experiments/validation/published_scores.jsonl). MAGUS runtime on this
machine: our reruns (tool "magus" in full.jsonl) and earlier measured sessions on the same
4-core machine type (cs581-timing-a/b branches, quoted in REPORT.md).
"""

import collections
import json
import os
import statistics

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
PUB = os.path.join(HERE, "..", "..", "experiments", "validation", "published_scores.jsonl")
TIE = 0.05  # points


def dtype(ds):
    if ds.startswith("1000"):
        return "ROSE"
    if ds.startswith("RNASim"):
        return "RNASim"
    if ds.startswith("BBA"):
        return "BAliBASE"
    return "16S"


def published():
    pub = {}
    for r in map(json.loads, open(PUB)):
        if r["method"] not in ("MAGUS(Fast)", "PASTA(3)") or "avgErr" not in r:
            continue
        ds = r["rep"].replace("RV100_", "") + "_R0" if r["dataset"] == "balibase" else r["dataset"] + "_" + r["rep"]
        pub[(ds, r["method"])] = r
    return pub


def paired(rows, a, b):
    """rows: {(ds, method): err%}; paired Δ = a - b over datasets with both."""
    ds = sorted({d for d, m in rows if m == a} & {d for d, m in rows if m == b})
    diffs = [rows[(d, a)] - rows[(d, b)] for d in ds]
    if not diffs:
        return None
    w = sum(x < -TIE for x in diffs)
    l = sum(x > TIE for x in diffs)
    t = len(diffs) - w - l
    try:
        p = wilcoxon(diffs).pvalue if len(diffs) >= 5 and any(diffs) else float("nan")
    except ValueError:
        p = float("nan")
    return {"n": len(diffs), "mean": statistics.mean(diffs), "W": w, "T": t, "L": l, "p": p}


def full():
    rows = [json.loads(l) for l in open(os.path.join(RES, "full.jsonl"))]
    pub = published()
    err, wall, status = {}, {}, {}
    for r in rows:
        key = (r["dataset"], r["tool"])
        status[key] = r["status"]
        wall[key] = r["wall"]
        if r["status"] == "ok":
            err[key] = 100 * r["avgErr"]
    datasets = list(dict.fromkeys(r["dataset"] for r in rows))
    alltools = list(dict.fromkeys(r["tool"] for r in rows))
    tools = [t for t in alltools if "truetree" not in t]
    diag = [t for t in alltools if "truetree" in t]
    for ds in datasets:
        for m, name in (("MAGUS(Fast)", "MAGUS(pub)"), ("PASTA(3)", "PASTA(pub)")):
            if (ds, m) in pub:
                err[(ds, name)] = 100 * pub[(ds, m)]["avgErr"]
                wall[(ds, name)] = pub[(ds, m)].get("published_seconds")
    runs = os.path.join(HERE, "..", "..", "experiments", "runs")
    for ds in datasets:  # earlier MAGUS(Fast) reruns, paper flags, 4 cores, same machine type (prep.json)
        pj = os.path.join(runs, ds, "prep.json")
        if os.path.exists(pj):
            r = json.load(open(pj))
            err[(ds, "MAGUS(4c)")] = 100 * r["avgErr"]
            wall[(ds, "MAGUS(4c)")] = r["seconds"]
    cols = ["MAGUS(pub)", "MAGUS(4c)", "PASTA(pub)"] + tools
    out = ["## Full datasets: average error (SPFN+SPFP)/2, %, and wall-clock seconds (4 threads)\n",
           "MAGUS(pub)/PASTA(pub) = FastSP on the authors' published alignment of the same replicate; their "
           "seconds are the paper's own timing (different hardware). MAGUS(4c) = earlier rerun of MAGUS(Fast) with the "
           "paper's flags on this 4-core machine type (cs581/experiments/runs/<rep>/prep.json; random backbones, so "
           "not identical to the published alignment).\n",
           "| dataset | " + " | ".join(cols) + " |", "|---" * (len(cols) + 1) + "|"]
    for ds in datasets:
        cells = []
        for c in cols:
            if (ds, c) in err:
                s = "{:.1f}".format(err[(ds, c)])
                if wall.get((ds, c)) is not None:
                    s += " ({:.0f}s)".format(wall[(ds, c)])
            elif (ds, c) in status:
                s = status[(ds, c)] + " ({:.0f}s)".format(wall[(ds, c)])
            else:
                s = ""
            cells.append(s)
        out.append("| {} | {} |".format(ds, " | ".join(cells)))
    out += ["", "## Paired comparison vs published MAGUS(Fast), same replicate (Δ = tool − MAGUS, points; "
                "W/T/L = tool better/tie(±{})/worse; Wilcoxon signed-rank)\n".format(TIE),
            "| data | tool | n | mean Δ | W/T/L | p |", "|---|---|---|---|---|---|"]
    for group in ("ROSE", "RNASim", "BAliBASE", "16S", "all"):
        sub = {k: v for k, v in err.items() if group == "all" or dtype(k[0]) == group}
        for t in tools + ["MAGUS(4c)", "PASTA(pub)"]:
            s = paired(sub, t, "MAGUS(pub)")
            if s:
                out.append("| {} | {} | {} | {:+.1f} | {}/{}/{} | {} |".format(
                    group, t, s["n"], s["mean"], s["W"], s["T"], s["L"],
                    "{:.2g}".format(s["p"]) if s["p"] == s["p"] else "–"))
    if diag:
        out += ["", "## Diagnostic: same aligner, true tree as guide tree (error %, seconds)\n",
                "| dataset | MAGUS(pub) | famsa | famsa-truetree | twilight-1 | twilight | twilight-truetree |",
                "|---|---|---|---|---|---|---|"]
        for ds in datasets:
            if not any((ds, t) in status for t in diag):
                continue
            cells = []
            for c in ("MAGUS(pub)", "famsa", "famsa-truetree", "twilight-1", "twilight", "twilight-truetree"):
                cells.append("{:.1f} ({:.0f}s)".format(err[(ds, c)], wall[(ds, c)]) if (ds, c) in err and wall.get((ds, c)) is not None
                             else (status.get((ds, c), "")))
            out.append("| {} | {} |".format(ds, " | ".join(cells)))
    open(os.path.join(RES, "full_table.md"), "w").write("\n".join(out) + "\n")

    # Pareto plots per data type: mean error vs mean wall-clock over the datasets of that type
    for group in ("ROSE", "RNASim", "BAliBASE", "16S"):
        dss = [d for d in datasets if dtype(d) == group]
        pts = {}
        for c in cols:
            have = [d for d in dss if (d, c) in err]
            if not have:
                continue
            ws = [wall.get((d, c)) for d in have]
            if any(w is None for w in ws):
                continue
            label = c if len(have) == len(dss) else "{} (n={})".format(c, len(have))
            pts[label] = (statistics.mean(ws), statistics.mean(err[(d, c)] for d in have))
        if not pts:
            continue
        fig, ax = plt.subplots(figsize=(6.4, 4.4))
        front, best = [], float("inf")
        for c, (x, y) in sorted(pts.items(), key=lambda kv: kv[1][0]):
            if "n=" in c:  # partial coverage: not comparable enough for the front
                continue
            if y < best:
                front.append((x, y))
                best = y
        ax.plot(*zip(*front), color="#999", lw=1, ls="--", zorder=1)
        for k, (c, (x, y)) in enumerate(sorted(pts.items(), key=lambda kv: kv[1])):
            pub_pt = "MAGUS" in c or "PASTA" in c
            ax.scatter(x, y, s=40, marker="s" if pub_pt else "o", color="#c44" if "MAGUS" in c else "#36a", zorder=2)
            ax.annotate(c, (x, y), textcoords="offset points", xytext=(4, 4 if k % 2 == 0 else -11), fontsize=8)
        ax.set_xscale("log")
        ax.set_xlabel("mean wall-clock, s (log; 4 threads here; '(pub)' = paper's own hardware)")
        ax.set_ylabel("mean error (SPFN+SPFP)/2, %")
        ax.set_title("{} ({} dataset{})".format(group, len(dss), "s" if len(dss) > 1 else ""))
        fig.tight_layout()
        fig.savefig(os.path.join(RES, "pareto_{}.png".format(group)), dpi=110)
        plt.close(fig)


def subsets():
    path = os.path.join(RES, "subsets.jsonl")
    if not os.path.exists(path):
        return
    rows = [json.loads(l) for l in open(path)]
    err = {(r["rep"], r["id"], r["tool"]): 100 * r["avgErr"] for r in rows if r["status"] == "ok"}
    kinds = {(r["rep"], r["id"]): r["kind"] for r in rows}
    wall = collections.defaultdict(list)
    for r in rows:
        if r["status"] == "ok" and "wall" in r:
            wall[(r["kind"], r["tool"])].append(r["wall"])
    tools = list(dict.fromkeys(r["tool"] for r in rows))
    out = ["## Subset level: Δ error vs MAFFT L-INS-i (MAGUS's command, rerun here), points; "
           "negative = more accurate than L-INS-i\n",
           "| kind | tool | n | mean err | mean Δ | W/T/L | p | mean s (1 core) |", "|---|---|---|---|---|---|---|---|"]
    base = "mafft-linsi-magus"
    for kind in ("magus40", "rand40", "bb200", "rand200"):
        for t in tools:
            keys = [k for k in kinds if kinds[k] == kind and k + (t,) in err and k + (base,) in err]
            if not keys:
                continue
            diffs = [err[k + (t,)] - err[k + (base,)] for k in keys]
            w, l = sum(x < -TIE for x in diffs), sum(x > TIE for x in diffs)
            try:
                p = wilcoxon(diffs).pvalue if len(diffs) >= 5 and any(diffs) else float("nan")
            except ValueError:
                p = float("nan")
            ws = wall.get((kind, t))
            out.append("| {} | {} | {} | {:.1f} | {:+.2f} | {}/{}/{} | {} | {} |".format(
                kind, t, len(keys), statistics.mean(err[k + (t,)] for k in keys), statistics.mean(diffs), w,
                len(diffs) - w - l, l, "{:.2g}".format(p) if p == p else "–",
                "{:.0f}".format(statistics.mean(ws)) if ws else "–"))
    # per data type for magus40
    out += ["", "### magus40 by data type (mean Δ vs L-INS-i, points)\n", "| tool | " + " | ".join(
        ("ROSE", "RNASim", "16S", "BAliBASE")) + " |", "|---|---|---|---|---|"]
    for t in tools:
        cells = []
        for g in ("ROSE", "RNASim", "16S", "BAliBASE"):
            keys = [k for k in kinds if kinds[k] in ("magus40",) and dtype(k[0]) == g and k + (t,) in err
                    and k + (base,) in err]
            cells.append("{:+.2f} (n={})".format(statistics.mean(err[k + (t,)] - err[k + (base,)] for k in keys),
                                                 len(keys)) if keys else "")
        out.append("| {} | {} |".format(t, " | ".join(cells)))
    open(os.path.join(RES, "subset_table.md"), "w").write("\n".join(out) + "\n")


if __name__ == "__main__":
    full()
    subsets()


def merge_pilot():
    path = os.path.join(RES, "merge_pilot.jsonl")
    if not os.path.exists(path):
        return
    rows = [r for r in map(json.loads, open(path)) if r["status"] == "ok"]
    for r in rows:
        if r["backbones"] == "tool":
            r["tool"] += " (subsets+backbones)"
    ctl = {r["rep"]: r for r in rows if r["tool"] == "cached-linsi"}
    out = ["## MAGUS with a different base method (merge pilot)\n",
           "Same MAGUS decomposition (25 subsets) and the same 10 MAFFT L-INS-i backbones as MAGUS's own run; "
           "only the subset aligner changes, then MAGUS's GCM merge (paper flags). Control = MAGUS's own "
           "L-INS-i subsets (reproduces MAGUS(4c)). Errors in %, Δ in points (negative = better than MAGUS).\n",
           "| rep | base method | subset err (mean of 25) | MAGUS final err | Δ final vs control | subset align s (1 core) |",
           "|---|---|---|---|---|---|"]
    deltas = collections.defaultdict(list)
    for r in sorted(rows, key=lambda r: (r["rep"], r["tool"] != "cached-linsi", r["tool"])):
        c = ctl.get(r["rep"])
        d = 100 * (r["avgErr"] - c["avgErr"]) if c else float("nan")
        if r["tool"] != "cached-linsi" and c:
            deltas[(dtype(r["rep"]), r["tool"])].append(d)
        out.append("| {} | {} | {:.2f} | {:.2f} | {} | {} |".format(
            r["rep"], "L-INS-i (MAGUS's own)" if r["tool"] == "cached-linsi" else r["tool"], 100 * r["subset_err_mean"],
            100 * r["avgErr"], "–" if r["tool"] == "cached-linsi" else "{:+.2f}".format(d),
            "–" if r["tool"] == "cached-linsi" else "{:.0f}".format(r["subset_wall"])))
    out += ["", "| data | base method | n | mean Δ | better/tie/worse |", "|---|---|---|---|---|"]
    for (g, t), ds in sorted(deltas.items()):
        w, l = sum(x < -TIE for x in ds), sum(x > TIE for x in ds)
        out.append("| {} | {} | {} | {:+.2f} | {}/{}/{} |".format(g, t, len(ds), statistics.mean(ds), w,
                                                                   len(ds) - w - l, l))
    open(os.path.join(RES, "merge_table.md"), "w").write("\n".join(out) + "\n")


merge_pilot()
