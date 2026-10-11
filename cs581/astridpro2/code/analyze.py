"""Tables for REPORT.md from results/*_runs.jsonl.
Usage: python analyze.py RESULTS_DIR > results/tables.md"""
import collections
import glob
import json
import os
import sys

import numpy as np
from scipy.stats import wilcoxon

R = sys.argv[1]
REF = "astrid-pro"
PRIMARY = ["astral-pro3", "astrid-multi", "astrid-disco", "asteroid"]
ORDER = ["astrid-pro", "astrid-pro-r0", "astrid-pro-s", "astrid-multi", "astrid-disco", "asteroid", "astral-pro3",
         "disco-astral", "fastmulrfs", "wqfm-gdl", "duploss2"]
KEYF = ("data", "cond", "rep", "sqln", "ngen")


def load(pattern):
    runs = collections.defaultdict(dict)
    for f in glob.glob(os.path.join(R, pattern)):
        for l in open(f):
            r = json.loads(l)
            if "error" in r:
                continue
            k = tuple(r.get(x) for x in KEYF)
            runs[k][r["method"]] = r
    return runs


def paired(runs, a, b):
    """FN-rate differences a - b over keys where both ran."""
    d = [(v[a]["FNrate"] - v[b]["FNrate"], v[a]["FN"] - v[b]["FN"]) for v in runs.values() if a in v and b in v]
    return np.array([x[0] for x in d]), np.array([x[1] for x in d])


def boot_ci(x, rng=np.random.default_rng(1), B=5000):
    if len(x) == 0:
        return (np.nan, np.nan)
    m = [rng.choice(x, len(x)).mean() for _ in range(B)]
    return np.percentile(m, 2.5), np.percentile(m, 97.5)


def wp(x):
    x = x[x != 0]
    if len(x) < 1:
        return 1.0
    return wilcoxon(x).pvalue


def comp_table(runs, title, methods, holm=False):
    print(f"\n#### {title}\n")
    print("| ASTRID-Pro vs | n | mean diff (FN rate) | 95% CI | W/T/L | p (Wilcoxon) |" + (" Holm p |" if holm else ""))
    print("|---|---|---|---|---|---|" + ("---|" if holm else ""))
    rows = []
    for m in methods:
        d, dfn = paired(runs, REF, m)
        if len(d) == 0:
            continue
        lo, hi = boot_ci(d)
        rows.append([m, len(d), d.mean(), lo, hi, (dfn < 0).sum(), (dfn == 0).sum(), (dfn > 0).sum(), wp(d)])
    if holm:
        ps = [r[-1] for r in rows]
        order = np.argsort(ps)
        adj, run = [0.0] * len(ps), 0.0
        for i, j in enumerate(order):
            run = max(run, min(1.0, (len(ps) - i) * ps[j]))
            adj[j] = run
    for i, r in enumerate(rows):
        s = f"| {r[0]} | {r[1]} | {r[2]:+.4f} | [{r[3]:+.4f}, {r[4]:+.4f}] | {r[5]}/{r[6]}/{r[7]} | {r[8]:.2g} |"
        if holm:
            s += f" {adj[i]:.2g} |"
        print(s)


def mean_table(runs, title, group, methods=ORDER):
    print(f"\n#### {title}\n")
    cols = [m for m in methods if any(m in v for v in runs.values())]
    print("| " + " | ".join(group) + " | n | " + " | ".join(cols) + " |")
    print("|" + "---|" * (len(group) + 1 + len(cols)))
    g = collections.defaultdict(list)
    for k, v in runs.items():
        g[tuple(k[KEYF.index(x)] for x in group)].append(v)
    for gk in sorted(g, key=lambda x: tuple(str(y) for y in x)):
        vs = g[gk]
        cells = []
        best = min((np.mean([v[m]["FNrate"] for v in vs if m in v]) for m in cols if any(m in v for v in vs)))
        for m in cols:
            xs = [v[m]["FNrate"] for v in vs if m in v]
            if not xs:
                cells.append("–")
                continue
            mu = np.mean(xs)
            c = f"{mu:.3f}" + ("" if len(xs) == len(vs) else f" ({len(xs)})")
            cells.append(f"**{c}**" if abs(mu - best) < 1e-12 else c)
        print("| " + " | ".join(str(x) for x in gk) + f" | {len(vs)} | " + " | ".join(cells) + " |")


def time_table(runs, title, group, methods=ORDER):
    print(f"\n#### {title}\n")
    cols = [m for m in methods if any(m in v for v in runs.values())]
    print("| " + " | ".join(group) + " | " + " | ".join(cols) + " |")
    print("|" + "---|" * (len(group) + len(cols)))
    g = collections.defaultdict(list)
    for k, v in runs.items():
        g[tuple(k[KEYF.index(x)] for x in group)].append(v)
    for gk in sorted(g, key=lambda x: tuple((y if isinstance(y, (int, float)) else str(y)) for y in x)):
        vs = g[gk]
        cells = []
        for m in cols:
            xs = [v[m]["sec"] for v in vs if m in v]
            cells.append(f"{np.mean(xs):.2f}" if xs else "–")
        print("| " + " | ".join(str(x) for x in gk) + " | " + " | ".join(cells) + " |")


def failures():
    print("\n### Failed runs (timeouts and crashes; excluded from all paired tables)\n")
    print("| file | method | condition | error | n |")
    print("|---|---|---|---|---|")
    c = collections.Counter()
    for f in sorted(glob.glob(os.path.join(R, "*_runs.jsonl"))):
        ok = set()
        rows = [json.loads(l) for l in open(f)]
        for r in rows:
            if "error" not in r:
                ok.add((tuple(r.get(x) for x in KEYF), r["method"]))
        for r in rows:
            if "error" in r and (tuple(r.get(x) for x in KEYF), r["method"]) not in ok:
                c[(os.path.basename(f), r["method"], r.get("cond"), r["error"])] += 1
    for k, n in sorted(c.items()):
        print("| " + " | ".join(str(x) for x in k) + f" | {n} |")


if __name__ == "__main__":
    fm, di = load("fmrfs_runs.jsonl"), load("disco_runs.jsonl")
    sim = dict(fm)
    sim.update(di)
    others = [m for m in ORDER if m not in (REF,) and m not in PRIMARY]
    print("### Pooled simulated data (FastMulRFS + DISCO, estimated gene trees)")
    comp_table(sim, "Primary comparisons (Holm over 4)", PRIMARY, holm=True)
    comp_table(sim, "Secondary comparisons (no correction)", others)
    for name, runs in (("FastMulRFS data", fm), ("DISCO data", di)):
        print(f"\n### {name}")
        comp_table(runs, f"{name}: ASTRID-Pro vs each method", [m for m in ORDER if m != REF])
    mean_table(fm, "FastMulRFS data: mean FN rate by sequence length and #genes", ["sqln", "ngen"])
    mean_table(fm, "FastMulRFS data: mean FN rate by condition (all sqln, ngen)", ["cond"])
    mean_table(di, "DISCO data: mean FN rate by condition (100 bp, 1000 genes)", ["cond"])
    time_table(fm, "FastMulRFS data: mean runtime (s, 1 thread)", ["ngen"])
    time_table(di, "DISCO data: mean runtime (s, 1 thread)", ["cond"])
    sc = load("scaling_runs.jsonl")
    if sc:
        print("\n### Scaling")
        tx = {k: v for k, v in sc.items() if k[0] == "scale_taxa"}
        gn = {k: v for k, v in sc.items() if k[0] == "scale_genes"}
        time_table(tx, "Runtime (s) vs #species (species_1000 rep 01, ~400-435 genes)", ["cond"])
        mean_table(tx, "FN rate vs #species", ["cond"])
        time_table(gn, "Runtime (s) vs #genes (gtrees_10000_l1, 100 species)", ["ngen"])
        mean_table(gn, "FN rate vs #genes", ["ngen"])
    failures()
    em = load("empirical_runs.jsonl")
    if em:
        print("\n### Empirical")
        mean_table(em, "FN rate vs reference tree", ["cond"])
        time_table(em, "Runtime (s)", ["cond"])
