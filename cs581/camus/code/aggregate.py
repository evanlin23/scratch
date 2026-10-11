"""Aggregate run_rep.py outputs: reproduction vs published CAMUS networks, and paired
comparisons of each variant against default CAMUS (Wilcoxon signed-rank, W/T/L).

usage: python3 aggregate.py RESDIR OUT_MD [OUT_CSV]
Error = (FN + FP)/2 of softwired clusters (PhyloNet CmpNets -m cluster), 1-reticulation network
unless noted. Tie band: |diff| < 1e-9 (cluster errors are discrete fractions).
"""
import csv
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import netmetric as nm  # noqa: E402

DATA = "/opt/data/camus/sim/camus-dataset"
PUB = "/opt/data/camus/inf/inferred-networks"
TRAIN = {("n25", "g_500")}  # the paper tuned t on 26-taxon FastTree data; we tune on the same
TIE = 1e-9
GTAG = {"g_500": "fasttree", "iqtree_500": "iqtree", "g_true": "true-gt"}


def load(resdir):
    recs = []
    for f in sorted(glob.glob(f"{resdir}/*.jsonl")):
        for l in open(f):
            r = json.loads(l)
            for k in ("k1", "ktrue"):
                r[f"err_{k}"] = (r[f"fn_{k}"] + r[f"fp_{k}"]) / 2
            recs.append(r)
    return recs


def pub_nets(cond, rep, gt, f="05"):
    """Published CAMUS networks (list indexed by k-1)."""
    base = f"{PUB}/{cond}/{rep}/camus_f{f}_astral_{GTAG[gt]}"
    if os.path.exists(base + ".csv"):
        rows = list(csv.reader(open(base + ".csv")))[1:]
        return [r[2] for r in rows if int(r[0]) > 0]
    if os.path.exists(base + ".nwk"):
        return [l.strip() for l in open(base + ".nwk") if l.strip()]
    return None


def fmt(x, d=3):
    return f"{x:.{d}f}"


def paired(recs, var, key="err_k1", ref="default"):
    by = defaultdict(dict)
    for r in recs:
        by[(r["cond"], r["gt"], r["rep"])][r["variant"]] = r
    rows = []
    groups = defaultdict(list)
    for (c, g, rp), d in by.items():
        if var in d and ref in d:
            groups[(c, g)].append((d[var][key] - d[ref][key], d[ref][key], d[var][key]))
    for (c, g), L in sorted(groups.items()):
        rows.append(((c, g), L))
    return rows


def stats(L):
    d = np.array([x[0] for x in L])
    w = int((d < -TIE).sum())
    l_ = int((d > TIE).sum())
    t = len(d) - w - l_
    nz = d[np.abs(d) > TIE]
    p = wilcoxon(nz).pvalue if len(nz) >= 1 and np.any(nz) else 1.0
    return len(d), np.mean([x[1] for x in L]), np.mean([x[2] for x in L]), d.mean(), w, t, l_, p


def main():
    resdir, out_md = sys.argv[1:3]
    recs = load(resdir)
    lines = ["# CAMUS pilot: aggregated results", "",
             "Error = (FN+FP)/2 of softwired clusters (equal to PhyloNet `CmpNets -m cluster`).",
             "Paired by replicate; diff = variant - default (negative = better).",
             "W/T/L = variant better / tie (|diff| < 1e-9) / worse. p = two-sided Wilcoxon signed-rank "
             "(zero differences dropped). Training condition (used to pick z settings): n25 FastTree (reps 20-39).",
             ""]

    # --- reproduction
    lines += ["## Reproduction against the published CAMUS networks (V2 data, t = 0.5)", "",
              "| condition | reps | our default = published (k=1 network identical) | published FN/FP (k=1) | ours FN/FP (k=1) |",
              "|---|---|---|---|---|"]
    by = defaultdict(list)
    for r in recs:
        if r["variant"] != "default":
            continue
        pn = pub_nets(r["cond"], r["rep"], r["gt"])
        if not pn:
            continue
        true = nm.softwired_clusters(nm.parse_enewick(open(f"{DATA}/{r['cond']}/{r['rep']}/true_net.nwk").readline())[0])
        pcl = nm.softwired_clusters(nm.parse_enewick(pn[0])[0])
        ocl = nm.softwired_clusters(nm.parse_enewick(r["net_k1"])[0])
        pfn, pfp = nm.cluster_error(true, pcl)
        by[(r["cond"], r["gt"])].append((pcl == ocl, pfn, pfp, r["fn_k1"], r["fp_k1"]))
    for (c, g), L in sorted(by.items()):
        a = np.array([x[1:] for x in L])
        lines.append(f"| {c} {GTAG.get(g, g)} | {len(L)} | {sum(x[0] for x in L)}/{len(L)} | "
                     f"{fmt(a[:, 0].mean())} / {fmt(a[:, 1].mean())} | {fmt(a[:, 2].mean())} / {fmt(a[:, 3].mean())} |")
    lines.append("")

    # --- paired comparisons
    variants = sorted({r["variant"] for r in recs} - {"default"})
    for key, title in (("err_k1", "1-reticulation network (paper's protocol)"),
                       ("err_ktrue", "network with the true number of reticulations")):
        lines += [f"## Variants vs default: {title}", "",
                  "| variant | condition | n | default err | variant err | mean diff | W/T/L | p |",
                  "|---|---|---|---|---|---|---|---|"]
        for v in variants:
            for (c, g), L in paired(recs, v, key):
                n, a, b, md, w, t, l_, p = stats(L)
                tag = " (train)" if (c, g) in TRAIN else ""
                lines.append(f"| {v} | {c} {GTAG.get(g, g)}{tag} | {n} | {fmt(a)} | {fmt(b)} | {md:+.4f} | {w}/{t}/{l_} | {p:.3g} |")
            # pooled held-out
            L = [x for (cg, LL) in paired(recs, v, key) if cg not in TRAIN for x in LL]
            if L:
                n, a, b, md, w, t, l_, p = stats(L)
                lines.append(f"| {v} | **held-out pooled** | {n} | {fmt(a)} | {fmt(b)} | {md:+.4f} | {w}/{t}/{l_} | {p:.3g} |")
        lines.append("")

    # --- base tree error
    lines += ["## Base-tree error (softwired = tree clusters vs true network clusters) and runtime", "",
              "| condition | variant | base FN | base FP | k=1 FN | k=1 FP | mean CAMUS time (s, 1 thread) |", "|---|---|---|---|---|---|---|"]
    agg = defaultdict(list)
    for r in recs:
        agg[(r["cond"], r["gt"], r["variant"])].append(r)
    for (c, g, v), L in sorted(agg.items()):
        if v not in ("default", "tqmc", "wastral", "true_major", "swap"):
            continue
        m = lambda k: np.mean([x[k] for x in L])  # noqa: E731
        lines.append(f"| {c} {GTAG.get(g, g)} | {v} | {fmt(m('base_fn'))} | {fmt(m('base_fp'))} | {fmt(m('fn_k1'))} | {fmt(m('fp_k1'))} | {m('time_s'):.1f} |")
    lines.append("")

    # --- score-based selection among base trees (no access to the truth)
    lines += ["## Selecting the base tree by the CAMUS objective (no truth used)", "",
              "For each replicate, among the candidate 1-reticulation networks (default, tqmc, wastral, swap), pick the "
              "one maximizing the number of filtered (t=0.5) gene-tree quartets displayed. Ties -> default.", "",
              "| condition | n | default err | selected err | mean diff | W/T/L | p | oracle-best err |", "|---|---|---|---|---|---|---|---|"]
    byrep = defaultdict(dict)
    for r in recs:
        byrep[(r["cond"], r["gt"], r["rep"])][r["variant"]] = r
    groups = defaultdict(list)
    cands = ["default", "tqmc", "wastral", "swap"]
    for (c, g, rp), d in byrep.items():
        if not all(x in d and d[x].get("score_k1") is not None for x in cands):
            continue
        best = max(cands, key=lambda x: (d[x]["score_k1"], x == "default"))
        groups[(c, g)].append((d[best]["err_k1"] - d["default"]["err_k1"], d["default"]["err_k1"], d[best]["err_k1"],
                               min(d[x]["err_k1"] for x in cands)))
    for (c, g), L in sorted(groups.items()):
        n, a, b, md, w, t, l_, p = stats(L)
        lines.append(f"| {c} {GTAG.get(g, g)} | {n} | {fmt(a)} | {fmt(b)} | {md:+.4f} | {w}/{t}/{l_} | {p:.3g} | {fmt(np.mean([x[3] for x in L]))} |")
    lines.append("")

    # --- threshold oracle: best fixed t per condition and per-replicate oracle t
    lines += ["## Threshold t: per-condition best fixed t and per-replicate oracle (upper bound for adaptive t)", "",
              "| condition | t=0 | t=0.2 | t=0.3 | t=0.5 | t=0.7 | t=0.8 | per-rep oracle t |", "|---|---|---|---|---|---|---|---|"]
    tv = ["t0", "t0.2", "t0.3", "default", "t0.7", "t0.8"]
    tg = defaultdict(list)
    for (c, g, rp), d in byrep.items():
        if all(x in d for x in tv):
            tg[(c, g)].append([d[x]["err_k1"] for x in tv])
    for (c, g), L in sorted(tg.items()):
        a = np.array(L)
        lines.append(f"| {c} {GTAG.get(g, g)} | " + " | ".join(fmt(x) for x in a.mean(0)) + f" | {fmt(a.min(1).mean())} |")
    lines.append("")
    open(out_md, "w").write("\n".join(lines) + "\n")
    if len(sys.argv) > 3:
        keys = ["cond", "gt", "rep", "variant", "base", "t", "z", "zmode", "r_true", "kmax", "time_s",
                "n_nontree_quartets", "base_fn", "base_fp", "fn_k1", "fp_k1", "fn_ktrue", "fp_ktrue", "score_k1"]
        with open(sys.argv[3], "w") as f:
            w = csv.writer(f)
            w.writerow(keys)
            for r in recs:
                w.writerow([r.get(k) for k in keys])


if __name__ == "__main__":
    main()
