"""Tables for REPORT.md from results/runs.jsonl (+ results/lnl.jsonl, + published Park et al. trees).

    python summarize.py > ../results/tables.md

Paired comparisons vs base_raxmlng on the replicates both arms finished: mean ΔFN (points), W/T/L with
a 0.5-point tie band (W = arm better by > 0.5), two-sided Wilcoxon signed-rank p (exact; zero
differences dropped), Holm-adjusted p across the pipeline arms of a dataset, CPU ratio (median [min-max]).
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import treeerr  # noqa: E402

R = os.path.join(HERE, "..", "results")
PARK = "/opt/data/park2021/1000M1_HF"
ROSE = "/opt/data/Datasets/ROSE/1000M1"
PUB = {"pub_raxmlng_24h": "RAxML-ng/R{r}/raxmlng-result.raxml.lastTree.TMP",
       "pub_iqtree2": "IQTree/R{r}/iqtree-result.treefile",
       "pub_gtm_iqtree": "GTM/500/iqtree/R{r}/branch_length."}
LABEL = {
    "base_fasttree": "FastTree 2", "base_iqfast": "IQ-TREE 3 --fast", "base_iqtree": "IQ-TREE 3 default",
    "base_raxmlng": "RAxML-NG (1 pars start)",
    "pub_raxmlng_24h": "published RAxML-NG (20 starts, 24 h cap)", "pub_iqtree2": "published IQ-TREE 2",
    "pub_gtm_iqtree": "published GTM (IQ-TREE start)",
}
BASE = "base_raxmlng"
PRIMARY = "place_ft_0.5_fix_rxfast"


def label(a):
    if a in LABEL:
        return LABEL[a]
    p = a.split("_")
    bb = {"ft": "FastTree", "iqf": "IQ-TREE--fast"}[p[1]]
    if p[0] == "constr":
        return "(A) constrained RAxML-NG, %s backbone, τ=%s" % (bb, p[2])
    pol = {"graft": "graft only", "rxfast": "+ RAxML-NG fast polish", "rxfull": "+ RAxML-NG full polish",
           "iqfast": "+ IQ-TREE --fast polish", "ft": "+ FastTree polish"}[p[4]]
    s = "(%s) %s backbone, τ=%s, %s EPA-ng, %s" % ("B" if p[3] == "fix" else "C", bb, p[2],
                                                  "patched" if p[3] == "fix" else "stock", pol)
    return ("**%s** (primary)" % s) if a == PRIMARY else s


def holm(ps):
    idx = sorted(range(len(ps)), key=lambda i: ps[i])
    adj, run = [None] * len(ps), 0.0
    for k, i in enumerate(idx):
        run = max(run, min(1.0, (len(ps) - k) * ps[i]))
        adj[i] = run
    return adj


def load():
    rows = [json.loads(l) for l in open(os.path.join(R, "runs.jsonl"))]
    lnl = {}
    if os.path.exists(os.path.join(R, "lnl.jsonl")):
        for l in open(os.path.join(R, "lnl.jsonl")):
            r = json.loads(l)
            lnl[(r["dataset"], r["rep"], r["arm"])] = r["lnl"]
    by = defaultdict(dict)  # (dataset, aln) -> arm -> rep -> row
    for r in rows:
        r["lnl"] = lnl.get((r["dataset"], r["rep"], r["arm"]))
        by[(r["dataset"], r["aln"])].setdefault(r["arm"], {})[r["rep"]] = r
    pubs = by[("M1HF", "true_align")]
    for a, pat in PUB.items():
        for rep in range(5):
            p = os.path.join(PARK, pat.format(r=rep))
            if os.path.exists(p):
                e = treeerr.error(os.path.join(ROSE, "R%d" % rep, "rose.tt"), p)
                pubs.setdefault(a, {})[rep] = {"fn": e["fn_rate"], "fp": e["fp_rate"], "cpu_s": None,
                                               "wall_s": None, "peak_rss_mb": None, "lnl": None}
    return by


def fmt_p(p):
    return "–" if p is None else ("%.3f" % p if p >= 0.001 else "%.1e" % p)


def table(arms, base_rows):
    out = ["| method | n | mean FN | mean FP | ΔFN vs RAxML-NG (W/T/L) | Wilcoxon p | Holm p | mean CPU min |"
           " CPU ratio vs RAxML-NG, median [range] | mean wall min | max peak RSS MB | mean ΔlnL vs RAxML-NG |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    stats = {}
    for a, reps in arms.items():
        common = sorted(set(reps) & set(base_rows))
        d = [100 * (reps[r]["fn"] - base_rows[r]["fn"]) for r in common]
        p = None
        if a != BASE and len(d) >= 2 and any(abs(x) > 1e-12 for x in d):
            try:
                p = wilcoxon(d).pvalue
            except ValueError:
                p = None
        ratio = [reps[r]["cpu_s"] / base_rows[r]["cpu_s"] for r in common if reps[r].get("cpu_s") is not None]
        dl = [reps[r]["lnl"] - base_rows[r]["lnl"] for r in common
              if reps[r].get("lnl") is not None and base_rows[r].get("lnl") is not None]
        stats[a] = (d, p, ratio, dl)
    pipe = [a for a in arms if a.startswith(("constr", "place")) and stats[a][1] is not None]
    hp = dict(zip(pipe, holm([stats[a][1] for a in pipe]))) if pipe else {}
    for a, reps in arms.items():
        d, p, ratio, dl = stats[a]
        v = list(reps.values())
        cpu = [x["cpu_s"] for x in v if x.get("cpu_s") is not None]
        wall = [x["wall_s"] for x in v if x.get("wall_s") is not None]
        rss = [x["peak_rss_mb"] for x in v if x.get("peak_rss_mb") is not None]
        w = sum(x < -0.5 for x in d)
        l_ = sum(x > 0.5 for x in d)
        out.append("| %s | %d | %.1f%% | %.1f%% | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            label(a), len(v), 100 * np.mean([x["fn"] for x in v]), 100 * np.mean([x["fp"] for x in v]),
            "–" if a == BASE else ("%+.2f (%d/%d/%d)" % (np.mean(d), w, len(d) - w - l_, l_) if d else "–"),
            fmt_p(p), fmt_p(hp.get(a)), "%.1f" % (np.mean(cpu) / 60) if cpu else "–",
            "%.2f [%.2f–%.2f]" % (np.median(ratio), min(ratio), max(ratio)) if ratio and a != BASE else "–",
            "%.1f" % (np.mean(wall) / 60) if wall else "–", "%.0f" % max(rss) if rss else "–",
            "%+.1f (n=%d)" % (np.mean(dl), len(dl)) if dl and a != BASE else "–"))
    return out


ORDER = ["base_fasttree", "base_iqfast", "base_iqtree", "base_raxmlng", "pub_raxmlng_24h", "pub_iqtree2",
         "pub_gtm_iqtree"]


def key(a):
    return (ORDER.index(a) if a in ORDER else 100, a != PRIMARY, a)


def main():
    by = load()
    print("<!-- generated by code/summarize.py -->")
    for (ds, aln), arms in sorted(by.items()):
        arms = {a: arms[a] for a in sorted(arms, key=key)}
        print("\n### %s (%s)\n" % (ds, aln))
        if BASE in arms:
            print("\n".join(table(arms, arms[BASE])))
        else:
            print("| method | n | mean FN | mean FP | per-rep FN | backbone FN | mean CPU min | mean wall min | max peak RSS MB |")
            print("|---|---|---|---|---|---|---|---|---|")
            for a, reps in arms.items():
                v = [reps[r] for r in sorted(reps)]
                print("| %s | %d | %.1f%% | %.1f%% | %s | %s | %.1f | %.1f | %.0f |" % (
                    label(a), len(v), 100 * np.mean([x["fn"] for x in v]), 100 * np.mean([x["fp"] for x in v]),
                    " ".join("%.1f" % (100 * x["fn"]) for x in v),
                    " ".join("%.1f" % (100 * x["backbone_fn"]) for x in v if "backbone_fn" in x) or "–",
                    np.mean([x["cpu_s"] for x in v]) / 60, np.mean([x["wall_s"] for x in v]) / 60,
                    max(x["peak_rss_mb"] for x in v)))
        # per-replicate FN of the main arms
        main_arms = [a for a in (BASE, PRIMARY, "constr_ft_0.5", "place_ft_0.5_fix_graft", "base_fasttree")
                     if a in arms]
        reps = sorted(set().union(*[arms[a].keys() for a in main_arms])) if main_arms else []
        if reps and BASE in arms:
            print("\nPer-replicate FN %% (backbone FN of the FastTree τ=0.5 backbone in the last row)\n")
            print("| method | " + " | ".join("R%d" % r for r in reps) + " |")
            print("|---|" + "---|" * len(reps))
            for a in main_arms:
                print("| %s | %s |" % (label(a), " | ".join(
                    "%.1f" % (100 * arms[a][r]["fn"]) if r in arms[a] else "–" for r in reps)))
            src = arms.get(PRIMARY) or arms.get("place_ft_0.5_fix_graft") or {}
            print("| backbone (FastTree, τ=0.5) | %s |" % " | ".join(
                "%.1f" % (100 * src[r]["backbone_fn"]) if r in src else "–" for r in reps))
        # placement identity, patched vs stock
        for pol in ("graft", "ft"):
            f, s = arms.get("place_ft_0.5_fix_" + pol, {}), arms.get("place_ft_0.5_stock_" + pol, {})
            c = sorted(set(f) & set(s))
            if c:
                print("\nPatched vs stock EPA-ng, FastTree backbone τ=0.5, %s: identical FN on %d/%d reps; "
                      "mean FN patched %.2f%% vs stock %.2f%%; per rep (patched/stock): %s" % (
                          pol, sum(abs(f[r]["fn"] - s[r]["fn"]) < 1e-12 for r in c), len(c),
                          100 * np.mean([f[r]["fn"] for r in c]), 100 * np.mean([s[r]["fn"] for r in c]),
                          ", ".join("R%d %.1f/%.1f" % (r, 100 * f[r]["fn"], 100 * s[r]["fn"]) for r in c)))


if __name__ == "__main__":
    main()
