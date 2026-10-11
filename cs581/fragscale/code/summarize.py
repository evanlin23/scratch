"""Tables and anytime curves for the fragscale pilot (implements PREREG.md).

    MLDATA=/opt/data/fscache python summarize.py [--plot]

Reads results/runs.jsonl and the anytime snapshot indexes under $MLDATA/<ds>/R<rep>/trees/anytime/.
Writes results/tables.md, results/anytime.tsv (and results/anytime.png with --plot).
"""
import collections
import json
import os
import statistics
import sys

from scipy.stats import wilcoxon

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import anytime  # noqa: E402
import runtrees as rt  # noqa: E402
import treeerr  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
PRIMARY = "place_ft_0.5_fix_rxfast"
FRACS = [0.05, 0.1, 0.15, 0.2, 0.25, 1 / 3, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0]
_cache = {}


def fn_of(true, tree):
    if tree is None:
        return None
    k = (true, tree)
    if k not in _cache:
        _cache[k] = 100 * treeerr.error(true, tree)["fn_rate"]
    return _cache[k]


def load():
    rows = {}
    for line in open(os.path.join(RES, "runs.jsonl")):
        r = json.loads(line)
        rows[(r["dataset"], r["rep"], r["aln"], r["arm"])] = r
    return rows


def snapdir(ds, rep, aln, m):
    return os.path.join(rt.MLDATA, ds, "R%d" % rep, "trees", "anytime", "%s.%s" % (aln, m))


def anytime_fn(ds, rep, aln, m, T):
    d = snapdir(ds, rep, aln, m)
    return fn_of(os.path.join(rt.MLDATA, ds, "R%d" % rep, "true_tree.tre"), anytime.at(d, T))


def paired(a, b, band=0.5):
    """a, b: dict rep -> FN (pp). Returns n, mean diff (a-b), W/T/L (a better = W), p."""
    reps = sorted(set(a) & set(b))
    reps = [r for r in reps if a[r] is not None and b[r] is not None]
    d = [a[r] - b[r] for r in reps]
    if not d:
        return None
    w = sum(x < -band for x in d)
    l_ = sum(x > band for x in d)
    t = len(d) - w - l_
    p = None
    if len(d) >= 2 and any(x != 0 for x in d):
        p = wilcoxon(d, zero_method="wilcox", alternative="two-sided", method="exact").pvalue
    return {"n": len(d), "mean": statistics.mean(d), "W": w, "T": t, "L": l_, "p": p, "d": dict(zip(reps, d))}


def fmt(pr):
    if pr is None:
        return "– | – | –"
    p = "–" if pr["p"] is None else "%.3f" % pr["p"]
    return "%+.2f | %d/%d/%d | %s" % (pr["mean"], pr["W"], pr["T"], pr["L"], p)


def section_equal_time(rows, ds, aln, out, curves):
    reps = sorted(r for (d, r, a, m) in rows if d == ds and a == aln and m == "base_raxmlng")
    if not reps:
        return
    full = {r: rows[(ds, r, aln, "base_raxmlng")]["cpu_s"] for r in reps}
    truef = lambda r: os.path.join(rt.MLDATA, ds, "R%d" % r, "true_tree.tre")  # noqa: E731
    out.append("\n### %s (%s): equal-CPU comparisons\n" % (ds, aln))
    out.append("RAxML-NG full search: mean FN %.2f%%, mean CPU %.1f min (n = %d).\n" % (
        statistics.mean(rows[(ds, r, aln, "base_raxmlng")]["fn"] * 100 for r in reps),
        statistics.mean(full.values()) / 60, len(reps)))
    out.append("| pipeline / method | n | mean FN | mean CPU min (fraction of RAxML-NG) | comparator | "
               "comparator mean FN | ΔFN (method − comparator) | W/T/L | Wilcoxon p |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    arms = [PRIMARY, "constr_ft_0.5", "base_iqfast", "base_fasttree", "place_ft_0.5_fix_ft", "place_ft_0.5_fix_graft"]
    for arm in arms:
        a, b, b_iq, fr = {}, {}, {}, []
        for r in reps:
            row = rows.get((ds, r, aln, arm))
            if not row:
                continue
            a[r] = row["fn"] * 100
            b[r] = anytime_fn(ds, r, aln, "raxmlng", row["cpu_s"])
            if (ds, r, aln, "base_iqtree") in rows:
                b_iq[r] = anytime_fn(ds, r, aln, "iqtree", row["cpu_s"])
            fr.append(row["cpu_s"] / full[r])
        if not a:
            continue
        cpu = statistics.mean(rows[(ds, r, aln, arm)]["cpu_s"] for r in a) / 60
        for lab, comp in (("RAxML-NG@same CPU", b), ("IQ-TREE@same CPU", b_iq)):
            pr = paired(a, comp)
            if pr is None:
                continue
            cm = statistics.mean(comp[r] for r in pr["d"])
            bold = "**" if arm == PRIMARY and lab.startswith("RAxML") else ""
            out.append("| %s%s%s | %d | %.2f%% | %.1f (%.2f) | %s | %.2f%% | %s |" % (
                bold, arm, bold, pr["n"], statistics.mean(a[r] for r in pr["d"]), cpu, statistics.mean(fr), lab, cm,
                fmt(pr)))
            if arm == PRIMARY and lab.startswith("RAxML"):
                out.append("")
                out.append("<!-- primary -->")
                per = ", ".join("R%d %.1f vs %.1f" % (r, a[r], comp[r]) for r in pr["d"])
                out.append("")
                out.append("Primary per replicate (pipeline vs RAxML-NG@T_pipe, FN %%): %s\n" % per)
                out.append("| pipeline / method | n | mean FN | mean CPU min (fraction of RAxML-NG) | comparator | "
                           "comparator mean FN | ΔFN (method − comparator) | W/T/L | Wilcoxon p |")
                out.append("|---|---|---|---|---|---|---|---|---|")
    # fixed fractions
    out.append("\nAnytime FN at fixed fractions of each replicate's full RAxML-NG CPU (paired vs RAxML-NG full):\n")
    out.append("| method @ budget | n | mean FN | ΔFN vs RAxML-NG full | W/T/L | p |")
    out.append("|---|---|---|---|---|---|")
    ref = {r: rows[(ds, r, aln, "base_raxmlng")]["fn"] * 100 for r in reps}
    for m in ("raxmlng", "iqtree"):
        for f in (1 / 3, 0.5, 1.0):
            v = {}
            for r in reps:
                if m == "iqtree" and (ds, r, aln, "base_iqtree") not in rows:
                    continue
                v[r] = anytime_fn(ds, r, aln, m, f * full[r] + (1e-6 if f < 1 else 1e9 if m == "raxmlng" else 0))
            pr = paired(v, ref)
            if pr is None:
                continue
            out.append("| %s @ %.2f× | %d | %.2f%% | %s |" % (m, f, pr["n"], statistics.mean(v[r] for r in pr["d"]),
                                                             fmt(pr)))
    # curves
    for m in ("raxmlng", "iqtree"):
        for f in FRACS:
            for r in reps:
                if m == "iqtree" and (ds, r, aln, "base_iqtree") not in rows:
                    continue
                T = f * full[r] if f < 1 else (1e9 if m == "raxmlng" else full[r])
                curves.append((ds, aln, m, round(f, 3), r, anytime_fn(ds, r, aln, m, T)))
    for arm in arms:
        for r in reps:
            row = rows.get((ds, r, aln, arm))
            if row:
                curves.append((ds, aln, arm, round(row["cpu_s"] / full[r], 3), r, row["fn"] * 100))


def section_all(rows, ds, aln, out):
    arms = sorted({m for (d, r, a, m) in rows if d == ds and a == aln})
    if not arms:
        return
    out.append("\n### %s (%s): all arms\n" % (ds, aln))
    out.append("| arm | n | mean FN | mean FP | per-rep FN | mean CPU min | mean wall min | max peak RSS MB | backbone FN | graft FN |")
    out.append("|---|---|---|---|---|---|---|---|---|---|")
    for m in arms:
        rs = [rows[k] for k in sorted(rows) if k[0] == ds and k[2] == aln and k[3] == m]
        g = lambda k: ("%.1f" % (100 * statistics.mean(x[k] for x in rs))) if all(k in x for x in rs) else "–"  # noqa
        out.append("| %s | %d | %.2f%% | %.2f%% | %s | %.1f | %.1f | %.0f | %s | %s |" % (
            m, len(rs), 100 * statistics.mean(x["fn"] for x in rs), 100 * statistics.mean(x["fp"] for x in rs),
            " ".join("R%d:%.1f" % (x["rep"], 100 * x["fn"]) for x in rs),
            statistics.mean(x["cpu_s"] for x in rs) / 60, statistics.mean(x["wall_s"] for x in rs) / 60,
            max(x["peak_rss_mb"] for x in rs), g("backbone_fn"), g("graft_fn")))


def section_10k(rows, out):
    ds, aln = "RNASim10KHF", "true_mask"
    for (d, r, a, m) in sorted(rows):
        if d == ds and a == aln and m.startswith("base_raxmlngcap"):
            sd = snapdir(d, r, a, m[5:])
            true = os.path.join(rt.MLDATA, d, "R%d" % r, "true_tree.tre")
            out.append("\nRAxML-NG anytime on %s R%d (cap %s CPU-s):\n" % (d, r, m[len("base_raxmlngcap"):]))
            out.append("| CPU h | FN |")
            out.append("|---|---|")
            for T in (600, 1800, 3600, 5400, 7200, 10800, 14400):
                t = anytime.at(sd, T)
                if t:
                    out.append("| %.2f | %.2f%% |" % (T / 3600, fn_of(true, t)))


def main():
    rows = load()
    out = ["<!-- generated by code/summarize.py -->"]
    curves = []
    section_equal_time(rows, "M1HF", "true_align", out, curves)
    for ds, aln in (("M1HF", "witch"), ("RNASimHF", "witch"), ("RNASimHF", "true_align")):
        section_equal_time(rows, ds, aln, out, curves)
    for ds, aln in sorted({(k[0], k[2]) for k in rows}):
        section_all(rows, ds, aln, out)
    section_10k(rows, out)
    open(os.path.join(RES, "tables.md"), "w").write("\n".join(out) + "\n")
    with open(os.path.join(RES, "anytime.tsv"), "w") as f:
        f.write("dataset\taln\tmethod\tcpu_frac\trep\tfn\n")
        for c in curves:
            f.write("\t".join("" if x is None else str(x) for x in c) + "\n")
    print("\n".join(out))
    if "--plot" in sys.argv:
        plot(curves)


def plot(curves):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    groups = sorted({(c[0], c[1]) for c in curves})
    fig, axes = plt.subplots(1, len(groups), figsize=(5.2 * len(groups), 4), squeeze=False)
    for ax, (ds, aln) in zip(axes[0], groups):
        for m, col in (("raxmlng", "#1f5fa8"), ("iqtree", "#c2571a")):
            pts = collections.defaultdict(list)
            for c in curves:
                if c[:3] == (ds, aln, m) and c[5] is not None:
                    pts[c[3]].append(c[5])
            if pts:
                xs = sorted(pts)
                ax.plot(xs, [statistics.mean(pts[x]) for x in xs], "-o", ms=3, color=col,
                        label={"raxmlng": "RAxML-NG (anytime)", "iqtree": "IQ-TREE 3 (anytime)"}[m])
        for arm, mk, lab in ((PRIMARY, "*", "place+graft+RAxML-NG fast polish"), ("constr_ft_0.5", "s", "constrained RAxML-NG"),
                             ("base_iqfast", "^", "IQ-TREE --fast"), ("base_fasttree", "v", "FastTree")):
            xs = [c[3] for c in curves if c[:3] == (ds, aln, arm)]
            ys = [c[5] for c in curves if c[:3] == (ds, aln, arm)]
            if xs:
                ax.plot([statistics.mean(xs)], [statistics.mean(ys)], mk, ms=10, color="#333", label=lab)
        ax.set_xlabel("CPU time / full RAxML-NG CPU (same replicate)")
        ax.set_ylabel("mean FN (%)")
        ax.set_title("%s, %s alignment" % (ds, aln))
        ax.set_ylim(15, 55)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "anytime.png"), dpi=120)


if __name__ == "__main__":
    main()
