"""Tables for REPORT.md from results/ (pre-registered analysis, ../PREREG.md).

    python3 summarize.py > ../results/summary.md
"""

import glob
import json
import os

import numpy as np
from scipy.stats import wilcoxon

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
CTRL = "linsi"
METHODS = [("linsi|cons0.7", "primary: L-INS-i \\| cons0.7"), ("linsi&fftns2-op3", "L-INS-i ∩ FFT-NS-2 --op 3"),
           ("linsi+clustalo|cons0.7", "(L-INS-i + Clustal) \\| cons0.7")]
SHORT = {"linsi|cons0.7": "L\\|cons", "linsi&fftns2-op3": "L∩F", "linsi+clustalo|cons0.7": "(L+C)\\|cons"}
TIE = 0.05


def family(n):
    if n.startswith("SIMMOD"):
        return "SIMMOD"
    if n.startswith("SIMHIGH"):
        return "SIMHIGH"
    if n.startswith("10AA"):
        return "10AA"
    if n.startswith("HF_"):
        return "HomFam"
    if n.startswith("BB_"):
        return "BAliBASE (fresh draws)"
    if n.startswith("16S"):
        return "16S.M (in-sample)"
    return "nucleotide (held out)"


PROTEIN = {"SIMMOD", "SIMHIGH", "10AA", "HomFam", "BAliBASE (fresh draws)"}


def load():
    data = {}
    for f in glob.glob(os.path.join(R, "*.results.jsonl")):
        n = os.path.basename(f)[:-len(".results.jsonl")]
        rows = {}
        for l in open(f):
            r = json.loads(l)
            rows[r["variant"]] = r
        if CTRL not in rows:
            continue
        d = {"rows": rows, "fam": family(n)}
        for k in ("gate", "magus", "trees"):
            p = os.path.join(R, "{}.{}.json".format(n, k))
            if os.path.exists(p):
                d[k] = json.load(open(p))
        tp = os.path.join(R, n + ".trees.jsonl")
        if os.path.exists(tp):
            d["trees"] = {r["method"]: r for r in map(json.loads, open(tp))}
        data[n] = d
    return data


def delta(d, m, key="avgErr"):
    r = d["rows"]
    if m not in r:
        return None
    return 100 * (r[m][key] - r[CTRL][key])


def stats(ds):
    ds = [x for x in ds if x is not None]
    if not ds:
        return "| 0 | | | | "
    w = sum(x < -TIE for x in ds)
    l = sum(x > TIE for x in ds)
    t = len(ds) - w - l
    try:
        p = wilcoxon(ds).pvalue if any(abs(x) > 0 for x in ds) and len(ds) > 1 else float("nan")
    except ValueError:
        p = float("nan")
    return "| {} | {:+.2f} | {}/{}/{} | {:.3g} |".format(len(ds), np.mean(ds), w, t, l, p)


def table_groups(data, groups, title):
    print("### " + title + "\n")
    print("| group | method | n | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | Δ TC |")
    print("|---|---|---|---|---|---|---|---|---|")
    for gname, names in groups:
        for m, lab in METHODS:
            ds = [delta(data[n], m) for n in names]
            fn = [delta(data[n], m, "SPFN") for n in names]
            fp = [delta(data[n], m, "SPFP") for n in names]
            tc = [delta(data[n], m, "TC") for n in names]
            mean = lambda v: "{:+.2f}".format(np.mean([x for x in v if x is not None])) if any(x is not None for x in v) else ""
            print("| {} | {} {} {} | {} | {} |".format(gname, lab, stats(ds), mean(fn), mean(fp), mean(tc)))
    print()


def main():
    data = load()
    fams = {}
    for n in sorted(data):
        fams.setdefault(data[n]["fam"], []).append(n)
    order = ["SIMMOD", "SIMHIGH", "10AA", "HomFam", "BAliBASE (fresh draws)", "nucleotide (held out)",
             "16S.M (in-sample)"]
    groups = [(f, fams[f]) for f in order if f in fams]
    prot = [n for f in order if f in PROTEIN for n in fams.get(f, [])]
    sim = fams.get("SIMMOD", []) + fams.get("SIMHIGH", [])
    nonbb = [n for n in prot if data[n]["fam"] != "BAliBASE (fresh draws)"]
    groups_pooled = [("**all held-out protein**", prot), ("simulated (SIMMOD+SIMHIGH)", sim),
                     ("protein excl. BAliBASE", nonbb)]
    print("## Paired Δ vs MAGUS's own merge (error points; negative = better)\n")
    table_groups(data, groups_pooled, "Pooled")
    table_groups(data, groups, "Per dataset family")

    print("### Per dataset\n")
    print("| dataset | family | MAGUS err | merge ctrl err | " + " | ".join("Δ " + SHORT[m] for m, _ in METHODS)
          + " | ΔSPFN / ΔSPFP (L\\|cons) | support | gate | frac L-only |")
    print("|---|---|---|---|" + "---|" * len(METHODS) + "---|---|---|---|")
    for f, names in groups:
        for n in names:
            d = data[n]
            mg = d.get("magus", {})
            me = mg.get("magus", mg).get("avgErr")
            g = d.get("gate", {})
            cells = []
            for m, _ in METHODS:
                x = delta(d, m)
                cells.append("" if x is None else "{:+.2f}".format(x))
            fn, fp = delta(d, METHODS[0][0], "SPFN"), delta(d, METHODS[0][0], "SPFP")
            print("| {} | {} | {} | {:.2f} | {} | {} | {} | {} | {} |".format(
                n, f, "" if me is None else "{:.2f}".format(100 * me), 100 * d["rows"][CTRL]["avgErr"],
                " | ".join(cells), "" if fn is None else "{:+.2f} / {:+.2f}".format(fn, fp),
                g.get("support_a_only", ""), {True: "filter", False: "keep"}.get(g.get("filter"), ""),
                g.get("frac_a_only", "")))
    print()

    # gate
    print("## Gate (support < 0.615 → filter) vs observed effect of the primary method\n")
    print("Observed help = Δ(L\\|cons0.7) < 0 (strict sign).\n")
    rows = [(n, data[n]) for f, names in groups for n in names if "gate" in data[n] and delta(data[n], METHODS[0][0]) is not None]
    for label, sel in (("all", rows), ("protein", [r for r in rows if r[1]["fam"] in PROTEIN]),
                       ("nucleotide (held out + 16S)", [r for r in rows if r[1]["fam"] not in PROTEIN])):
        tp = fp = tn = fn = 0
        for n, d in sel:
            pred = d["gate"]["filter"]
            obs = delta(d, METHODS[0][0]) < 0
            tp += pred and obs; fp += pred and not obs; tn += (not pred) and (not obs); fn += (not pred) and obs
        tot = tp + fp + tn + fn
        if tot:
            print("- {}: n = {}; accuracy {}/{} = {:.0%} (filter & helped {}, filter & hurt {}, keep & would-hurt {}, "
                  "keep & would-help {})".format(label, tot, tp + tn, tot, (tp + tn) / tot, tp, fp, tn, fn))
    print()
    print("Gated policy (primary method if the gate says filter, else MAGUS's own merge), Δ vs MAGUS:\n")
    print("| group | n | mean Δ | W/T/L | p | n filtered | mean Δ if always filtering |")
    print("|---|---|---|---|---|---|---|")
    for gname, names in groups_pooled + groups:
        ds, al, nf = [], [], 0
        for n in names:
            d = data[n]
            if "gate" not in d or delta(d, METHODS[0][0]) is None:
                continue
            x = delta(d, METHODS[0][0])
            al.append(x)
            nf += d["gate"]["filter"]
            ds.append(x if d["gate"]["filter"] else 0.0)
        if ds:
            s = stats(ds)
            print("| {} {} {} | {:+.2f} |".format(gname, s, nf, np.mean(al)))
    print()

    # trees
    tr = [(n, data[n]["trees"]) for n in sim if "trees" in data[n]]
    if tr:
        print("## Tree error (FastTree -lg -gamma, nRF to the true tree, %)\n")
        ms = ["true", "magus", "cons", "fftint", "lccons"]
        print("| dataset | " + " | ".join(ms) + " |")
        print("|---|" + "---|" * len(ms))
        for n, t in tr:
            print("| {} | ".format(n) + " | ".join("{:.2f}".format(100 * t[m]["RF"]) if m in t else "" for m in ms) + " |")
        for m in ms[2:]:
            ds = [100 * (t[m]["RF"] - t["magus"]["RF"]) for n, t in tr if m in t and "magus" in t]
            print("\nΔ nRF {} − MAGUS: n = {}, mean {:+.2f}, W/T/L {}/{}/{}{}".format(
                m, len(ds), np.mean(ds), sum(x < -TIE for x in ds), sum(abs(x) <= TIE for x in ds),
                sum(x > TIE for x in ds),
                ", p = {:.3g}".format(wilcoxon(ds).pvalue) if len(ds) > 1 and any(ds) else ""))
        print()

    # runtime
    print("## Runtime (seconds, 4 threads; one job at a time for proteins)\n")
    print("Extra = prep (new backbone alignments + masking/intersection) + (merge − control merge). "
          "Gate = Clustal backbones + support statistic.\n")
    print("| dataset | MAGUS e2e | ctrl merge | " + " | ".join("extra " + SHORT[m] for m, _ in METHODS) + " | gate |")
    print("|---|---|---|" + "---|" * len(METHODS) + "---|")
    for f, names in groups:
        for n in names:
            d = data[n]
            mg = d.get("magus", {})
            wall = mg.get("magus", {}).get("wall") if "magus" in mg else mg.get("seconds")
            c = d["rows"][CTRL]
            ex = []
            for m, _ in METHODS:
                r = d["rows"].get(m)
                ex.append("" if r is None else "{:.0f}".format(r["prep_wall"] + r["merge_wall"] - c["merge_wall"]))
            g = d.get("gate", {})
            gt = (g.get("clustalo_bb_wall") or 0) + g.get("gate_wall", 0) if g else ""
            print("| {} | {} | {:.0f} | {} | {} |".format(n, "" if wall is None else "{:.0f}".format(wall),
                                                       c["merge_wall"], " | ".join(ex),
                                                       "" if gt == "" else "{:.0f}".format(gt)))
    print("\nNucleotide MAGUS times are from the cached worker runs (contended); not comparable.")


if __name__ == "__main__":
    main()
