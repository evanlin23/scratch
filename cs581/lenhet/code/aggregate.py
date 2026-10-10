"""Collect score.json files into results/*.csv and markdown tables.

    python aggregate.py [ROOT=/opt/runs/lenhet] [OUT=../results]

Primary metric for UPP/WITCH (and the -tfa/-trim variants built on them) is
FastSP with lower-case insertion letters masked (-ml): those letters are
declared unaligned by the method. The unmasked ("raw") number is also kept.
"""
import csv
import glob
import json
import os
import statistics
import sys

from scipy.stats import wilcoxon

ROOT = sys.argv[1] if len(sys.argv) > 1 else "/opt/runs/lenhet"
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
COND_ORDER = ["rand_f0.1_m0", "rand_f0.1_m0.5", "rand_f0.1_m1", "rand_f0.1_m3", "dom_f0.1_m0"]
COND_LABEL = {"rand_f0.1_m0": "control (no flank)", "rand_f0.1_m0.5": "random flank +0.5x",
              "rand_f0.1_m1": "random flank +1x", "rand_f0.1_m3": "random flank +3x",
              "dom_f0.1_m0": "second domain (+~1x)"}
METH_ORDER = ["upp", "witch", "emma", "mafft-add", "mafft-addlong", "upp-trim", "upp-tfa", "emma-trim",
              "emma-tfa", "witch-tfa", "mafft", "magus"]


def pick(sc, part):
    return sc.get(part + "_ml", sc.get(part))


def load():
    rows = []
    for p in glob.glob(os.path.join(ROOT, "*", "R*", "*", "score.json")):
        d = json.load(open(p))
        cond, rep, meth = p.split(os.sep)[-4:-1]
        if cond.startswith("valid"):
            continue
        r = {"cond": cond, "rep": rep, "method": meth, "time": d["time"]}
        for part in ("all", "long"):
            s = pick(d["score"], part)
            raw = d["score"][part]
            r[part + "_SPFN"], r[part + "_SPFP"] = s["SPFN"], s["SPFP"]
            r[part + "_err"] = (s["SPFN"] + s["SPFP"]) / 2
            r[part + "_rawSPFP"] = raw["SPFP"]
        rows.append(r)
    return rows


def table(rows, part):
    conds = [c for c in COND_ORDER if any(r["cond"] == c for r in rows)]
    meths = [m for m in METH_ORDER if any(r["method"] == m for r in rows)]
    lines = ["| method | " + " | ".join(COND_LABEL[c] for c in conds) + " |",
             "|---|" + "---|" * len(conds)]
    for m in meths:
        cells = []
        for c in conds:
            v = [r for r in rows if r["cond"] == c and r["method"] == m]
            if not v:
                cells.append("–")
                continue
            fn = statistics.mean(r[part + "_SPFN"] for r in v)
            fp = statistics.mean(r[part + "_SPFP"] for r in v)
            cells.append("%.3f / %.3f (n=%d)" % (fn, fp, len(v)))
        lines.append("| %s | %s |" % (m, " | ".join(cells)))
    return "\n".join(lines)


def rawtable(rows):
    conds = [c for c in COND_ORDER if any(r["cond"] == c for r in rows)]
    lines = ["| method | " + " | ".join(COND_LABEL[c] for c in conds) + " |", "|---|" + "---|" * len(conds)]
    for m in ("upp", "witch", "upp-tfa"):
        cells = []
        for c in conds:
            v = [r["long_rawSPFP"] for r in rows if r["cond"] == c and r["method"] == m]
            cells.append("%.3f" % statistics.mean(v) if v else "–")
        lines.append("| %s | %s |" % (m, " | ".join(cells)))
    return "\n".join(lines)


def runtime(rows):
    conds = [c for c in COND_ORDER if any(r["cond"] == c for r in rows)]
    meths = [m for m in METH_ORDER if any(r["method"] == m for r in rows)]
    lines = ["| method | " + " | ".join(COND_LABEL[c] for c in conds) + " |", "|---|" + "---|" * len(conds)]
    for m in meths:
        cells = []
        for c in conds:
            v = [r["time"] for r in rows if r["cond"] == c and r["method"] == m]
            v1 = [r["time"] for r in rows if r["cond"] == c and r["method"] == m and r["rep"] != "R0"]
            v = v1 or v  # R1/R2 ran 4 single-thread jobs at once; part of R0 ran 2x2 threads
            cells.append("%.0f" % statistics.mean(v) if v else "–")
        lines.append("| %s | %s |" % (m, " | ".join(cells)))
    return "\n".join(lines)


def paired(rows, new, base, part, conds=None):
    idx = {(r["cond"], r["rep"], r["method"]): r for r in rows}
    diffs = []
    for (c, rep, m), r in idx.items():
        if m == new and (conds is None or c in conds) and (c, rep, base) in idx:
            diffs.append(r[part + "_err"] - idx[(c, rep, base)][part + "_err"])
    if not diffs:
        return None
    w = sum(d < -1e-9 for d in diffs)
    l_ = sum(d > 1e-9 for d in diffs)
    t = len(diffs) - w - l_
    nz = [d for d in diffs if abs(d) > 1e-9]
    p = wilcoxon(nz).pvalue if len(nz) >= 1 else float("nan")
    return {"new": new, "base": base, "part": part, "n": len(diffs), "mean_diff": statistics.mean(diffs),
            "W": w, "T": t, "L": l_, "p": p}


def main():
    rows = load()
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "scores.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0]))
        wr.writeheader()
        for r in sorted(rows, key=lambda r: (r["cond"], r["method"], r["rep"])):
            wr.writerow(r)
    md = ["# Results (auto-generated by code/aggregate.py)", "",
          "Cells: mean SPFN / SPFP over replicates (UPP/WITCH-based: insertion letters masked).", "",
          "## Long (lengthened) sequences only", "", table(rows, "long"), "",
          "## All sequences", "", table(rows, "all"), "",
          "## Long sequences, SPFP without masking lower-case insertion letters", "",
          "(what a user gets if the UPP/WITCH output is used as-is, e.g. upper-cased for tree estimation)", "",
          rawtable(rows), "",
          "## Runtime (s, mean; add step only, backbone given; R1-R2: 1 thread, 4 jobs on 4 cores. "
          "mafft-addlong / mafft / magus: R0 only)", "", runtime(rows), "",
          "## Paired comparisons (error = (SPFN+SPFP)/2; W = new better)", "",
          "| new vs base | scored | conditions | n | mean diff | W/T/L | Wilcoxon p |", "|---|---|---|---|---|---|---|"]
    stats = []
    for new, base in (("upp-tfa", "upp"), ("emma-tfa", "emma"), ("upp-trim", "upp"), ("emma-trim", "emma"),
                      ("emma", "upp"), ("witch", "upp"), ("magus", "upp-tfa"), ("mafft", "upp-tfa")):
        for part in ("long", "all"):
            for label, conds in (("all", None), ("random flank", [c for c in COND_ORDER if c.startswith("rand") and c != "rand_f0.1_m0"]),
                                 ("second domain", ["dom_f0.1_m0"]), ("control", ["rand_f0.1_m0"])):
                s = paired(rows, new, base, part, conds)
                if s:
                    s["conds"] = label
                    stats.append(s)
                    md.append("| %s vs %s | %s | %s | %d | %+.4f | %d/%d/%d | %.3g |" % (
                        new, base, part, label, s["n"], s["mean_diff"], s["W"], s["T"], s["L"], s["p"]))
    json.dump(stats, open(os.path.join(OUT, "paired.json"), "w"), indent=1)
    open(os.path.join(OUT, "tables.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    main()
