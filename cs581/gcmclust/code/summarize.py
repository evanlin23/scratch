"""Tables for REPORT.md from results/*.results.jsonl.

    python3 summarize.py [train|held|all] [--select]   # markdown to stdout
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
SPLIT = {}
for line in open(os.path.join(HERE, "reps.txt")):
    if line.strip() and not line.startswith("#"):
        n, s = line.split()[:2]
        SPLIT[n] = s
PROT = lambda r: r.startswith("BBA") or r.startswith("SIM")  # noqa: E731
BASE = "raw:mcl:4"
for _i, _a in enumerate(sys.argv):
    if _a == "--base":
        BASE = sys.argv[_i + 1]


def load():
    rows = defaultdict(dict)  # rep -> variant -> row
    for f in glob.glob(os.path.join(RES, "*.results.jsonl")):
        for line in open(f):
            r = json.loads(line)
            rows[r["rep"]][r["variant"]] = r
    return rows


def family(v):
    g, m, *p = v.split(":")
    return g, m


def table(rows, reps, variants, title):
    out = ["#### " + title, "",
           "| variant | n | mean Δ err vs `{}` | W/T/L | p | ΔSPFN / ΔSPFP | cluster s + trace s (base total s) | clusters / singleton frac |".format(BASE),
           "|---|---|---|---|---|---|---|---|"]
    for v in variants:
        d, fn, fp, mw, bw, nc, sg, fail = [], [], [], [], [], [], [], 0
        for r in reps:
            x, b = rows[r].get(v), rows[r].get(BASE)
            if not x or not b:
                continue
            if "failed" in x:
                fail += 1
                continue
            d.append(x["err"] - b["err"])
            fn.append(100 * (x["SPFN"] - b["SPFN"]))
            fp.append(100 * (x["SPFP"] - b["SPFP"]))
            mw.append((x.get("cluster_s", 0), x["merge_wall"]))
            bw.append(b["merge_wall"] + b.get("cluster_s", 0))
            nc.append(x["n_clusters"])
            sg.append(x["frac_singleton"])
        if not d and not fail:
            continue
        d = np.array(d)
        w, t, l = int((d < -0.05).sum()), int((abs(d) <= 0.05).sum()), int((d > 0.05).sum())
        try:
            p = "{:.3f}".format(wilcoxon(d).pvalue) if len(d) >= 5 and np.any(d != 0) else "–"
        except ValueError:
            p = "–"
        f = " ({} failed)".format(fail) if fail else ""
        out.append("| `{}` | {}{} | {:+.2f} | {}/{}/{} | {} | {:+.2f} / {:+.2f} | {} ({:.0f}) | {:.0f} / {:.3f} |".format(
            v, len(d), f, d.mean() if len(d) else float("nan"), w, t, l, p, np.mean(fn) if fn else 0,
            np.mean(fp) if fp else 0, "{:.0f} + {:.0f}".format(*np.mean(mw, axis=0)) if mw else "–", np.mean(bw) if bw else 0,
            np.mean(nc) if nc else 0, np.mean(sg) if sg else 0))
    return "\n".join(out) + "\n"


def per_rep(rows, reps, variants):
    out = ["| replicate | MAGUS err | " + " | ".join("`{}`".format(v) for v in variants) + " |",
           "|---|---|" + "---|" * len(variants)]
    for r in reps:
        b = rows[r].get(BASE)
        if not b:
            continue
        cells = []
        for v in variants:
            x = rows[r].get(v)
            cells.append("–" if not x else ("fail" if "failed" in x else "{:+.2f}".format(x["err"] - b["err"])))
        out.append("| {} | {:.2f} | {} |".format(r, b["err"], " | ".join(cells)))
    return "\n".join(out) + "\n"


def select(rows, reps):
    """PREREG selection: per family x graph, lowest mean Δ over training reps (complete, no failures)."""
    allv = sorted({v for r in reps for v in rows[r]} - {BASE})
    best = {}
    for v in allv:
        d = []
        ok = True
        for r in reps:
            x, b = rows[r].get(v), rows[r].get(BASE)
            if not x or "failed" in x or not b:
                ok = False
                break
            d.append(x["err"] - b["err"])
        if not ok:
            continue
        key = family(v)
        m = float(np.mean(d))
        if key not in best or m < best[key][1] - 0.01:
            best[key] = (v, m)
    return best


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "train"
    rows = load()
    reps = [r for r in SPLIT if (which == "all" or SPLIT[r] == which) and r in rows]
    variants = sorted({v for r in reps for v in rows[r]}, key=lambda v: (v.split(":")[1], v.split(":")[0], v))
    print(table(rows, reps, variants, "{}: all replicates (n = {})".format(which, len(reps))))
    print(table(rows, [r for r in reps if PROT(r)], variants, "{}: proteins".format(which)))
    print(table(rows, [r for r in reps if not PROT(r)], variants, "{}: DNA/RNA".format(which)))
    if "--select" in sys.argv:
        best = select(rows, reps)
        print("#### selected (PREREG rule)\n")
        for k, (v, m) in sorted(best.items()):
            print("- {} / {}: `{}` (training mean Δ {:+.2f})".format(k[1], k[0], v, m))
        print()
    if "--perrep" in sys.argv:
        vs = sys.argv[sys.argv.index("--perrep") + 1].split(",")
        print(per_rep(rows, reps, vs))


def cluster_stats(rows, reps, variants):
    out = ["| variant | clusters (size ≥ 2) | mean size | max size | singleton nodes | clusters violating 1-col/subset | "
           "aln length / ref |", "|---|---|---|---|---|---|---|"]
    for v in variants:
        xs = [rows[r][v] for r in reps if v in rows[r] and "err" in rows[r][v]]
        if not xs:
            continue
        m = lambda k: np.mean([x[k] for x in xs])  # noqa: E731
        out.append("| `{}` | {:.0f} | {:.1f} | {:.0f} | {:.3f} | {:.2f} | {:.2f} |".format(
            v, m("n_clusters"), m("mean_size"), m("max_size"), m("frac_singleton"), m("frac_violating"),
            np.mean([x["LenEst"] / x["LenRef"] for x in xs])))
    return "\n".join(out) + "\n"


if __name__ == "__main__" and "--stats" in sys.argv:
    vs = sys.argv[sys.argv.index("--stats") + 1].split(",")
    for lab, sel in (("proteins", PROT), ("DNA/RNA", lambda r: not PROT(r))):
        print("#### cluster statistics, {} ({})\n".format(lab, which))
        print(cluster_stats(rows, [r for r in reps if sel(r)], vs))
