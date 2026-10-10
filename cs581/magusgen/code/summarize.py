"""Tables for REPORT.md from results/*.results.jsonl (paired vs the `linsi` control = MAGUS).

    python3 summarize.py > ../results/tables.md
"""

import glob
import json
import os

import numpy as np
from scipy.stats import wilcoxon

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
CTRL = "linsi"
TIE = 0.05


def dtype(n):
    if n.startswith(("BBA", "SIM")):
        return "protein"
    return "DNA/RNA"


def load():
    data = {}
    for f in sorted(glob.glob(os.path.join(R, "*.results.jsonl"))):
        n = os.path.basename(f)[:-len(".results.jsonl")]
        rows = {}
        for l in open(f):
            r = json.loads(l)
            rows[r["variant"]] = r
        if CTRL in rows:
            data[n] = rows
    return data


def overhead(rows, v):
    """Extra wall seconds over MAGUS's own merge (MAGUS's backbones/subsets are reused)."""
    r = rows[v]
    if v.startswith(("ss:", "ssu:")):
        base = v.split(":", 1)[1]
        return overhead(rows, base) + r.get("ev_wall", 0) + r.get("split_wall", 0) + r["merge_wall"]
    if v == CTRL:
        return 0.0
    return (r.get("bb_wall") or 0) + (r.get("prep_wall") or 0) + r["merge_wall"] - rows[CTRL]["merge_wall"]


def stats(ds, key="avgErr"):
    """Mean over scored pairs; a timed-out run (no score) counts as a loss in W/T/L and is left out of
    the mean and the Wilcoxon test."""
    fails = sum(1 for x in ds if key not in x[1])
    d = np.array([100 * (x[1][key] - x[0][key]) for x in ds if key in x[1]])
    if len(d) == 0:
        return None
    w = int((d < -TIE).sum()); l = int((d > TIE).sum()) + fails; t = len(d) + fails - w - l
    try:
        p = wilcoxon(d).pvalue if np.any(d != 0) and len(d) > 1 else 1.0
    except ValueError:
        p = 1.0
    return d.mean(), w, t, l, p, len(d) + fails


def main():
    data = load()
    variants = []
    for rows in data.values():
        for v in rows:
            if v != CTRL and v not in variants:
                variants.append(v)
    print("## Per dataset: MAGUS error and Δ error (points) of each variant\n")
    print("| dataset | type | MAGUS err | " + " | ".join("`{}`".format(v.replace("|", "\\|")) for v in variants) + " |")
    print("|---|---|---|" + "---|" * len(variants))
    for n, rows in data.items():
        c = rows[CTRL]["avgErr"]
        cells = [("timeout" if "avgErr" not in rows[v] else "{:+.2f}".format(100 * (rows[v]["avgErr"] - c)))
                 if v in rows else "" for v in variants]
        print("| {} | {} | {:.2f} | {} |".format(n, dtype(n), 100 * c, " | ".join(cells)))
    print("\n## Summary per variant (Δ = variant − MAGUS, points; W/T/L with |Δ| < {} a tie; two-sided Wilcoxon)\n".format(TIE))
    print("| variant | group | n | Δ error | W/T/L | p | Δ SPFN | Δ SPFP | extra wall s over MAGUS (median) |")
    print("|---|---|---|---|---|---|---|---|---|")
    for v in variants:
        for g in ("DNA/RNA", "protein", "pooled"):
            ds = [(rows[CTRL], rows[v]) for n, rows in data.items() if v in rows and (g == "pooled" or dtype(n) == g)]
            s = stats(ds)
            if not s:
                continue
            sn, sp = stats(ds, "SPFN"), stats(ds, "SPFP")
            ns = [n for n, rows in data.items() if v in rows and (g == "pooled" or dtype(n) == g)]
            oh = [overhead(data[n], v) for n in ns]
            print("| `{}` | {} | {} | {:+.2f} | {}/{}/{} | {:.3g} | {:+.2f} | {:+.2f} | {:.0f} |".format(
                v.replace("|", "\\|"), g, s[5], s[0], s[1], s[2], s[3], s[4], sn[0], sp[0], np.median(oh)))


if __name__ == "__main__":
    main()
