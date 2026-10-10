"""Paired summary of bbe.py results: every variant vs the control (MAGUS's own backbones, `linsi`).

    python3 summarize.py RESULTS_DIR [--reps BBA0039_R0,...]   -> markdown on stdout

RESULTS_DIR holds <rep>.results.jsonl and <rep>.diag.jsonl (copied from the work dirs).
Delta = variant - control in error points (x100); W/T/L = variant better / within 0.05 / worse.
"""

import collections
import glob
import json
import os
import sys

from scipy.stats import wilcoxon

TIE = 0.05


def load(d, kind):
    rows = collections.defaultdict(dict)
    for f in sorted(glob.glob(os.path.join(d, "*.{}.jsonl".format(kind)))):
        for line in open(f):
            r = json.loads(line)
            rows[r["rep"]][r["variant"]] = r
    return rows


def paired(res, variant, key="avgErr", reps=None):
    ds = []
    for rep, vs in sorted(res.items()):
        if reps and rep not in reps:
            continue
        if variant in vs and "linsi" in vs:
            ds.append((rep, 100 * (vs[variant][key] - vs["linsi"][key])))
    return ds


def wtl(ds):
    w = sum(d < -TIE for _, d in ds)
    l = sum(d > TIE for _, d in ds)
    return w, len(ds) - w - l, l


def pval(ds):
    xs = [d for _, d in ds]
    if len(xs) < 3 or all(abs(x) < 1e-12 for x in xs):
        return float("nan")
    try:
        return wilcoxon(xs).pvalue
    except ValueError:
        return float("nan")


def table(res, variants, reps=None, title=""):
    out = ["", "#### " + title if title else "", "",
           "| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for v in variants:
        ds = paired(res, v, reps=reps)
        if not ds:
            continue
        fn = paired(res, v, "SPFN", reps)
        fp = paired(res, v, "SPFP", reps)
        rows = [res[r][v] for r, _ in ds]
        mean_err = sum(100 * r["avgErr"] for r in rows) / len(rows)
        bbw = [r["bb_wall"] for r in rows if r.get("bb_wall")]
        mw = sum(r["merge_wall"] for r in rows) / len(rows)
        w, t, l = wtl(ds)
        out.append("| {} | {} | {:.2f} | {:+.2f} | {}/{}/{} | {:.3g} | {:+.2f} | {:+.2f} | {} | {:.0f} |".format(
            v.replace("|", "\\|"), len(ds), mean_err, sum(d for _, d in ds) / len(ds), w, t, l, pval(ds),
            sum(d for _, d in fn) / len(fn), sum(d for _, d in fp) / len(fp),
            "{:.0f}".format(sum(bbw) / len(bbw)) if bbw else "–", mw))
    return "\n".join(out)


def per_rep(res, variants, reps=None):
    reps = sorted(r for r in res if not reps or r in reps)
    out = ["", "| variant | " + " | ".join(r.replace("_R0", "") for r in reps) + " |",
           "|---|" + "---|" * len(reps)]
    for v in ["linsi"] + variants:
        cells = []
        for r in reps:
            if v in res[r]:
                e = 100 * res[r][v]["avgErr"]
                cells.append("{:.2f}".format(e) if v == "linsi" else "{:+.2f}".format(e - 100 * res[r]["linsi"]["avgErr"]))
            else:
                cells.append("")
        out.append("| {} | {} |".format(v.replace("|", "\\|"), " | ".join(cells)))
    return "\n".join(out)


if __name__ == "__main__":
    d = sys.argv[1]
    reps = None
    if "--reps" in sys.argv:
        reps = set(sys.argv[sys.argv.index("--reps") + 1].split(","))
    res = load(d, "results")
    variants = []
    for vs in res.values():
        for v in vs:
            if v != "linsi" and v not in variants:
                variants.append(v)
    print(table(res, variants, reps))
    print(per_rep(res, variants, reps))
