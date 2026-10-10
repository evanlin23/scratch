"""Tables for the DISCO-data runs (estimated gene trees g_100, or true trees for the
gtrees_10000 curve). Usage: python analyze_disco.py out.md runs.jsonl [curve.jsonl]"""
import collections
import json
import sys

import numpy as np
from scipy.stats import wilcoxon

out = sys.argv[1]
rows, _seen = [], set()
for fn in sys.argv[2].split(","):
    for l in open(fn):
        r = json.loads(l)
        k = (r["cond"], r["rep"], r["sqln"], r["ngen"])
        if k not in _seen:
            _seen.add(k)
            rows.append(r)
L = []
SHOW = ["astrid-multi", "astrid-multi-w", "astrid-pro", "astrid-pro-min", "astrid-pro-w", "ortho-allnodes",
        "astrid-disco", "astral-pro"]


def fnr(r, m):
    v = r["methods"].get(m)
    return None if v is None else v["FN"] / v["nI"]


L.append("# DISCO-data comparison (doi:10.13012/B2IDB-4050038_V1), estimated gene trees from 100 bp\n")
L.append("Mean species-tree FN rate (100 species). `-w` = support-weighted. Conditions: default (dup 5e-10, "
         "loss/dup 1, 1000 genes), gdl_1e-9_1 (dup 1e-9, loss/dup 1, 1000 genes), gdl_1e-9_05 (dup 1e-9, "
         "loss/dup 0.5, i.e. supercritical; only the first 100 genes because families average ~1000 leaves).\n")
L.append("| condition | genes | n | " + " | ".join(SHOW) + " |")
L.append("|---|---|---|" + "---|" * len(SHOW))
for c in ["default", "gdl_1e-9_1", "gdl_1e-9_05"]:
    sub = [r for r in rows if r["cond"] == c]
    if not sub:
        continue
    vals = []
    for m in SHOW:
        v = [fnr(r, m) for r in sub if fnr(r, m) is not None]
        vals.append("%.3f" % np.mean(v) if v else "-")
    L.append("| %s | %d | %d | " % (c, sub[0]["ngen"], len(sub)) + " | ".join(vals) + " |")

L.append("\n## Paired tests, all conditions pooled (W/T/L from the first method's view; tie = equal FN)\n")
L.append("| method | baseline | n | mean diff | W/T/L | p |")
L.append("|---|---|---|---|---|---|")
for a in ["astrid-pro-min", "astrid-pro", "astrid-pro-w"]:
    for b in ["astrid-multi", "astrid-multi-w", "astrid-disco", "astral-pro"]:
        x, y = [], []
        for r in rows:
            fa, fb = fnr(r, a), fnr(r, b)
            if fa is not None and fb is not None:
                x.append(fa); y.append(fb)
        if not x:
            continue
        x, y = np.array(x), np.array(y)
        d = x - y
        w, t, lo = int((d < -1e-9).sum()), int((abs(d) <= 1e-9).sum()), int((d > 1e-9).sum())
        p = wilcoxon(x, y).pvalue if w + lo else 1.0
        L.append("| %s | %s | %d | %+.4f | %d/%d/%d | %.2g |" % (a, b, len(x), d.mean(), w, t, lo, p))

secs = collections.defaultdict(list)
for r in rows:
    for k, v in r["secs"].items():
        secs[(r["cond"], k)].append(v)
L.append("\n## Runtime (s, mean per run, single core)\n")
L.append("| condition | " + " | ".join(["root+tag", "astrid-multi", "astrid-pro", "astrid-disco", "astral-pro"]) + " |")
L.append("|---|---|---|---|---|---|")
for c in ["default", "gdl_1e-9_1", "gdl_1e-9_05"]:
    if secs[(c, "astrid-multi")]:
        L.append("| %s | " % c + " | ".join("%.1f" % np.mean(secs[(c, k)]) for k in
                                         ["root+tag", "astrid-multi", "astrid-pro", "astrid-disco", "astral-pro"]) + " |")

if len(sys.argv) > 3:
    cr = [json.loads(l) for l in open(sys.argv[3])]
    L.append("\n## True gene trees with GDL + ILS (gtrees_10000_l1): FN rate vs number of genes\n")
    ms = ["astrid-multi", "astrid-pro", "astrid-pro-min", "astrid-disco", "astral-pro"]
    L.append("| genes | n reps | " + " | ".join(ms) + " |")
    L.append("|---|---|" + "---|" * len(ms))
    for n in sorted({r["ngen"] for r in cr}):
        sub = [r for r in cr if r["ngen"] == n]
        L.append("| %d | %d | " % (n, len(sub)) + " | ".join(
            ("%.3f" % np.mean([fnr(r, m) for r in sub if fnr(r, m) is not None]))
            if any(fnr(r, m) is not None for r in sub) else "-" for m in ms) + " |")
open(out, "w").write("\n".join(L) + "\n")
print("\n".join(L))
