"""Pre-registered DISCO-data held-out test (results/PREREG_disco_heldout.md).
Usage: python analyze_prereg.py out.md setA1.jsonl,setA2.jsonl [setB1.jsonl,setB2.jsonl]"""
import json
import sys

import numpy as np
from scipy.stats import wilcoxon


def load(spec, keep):
    rows, seen = [], set()
    for fn in spec.split(","):
        try:
            lines = open(fn).readlines()
        except FileNotFoundError:
            continue
        for l in lines:
            r = json.loads(l)
            k = (r["cond"], r["rep"], r["sqln"], r["ngen"])
            if k not in seen and keep(r):
                seen.add(k)
                rows.append(r)
    return rows


A = load(sys.argv[2], lambda r: r["cond"] in ("default", "gdl_1e-9_1", "gdl_1e-9_05") and r["rep"] in
         ("05", "06", "07", "08", "09", "10"))
B = load(sys.argv[3], lambda r: r["cond"] in ("gdl_5e-10_05", "ils_1e4", "ils_2e8")) if len(sys.argv) > 3 else []


def fnr(r, m):
    v = r["methods"].get(m)
    return None if v is None else v["FN"] / v["nI"]


def test(rows, a, b):
    x = np.array([fnr(r, a) for r in rows]); y = np.array([fnr(r, b) for r in rows])
    d = x - y
    w, t, lo = int((d < -1e-9).sum()), int((abs(d) <= 1e-9).sum()), int((d > 1e-9).sum())
    p = wilcoxon(x, y).pvalue if w + lo else 1.0
    return d.mean(), w, t, lo, p, x.mean(), y.mean()


L = ["# Pre-registered held-out test on DISCO data (see PREREG_disco_heldout.md)\n"]
for name, rows in [("A (reps 05-10: default, gdl_1e-9_1, gdl_1e-9_05)", A), ("B (new conditions, reps 01-05)", B),
                   ("A ∪ B (primary)", A + B)]:
    if not rows:
        continue
    L.append("\n## Set %s: n = %d runs\n" % (name, len(rows)))
    L.append("| comparison | mean FN (Pro) | mean FN (base) | mean diff | W/T/L | p | Holm-adjusted p |")
    L.append("|---|---|---|---|---|---|---|")
    prim = [("astrid-pro", b) for b in ("astrid-multi", "astrid-disco", "astral-pro")]
    res = [(c, test(rows, *c)) for c in prim]
    ps = sorted(range(len(res)), key=lambda i: res[i][1][4])
    adj = [None] * len(res)
    run = 0.0
    for rank, i in enumerate(ps):
        run = max(run, min(1.0, (len(res) - rank) * res[i][1][4]))
        adj[i] = run
    for (c, (dm, w, t, lo, p, ma, mb)), pa in zip(res, adj):
        L.append("| %s vs %s | %.4f | %.4f | %+.4f | %d/%d/%d | %.3g | %.3g |" % (c[0], c[1], ma, mb, dm, w, t, lo, p, pa))
    for a in ("astrid-pro-w", "astrid-pro-min"):
        for b in ("astrid-multi", "astrid-disco", "astral-pro"):
            if all(fnr(r, a) is not None for r in rows):
                dm, w, t, lo, p, ma, mb = test(rows, a, b)
                L.append("| (secondary) %s vs %s | %.4f | %.4f | %+.4f | %d/%d/%d | %.3g | - |" % (a, b, ma, mb, dm, w, t, lo, p))
    L.append("\nPer-condition mean FN rate:\n")
    ms = ["astrid-multi", "astrid-pro", "astrid-pro-w", "astrid-pro-min", "astrid-disco", "astral-pro"]
    L.append("| condition | n | " + " | ".join(ms) + " |")
    L.append("|---|---|" + "---|" * len(ms))
    for c in sorted({r["cond"] for r in rows}):
        sub = [r for r in rows if r["cond"] == c]
        L.append("| %s | %d | " % (c, len(sub)) + " | ".join(
            "%.3f" % np.mean([fnr(r, m) for r in sub if fnr(r, m) is not None]) if any(fnr(r, m) is not None for r in sub) else "-"
            for m in ms) + " |")
open(sys.argv[1], "w").write("\n".join(L) + "\n")
print("\n".join(L))
