"""Paired comparison tables for the FastMulRFS-data runs.
Dev replicates 01-03 pick the ASTRID-Pro variant (mean vs min); held-out replicates
04-10 are used for the final comparison. Metric: species-tree FN rate (FN / internal
branches of the true tree; trees are binary so FN rate = RF rate). Ties: equal FN count.
Usage: python analyze.py out.md in1.jsonl [in2.jsonl ...]"""
import collections
import json
import sys

import numpy as np
from scipy.stats import wilcoxon

out = sys.argv[1]
rows = []
for fn in sys.argv[2:]:
    for l in open(fn):
        rows.append(json.loads(l))
DEV = {"01", "02", "03"}


def fnr(r, m):
    v = r["methods"].get(m)
    return None if v is None else v["FN"] / v["nI"]


def paired(rs, a, b):
    x, y = [], []
    for r in rs:
        fa, fb = fnr(r, a), fnr(r, b)
        if fa is not None and fb is not None:
            x.append(fa); y.append(fb)
    x, y = np.array(x), np.array(y)
    if len(x) == 0:
        return None
    d = x - y
    w = int((d < -1e-9).sum()); t = int((np.abs(d) <= 1e-9).sum()); lo = int((d > 1e-9).sum())
    p = wilcoxon(x, y, zero_method="wilcox", alternative="two-sided").pvalue if (w + lo) > 0 else 1.0
    return dict(n=len(x), mean_a=x.mean(), mean_b=y.mean(), diff=d.mean(), W=w, T=t, L=lo, p=p)


lines = []
dev = [r for r in rows if r["rep"] in DEV]
hold = [r for r in rows if r["rep"] not in DEV]
lines.append("# FastMulRFS-data comparison (estimated gene trees)\n")
lines.append("Runs: %d (dev %d, held-out %d). Metric: species-tree FN rate (= RF rate; all trees binary). "
             "W/T/L = ASTRID-Pro better / tie (equal FN count) / worse. Two-sided Wilcoxon signed-rank, paired by "
             "(condition, replicate, sequence length, number of genes).\n" % (len(rows), len(dev), len(hold)))
methods = sorted({m for r in rows for m in r["methods"]})


def mean_table(rs, title):
    lines.append("\n## %s: mean FN rate by (sequence length, #genes)\n" % title)
    keys = sorted({(int(r["sqln"]), r["ngen"]) for r in rs})
    show = [m for m in ["astrid-multi", "astrid-pro", "astrid-pro-min", "ortho-allnodes", "spec-allpairs",
                        "astrid-disco", "astral-pro", "fastmulrfs", "astral-multi(pub)", "mulrf(pub)"] if m in methods]
    lines.append("| sqln | ngen | n | " + " | ".join(show) + " |")
    lines.append("|---|---|---|" + "---|" * len(show))
    for k in keys:
        sub = [r for r in rs if (int(r["sqln"]), r["ngen"]) == k]
        vals = []
        for m in show:
            v = [fnr(r, m) for r in sub if fnr(r, m) is not None]
            vals.append("%.3f" % np.mean(v) if v else "-")
        lines.append("| %d | %d | %d | " % (k[0], k[1], len(sub)) + " | ".join(vals) + " |")
    allv = []
    for m in show:
        v = [fnr(r, m) for r in rs if fnr(r, m) is not None]
        allv.append("%.4f" % np.mean(v) if v else "-")
    lines.append("| all | all | %d | " % len(rs) + " | ".join(allv) + " |")


mean_table(dev, "Dev replicates 01-03")
dm = np.mean([fnr(r, "astrid-pro") for r in dev]); dmin = np.mean([fnr(r, "astrid-pro-min") for r in dev])
chosen = "astrid-pro" if dm <= dmin else "astrid-pro-min"
lines.append("\nDev choice: mean ASTRID-Pro FN %.4f vs closest-copy %.4f -> **%s** used as ASTRID-Pro on held-out.\n"
             % (dm, dmin, chosen))
mean_table(hold, "Held-out replicates 04-10")


def cmp_table(rs, title):
    lines.append("\n## %s: %s vs baselines\n" % (title, chosen))
    lines.append("| baseline | n | mean FN (Pro) | mean FN (base) | mean diff | W/T/L | p |")
    lines.append("|---|---|---|---|---|---|---|")
    for b in ["astrid-multi", "astrid-disco", "astral-pro", "fastmulrfs", "astral-multi(pub)", "mulrf(pub)",
              "astrid-pro-min" if chosen == "astrid-pro" else "astrid-pro"]:
        s = paired(rs, chosen, b)
        if s:
            lines.append("| %s | %d | %.4f | %.4f | %+.4f | %d/%d/%d | %.2g |" % (
                b, s["n"], s["mean_a"], s["mean_b"], s["diff"], s["W"], s["T"], s["L"], s["p"]))


cmp_table(hold, "Held-out (all strata pooled)")
for ng in [25, 100, 500]:
    cmp_table([r for r in hold if r["ngen"] == ng], "Held-out, %d genes" % ng)
for sq in ["25", "100", "250"]:
    cmp_table([r for r in hold if r["sqln"] == sq], "Held-out, %s bp" % sq)
# runtime
secs = collections.defaultdict(list)
for r in rows:
    for k, v in r["secs"].items():
        secs[(k, r["ngen"])].append(v)
lines.append("\n## Runtime (seconds per run, mean; 100 species; single core)\n")
lines.append("| step | 25 genes | 100 genes | 500 genes |")
lines.append("|---|---|---|---|")
for k in sorted({k for k, _ in secs}):
    lines.append("| %s | " % k + " | ".join("%.2f" % np.mean(secs[(k, n)]) if secs[(k, n)] else "-" for n in [25, 100, 500]) + " |")
open(out, "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
