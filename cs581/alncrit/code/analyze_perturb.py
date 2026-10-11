"""Tree cost of FN-only alignment error (split/shift perturbations of TRUE) vs estimated alignments.

    python analyze_perturb.py OUT.md PLOT.png
Reads /opt/runs/perturb/{scores,trees}.jsonl, /opt/runs/soft/out/results.jsonl,
/opt/runs/alncrit/trees_soft.jsonl. For each alignment type: mean SPFN, SPFP, compression and
FN excess over the TRUE-alignment tree of the same replicate; slope = excess / SPFN (through origin).
"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

rows = []
ps = {r["key"]: r for r in map(json.loads, open("/opt/runs/perturb/scores.jsonl"))}
for r in map(json.loads, open("/opt/runs/perturb/trees.jsonl")):
    s = ps[r["key"]]
    rep, typ = r["key"].split("/")
    rows.append(dict(rep=rep, type=typ, fn=100 * r["fn_rate"], SPFN=100 * s["SPFN"], SPFP=100 * s["SPFP"],
                     comp=s["Compression"]))
ss = {(r["dataset"], r["variant"]): r for r in map(json.loads, open("/opt/runs/soft/out/results.jsonl")) if "avgErr" in r}
for r in map(json.loads, open("/opt/runs/alncrit/trees_soft.jsonl")):
    if r["method"] != "ft":
        continue
    rep, v = r["key"].split("/")
    s = ss.get((rep, v), {"SPFN": 0, "SPFP": 0, "Compression": 1})
    rows.append(dict(rep=rep, type={"default": "MAGUS", "slow": "MAGUS(Slow)", "true": "TRUE"}.get(v, v),
                     fn=100 * r["fn_rate"], SPFN=100 * s["SPFN"], SPFP=100 * s["SPFP"], comp=s["Compression"]))
d = pd.DataFrame(rows)
k = d.groupby("rep").type.transform("count")
d = d[k == k.max()].copy()
true = d[d.type == "TRUE"].set_index("rep").fn
d["excess"] = d.fn - d.rep.map(true)
L = ["Replicates with all alignment types: %d\n" % d.rep.nunique(),
     "| alignment | SPFN % | SPFP % | compression | tree FN % | FN excess over TRUE (pts) | Wilcoxon p (excess ≠ 0) | excess per SPFN point |",
     "|---|---|---|---|---|---|---|---|"]
order = ["TRUE", "split_0.1", "shift_0.1", "split_0.2", "MAGUS", "MAGUS(Slow)", "slow-soft-m3"]
for t in order:
    g = d[d.type == t]
    if g.empty:
        continue
    p = stats.wilcoxon(g.excess).pvalue if t != "TRUE" else float("nan")
    slope = (g.excess * g.SPFN).sum() / (g.SPFN ** 2).sum() if t != "TRUE" else float("nan")
    L.append("| %s | %.2f | %.3f | %.3f | %.2f | %+.2f | %.3g | %.3f |" % (t, g.SPFN.mean(), g.SPFP.mean(), g.comp.mean(),
                                                                       g.fn.mean(), g.excess.mean(), p, slope))
piv = d.pivot_table(index="rep", columns="type", values="fn")
L.append("")
for a, b in (("split_0.1", "MAGUS"), ("split_0.2", "MAGUS"), ("shift_0.1", "MAGUS"), ("split_0.1", "slow-soft-m3")):
    x = piv[a] - piv[b]
    L.append("- %s − %s: tree FN %+.2f pts, W/T/L %d/%d/%d (tie 0.1), Wilcoxon p = %.3g" % (
        a, b, x.mean(), (x < -0.1).sum(), (x.abs() <= 0.1).sum(), (x > 0.1).sum(), stats.wilcoxon(x).pvalue))
open(sys.argv[1], "w").write("\n".join(L) + "\n")
print("\n".join(L))
fig, ax = plt.subplots(figsize=(6, 4.2))
for t in order[1:]:
    g = d[d.type == t]
    ax.scatter(g.SPFN, g.excess, s=16, label=t)
ax.axhline(0, color="#999", lw=0.5)
ax.set_xlabel("SPFN % of the alignment")
ax.set_ylabel("tree FN excess over TRUE-alignment tree (pts)")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(sys.argv[2], dpi=130)
