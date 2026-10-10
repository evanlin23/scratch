"""Figure + markdown table for the 4-taxon adversarial pools: fraction of disjoint
datasets with a wrong species tree vs number of families (curve4 + atlas4 runs)."""
import collections
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

R = "/opt/runs/gdlcons"
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MARK = ["o", "s", "^", "D", "v", "P", "X", "*"]
LABEL = {"true:0": "ASTRAL-Pro, true tags", "ovl:0": "ASTRAL-Pro, overlap tags (true root)",
         "own:0": "ASTRAL-Pro3 (own root+tags)", "multi:0": "ASTRAL-multi (astral4)",
         "d2s:0.02": "ASTRAL-Pro, naive D->S q=0.02", "d2s:0.1": "ASTRAL-Pro, naive D->S q=0.1",
         "flip:0.1": "ASTRAL-Pro, naive flip q=0.1", "d2s:0.25": "naive D->S q=0.25", "d2s:1": "naive D->S q=1",
         "flipv:0.1": "ASTRAL-Pro, valid flips q=0.1"}


def load(path, key):
    d = collections.defaultdict(list)
    if os.path.exists(path):
        for l in open(path):
            r = json.loads(l)
            if "FN" in r:
                d[(r[key], r["K"])].append(r["FN"] > 0)
    return d


md = ["# Adversarial 4-taxon pools: fraction of datasets with a wrong species tree", "",
      "True gene trees; disjoint blocks of K families from one pool per configuration.", ""]
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
for ax, cand in zip(axes, ["cand2", "cand1"]):
    d = load("%s/curve4_%s.jsonl" % (R, cand), "model")
    a = load("%s/atlas4_%s.jsonl" % (R, cand), "method")
    series = [(LABEL.get(m, m), {k: v for (mm, k), v in d.items() if mm == m})
              for m in ["true:0", "own:0", "ovl:0", "multi:0", "d2s:0.02", "d2s:0.1"]]
    series += [(m, {k: v for (mm, k), v in a.items() if mm == m})
               for m in ["astrid-multi", "astrid-disco", "astral-disco", "fastmulrfs", "stag"]]
    series = [(n, s) for n, s in series if s]
    md += ["## %s" % cand, "", "| method | " + " | ".join(str(k) for k in [100, 500, 2000, 10000]) + " |",
           "|---|---|---|---|---|"]
    for i, (name, s) in enumerate(series):
        ks = sorted(s)
        ax.plot(ks, [sum(s[k]) / len(s[k]) for k in ks], color=COLORS[i % 8], marker=MARK[i % 8],
                ms=4, lw=1.5, ls="-" if i < 6 else "--", label=name)
        md.append("| %s | %s |" % (name, " | ".join("%d/%d" % (sum(s[k]), len(s[k])) if k in s else ""
                                                      for k in [100, 500, 2000, 10000])))
    md.append("")
    ax.set_xscale("log")
    ax.set_title({"cand2": "cand2: overlap tagging is inconsistent", "cand1": "cand1: only naive flips fail"}[cand], fontsize=9)
    ax.set_xlabel("gene families", fontsize=8)
    ax.grid(alpha=0.25, lw=0.5)
    ax.tick_params(labelsize=7)
axes[0].set_ylabel("fraction of datasets wrong", fontsize=8)
h, l = [], []
for ax in axes:
    for a_, b_ in zip(*ax.get_legend_handles_labels()):
        if b_ not in l:
            h.append(a_)
            l.append(b_)
fig.legend(h, l, loc="lower center", ncol=3, fontsize=7, frameon=False)
fig.tight_layout(rect=(0, 0.22, 1, 1))
fig.savefig(os.path.join(RES, "counterexample_curves.png"), dpi=130)
open(os.path.join(RES, "counterexample_summary.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md))
