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
# more than 8 series: the line style (solid = ASTRAL-Pro variants, dashed = atlas) and marker carry identity too
MARK = ["o", "s", "^", "D", "v", "P", "X", "*"]
LABEL = {"true:0": "ASTRAL-Pro, true tags", "ovl:0": "ASTRAL-Pro, overlap tags (true root)",
         "own:0": "ASTRAL-Pro3 (own root+tags)", "multi:0": "ASTRAL-multi (astral4)",
         "d2s:0.02": "ASTRAL-Pro, naive D->S q=0.02", "d2s:0.1": "ASTRAL-Pro, naive D->S q=0.1",
         "flip:0.1": "ASTRAL-Pro, naive flip q=0.1", "d2s:0.25": "naive D->S q=0.25", "d2s:1": "naive D->S q=1",
         "flipv:0.1": "ASTRAL-Pro, valid flips q=0.1",
         "rovl:0.3": "ASTRAL-Pro, random re-root p=0.3", "rovl:1": "ASTRAL-Pro, random re-root p=1",
         "rsp-A:1": "ASTRAL-Pro, root on an A leaf"}


ORDER = [LABEL[m] for m in ["true:0", "own:0", "ovl:0", "multi:0", "rovl:0.3", "rovl:1", "rsp-A:1", "d2s:0.02", "d2s:0.1"]] + \
    ["astrid-multi", "astrid-disco", "astral-disco", "fastmulrfs", "stag"]


def load(path, key):
    d = collections.defaultdict(list)
    if os.path.exists(path):
        for l in open(path):
            r = json.loads(l)
            if "FN" in r:
                d[(r[key], r["K"])].append(r["FN"] > 0)
    return d


ROW1 = ["true:0", "own:0", "ovl:0", "rovl:0.3", "rovl:1", "rsp-A:1", "d2s:0.02", "d2s:0.1"]
ROW2 = ["own:0", "multi:0", "astrid-multi", "astrid-disco", "astral-disco", "fastmulrfs", "stag"]
CANDS = ["cand2", "cand3", "cand4", "cand1"]
TITLE = {"cand2": "cand2: hidden-paralog tagging", "cand3": "cand3: random re-rooting",
         "cand4": "cand4: ASTRAL-Pro's own rooting", "cand1": "cand1: only naive flips fail"}
KS = [100, 500, 1000, 2000, 5000, 10000]
md = ["# Adversarial 4-taxon pools: fraction of datasets with a wrong species tree", "",
      "True gene trees; disjoint blocks of K families from one pool per configuration.", ""]
data = {}
for cand in CANDS:
    d = load("%s/curve4_%s.jsonl" % (R, cand), "model")
    a = load("%s/atlas4_%s.jsonl" % (R, cand), "method")
    ser = {}
    for (m, k), v in list(d.items()) + list(a.items()):
        ser.setdefault(m, {})[k] = v
    data[cand] = ser
    md += ["## %s (%s)" % (cand, TITLE[cand].split(": ")[1]), "", "| method | " + " | ".join(map(str, KS)) + " |",
           "|---|" + "---|" * len(KS)]
    for m in dict.fromkeys(ROW1 + ROW2):
        if m in ser:
            md.append("| %s | %s |" % (LABEL.get(m, m), " | ".join("%d/%d" % (sum(ser[m][k]), len(ser[m][k])) if k in ser[m] else "" for k in KS)))
    md.append("")

fig, axes = plt.subplots(2, 4, figsize=(15, 6.6), sharey=True)
for row, members in enumerate([ROW1, ROW2]):
    for col, cand in enumerate(CANDS):
        ax = axes[row][col]
        for j, m in enumerate(members):
            if m not in data[cand]:
                continue
            ks = sorted(data[cand][m])
            ax.plot(ks, [sum(data[cand][m][k]) / len(data[cand][m][k]) for k in ks], color=COLORS[j], marker=MARK[j],
                    ms=4, lw=1.5, label=LABEL.get(m, m))
        ax.set_xscale("log")
        ax.set_xlim(70, 14000)
        ax.grid(alpha=0.25, lw=0.5)
        ax.tick_params(labelsize=7)
        if row == 0:
            ax.set_title(TITLE[cand], fontsize=9)
        else:
            ax.set_xlabel("gene families", fontsize=8)
        if col == 0:
            ax.set_ylabel(["ASTRAL-Pro rooting/tagging variants", "method atlas"][row] + "\nfraction of datasets wrong", fontsize=8)
    hs = [plt.Line2D([], [], color=COLORS[j], marker=MARK[j], lw=1.5, ms=4) for j in range(len(members))]
    fig.legend(hs, [LABEL.get(m, m) for m in members], loc="center right", bbox_to_anchor=(1.0, 0.75 - 0.47 * row),
               fontsize=7, frameon=False)
fig.tight_layout(rect=(0, 0, 0.8, 1))
fig.savefig(os.path.join(RES, "counterexample_curves.png"), dpi=130)
open(os.path.join(RES, "counterexample_summary.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md))
