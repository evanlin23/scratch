"""results/scaling.png: wall time vs #species and vs #genes (idle machine, one job at a time)."""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
rows = [json.loads(l) for l in open(os.path.join(R, "scaling.jsonl"))]
ok = {}
for r in rows:
    if "error" in r:
        continue
    ok[(r["cond"], r["ngen"], r["method"], r["threads"])] = r["sec"]
style = {("astrid-pro", 1): ("ASTRID-Pro 1t", "#1f77b4", "-"), ("astrid-pro", 4): ("ASTRID-Pro 4t", "#1f77b4", "--"),
         ("astrid-disco", 1): ("ASTRID-DISCO 1t", "#2ca02c", "-"), ("asteroid", 1): ("Asteroid 1t", "#9467bd", "-"),
         ("astral-pro3", 4): ("ASTRAL-Pro3 4t", "#d62728", "-"), ("wqfm-gdl", 4): ("wQFM-GDL", "#ff7f0e", "-")}
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for (m, t), (lab, c, ls) in style.items():
    xs = [k for k in (100, 200, 500, 1000) if (f"taxa_{k}", 1000, m, t) in ok]
    if xs:
        ax[0].plot(xs, [ok[(f"taxa_{k}", 1000, m, t)] for k in xs], ls, color=c, marker="o", label=lab)
    xs = [n for n in (100, 1000, 10000) if ("genes", n, m, t) in ok]
    if xs:
        ax[1].plot(xs, [ok[("genes", n, m, t)] for n in xs], ls, color=c, marker="o", label=lab)
for a, xl in zip(ax, ["#species (1000 genes, DISCO species_1000 subsets)", "#genes (100 species, gtrees_10000_l1)"]):
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel(xl); a.set_ylabel("wall time (s)"); a.grid(alpha=.3)
ax[0].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(R, "scaling.png"), dpi=120)
