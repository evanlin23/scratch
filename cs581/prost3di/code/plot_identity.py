"""Scatter: per-set SP gain of 3Di L-INS-i over AA L-INS-i vs. estimated identity (results/identity_gain.png)."""
import csv, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
fig, ax = plt.subplots(figsize=(6.4, 4), dpi=150)
for grp, col, marker in (("BB11", "#2a78d6", "o"), ("BBS1", "#2a78d6", "^"), ("BB12", "#eb6834", "s")):
    r = [x for x in rows if x["set"].startswith(grp) and not (grp == "BB11" and x["set"].startswith("BBS"))]
    lab = {"BB11": "RV11 full-length", "BBS1": "RV11 truncated", "BB12": "RV12 full-length"}[grp]
    ax.scatter([float(x["identity"]) for x in r], [float(x["linsi3di"]) - float(x["linsi"]) for x in r],
               s=22, c=col, marker=marker, alpha=0.75, edgecolors="white", linewidths=0.6, label=lab)
ax.axhline(0, color="#888", lw=1)
ax.axvline(0.27, color="#bbb", lw=1, ls="--")
ax.set_xlabel("mean pairwise identity of the L-INS-i alignment (no reference used)")
ax.set_ylabel("SP(3Di L-INS-i) - SP(AA L-INS-i)")
ax.set_title("Predicted-3Di alignment helps below ~25-30% identity", loc="left", fontsize=10)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", color="#eee"); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(sys.argv[2])
