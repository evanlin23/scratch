"""Aggregate scores_*.jsonl -> results TSVs + markdown tables with paired statistics.
Usage: python3 summarize.py RUNROOT RESULTS_DIR
"""
import glob
import json
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import wilcoxon

root, outd = sys.argv[1], sys.argv[2]
os.makedirs(outd, exist_ok=True)
rows = []
for p in glob.glob(f"{root}/*/*/scores_sub*.jsonl") + glob.glob(f"{root}/val_est/*/*/scores_est.jsonl"):
    rows += [json.loads(l) for l in open(p)]
CONDS = ["1000M2", "1000L1", "1000S3", "RNASim1000"]


def paired(a, b):
    """d = b - a (negative = b better). Returns mean diff, W/T/L of b vs a, two-sided Wilcoxon p."""
    a, b = np.asarray(a), np.asarray(b)
    d = b - a
    w, t, l = int((d < -1e-12).sum()), int((abs(d) <= 1e-12).sum()), int((d > 1e-12).sum())
    p = float("nan")
    if (abs(d) > 1e-12).sum() >= 1:
        try:
            p = float(wilcoxon(a, b, zero_method="wilcox").pvalue)
        except ValueError:
            pass
    return d.mean(), w, t, l, p


out = []
# --- full trees
full = [r for r in rows if r["kind"] == "full"]
with open(f"{outd}/full_trees.tsv", "w") as f:
    f.write("cond\trep\tmaxsub\tmethod\tFN\tFP\tseconds\n")
    for r in sorted(full, key=lambda r: (r["cond"], r["rep"], r["maxsub"], r["method"])):
        f.write("%s\t%s\t%d\t%s\t%.5f\t%.5f\t%.1f\n" % (r["cond"], r["rep"], r["maxsub"], r["method"], r["FN"], r["FP"],
                                                     r["seconds"]))
F = defaultdict(dict)  # (cond, label) -> rep -> FN
S = defaultdict(dict)
for r in full:
    lab = r["method"] if r["method"] in ("FastTree", "IQ-TREE") else f"{r['method']}@{r['maxsub']}"
    F[(r["cond"], lab)][r["rep"]] = r["FN"]
    S[(r["cond"], lab)][r["rep"]] = r["seconds"]
labels = sorted({k[1] for k in F}, key=lambda x: (x not in ("FastTree", "IQ-TREE"), x))
out.append("## Full-tree FN error (%), mean over replicates (n = replicates)\n")
out.append("| method | " + " | ".join(CONDS) + " |")
out.append("|---|" + "---|" * len(CONDS))
for lab in labels:
    cells = []
    for c in CONDS:
        v = F.get((c, lab), {})
        cells.append("%.2f (n=%d)" % (100 * np.mean(list(v.values())), len(v)) if v else "–")
    out.append(f"| {lab} | " + " | ".join(cells) + " |")
out.append("\n## Mean wall-clock seconds per replicate (1 thread; GTM rows include the FastTree guide + all subset trees + merge)\n")
out.append("| method | " + " | ".join(CONDS) + " |")
out.append("|---|" + "---|" * len(CONDS))
for lab in labels:
    cells = []
    for c in CONDS:
        v = S.get((c, lab), {})
        cells.append("%.0f" % np.mean(list(v.values())) if v else "–")
    out.append(f"| {lab} | " + " | ".join(cells) + " |")


def cmp_table(title, pairs, D):
    out.append(f"\n## {title}\n")
    out.append("Difference = B − A in FN percentage points (negative = B better). W/T/L = replicates where B is better/tied/worse. Two-sided Wilcoxon signed-rank p (pooled row pools all conditions).\n")
    out.append("| A | B | condition | n | mean A | mean B | mean diff | W/T/L | p |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    for A, B in pairs:
        pa, pb = [], []
        for c in CONDS + ["pooled"]:
            if c == "pooled":
                a, b = pa, pb
            else:
                reps = sorted(set(D.get((c, A), {})) & set(D.get((c, B), {})))
                if not reps:
                    continue
                a = [D[(c, A)][r] for r in reps]
                b = [D[(c, B)][r] for r in reps]
                pa += a
                pb += b
            if not a:
                continue
            md, w, t, l, p = paired(a, b)
            out.append("| %s | %s | %s | %d | %.2f | %.2f | %+.2f | %d/%d/%d | %.3g |" %
                       (A, B, c, len(a), 100 * np.mean(a), 100 * np.mean(b), 100 * md, w, t, l, p))


ms = sorted({r["maxsub"] for r in full})
pairs = []
for m in ms:
    for dl in ("PF", "NNJ"):
        for base in ("FT", "IQ", "BME"):
            pairs.append((f"GTM+{base}@{m}", f"GTM+{dl}@{m}"))
        pairs.append(("FastTree", f"GTM+{dl}@{m}"))
        pairs.append(("IQ-TREE", f"GTM+{dl}@{m}"))
    pairs.append(("FastTree", f"GTM+IQ@{m}"))
    pairs.append(("IQ-TREE", f"GTM+IQ@{m}"))
pairs = [p for p in pairs if any(k[1] == p[0] for k in F) and any(k[1] == p[1] for k in F)]
cmp_table("Paired comparisons, full trees", pairs, F)

# --- subset trees (true and estimated alignment): per-subset pairing
for kind, title in (("subset", "Subset trees on the TRUE alignment (pipeline subsets, <=50 taxa)"),
                    ("subset_est", "Subset trees on a MAFFT L-INS-i ESTIMATED alignment (2 largest true-tree centroid subsets per replicate)")):
    sub = [r for r in rows if r["kind"] == kind and r.get("maxsub", 50) == 50]
    if not sub:
        continue
    D = defaultdict(dict)
    T = defaultdict(list)
    for r in sub:
        D[(r["cond"], r["method"])][(r["rep"], r["subset"])] = r["FN"]
        T[r["method"]].append(r["seconds"])
    meths = sorted({r["method"] for r in sub})
    out.append(f"\n## {title}: mean FN (%) over subsets\n")
    out.append("| method | " + " | ".join(CONDS) + " | mean s/subset |")
    out.append("|---|" + "---|" * (len(CONDS) + 1))
    for m in meths:
        cells = []
        for c in CONDS:
            v = D.get((c, m), {})
            cells.append("%.2f (n=%d)" % (100 * np.mean(list(v.values())), len(v)) if v else "–")
        out.append(f"| {m} | " + " | ".join(cells) + " | %.1f |" % np.mean(T[m]))
    pp = [(b, dl) for dl in ("PF", "NNJ") for b in ("FT", "IQ", "BME", "NJ") if b in meths and dl in meths]
    cmp_table(f"Paired, {title}", pp, D)

open(f"{outd}/SUMMARY.md", "w").write("\n".join(out) + "\n")
print("\n".join(out))
