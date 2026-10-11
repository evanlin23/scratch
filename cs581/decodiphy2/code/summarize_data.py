"""Summarize results/data_runs.jsonl -> results/data_summary.md"""
import json, sys, collections
import numpy as np
from scipy.stats import wilcoxon
rows = [json.loads(l) for l in open(sys.argv[1])]
rows = [r for r in rows if "error" not in r and r.get("class_ok")]
out = []
P = out.append
P(f"Runs analysed: {len(rows)} (authors' stored E1 runs, 12 trees, k in {{2,3,5,7,10}}, seeds 1-10, noise0 + 6 noisy settings)\n")
def pct(x): return f"{100*np.mean(x):.1f}%" if len(x) else "-"
P("## 1. How often are configurations non-identifiable?\n")
P("| noise | k | runs | TRUE: adjacent | TRUE: kernel (claw etc.) | FIT: adjacent | FIT: kernel | FIT: boundary x | FIT: class width > 0 | FIT: any flag |")
P("|---|---|---|---|---|---|---|---|---|---|")
groups = collections.defaultdict(list)
for r in rows:
    groups[(r["noise"], r["k"])].append(r); groups[(r["noise"], "all")].append(r); groups[("all", "all")].append(r)
for key in sorted(groups, key=lambda t: (str(t[0]), str(t[1]).zfill(3))):
    g = groups[key]
    anyf = [r["fit_adjacent"] or r["fit_kernel"] or r["fit_boundary"] for r in g]
    P(f"| {key[0]} | {key[1]} | {len(g)} | {pct([r['true_adjacent'] for r in g])} | {pct([r['true_kernel'] for r in g])} | "
      f"{pct([r['fit_adjacent'] for r in g])} | {pct([r['fit_kernel'] for r in g])} | {pct([r['fit_boundary'] for r in g])} | "
      f"{pct([r['nonid_class'] for r in g])} | {pct(anyf)} |")
P("")
fl = [r for r in rows if r["nonid_class"]]
P(f"## 2. Fits whose equivalence class is non-trivial: {len(fl)} of {len(rows)} ({100*len(fl)/len(rows):.1f}%)\n")
w = np.array([r["p_width_max"] for r in fl])
P(f"Max abundance range within the class (max_q hi_q - lo_q): median {np.median(w):.3f}, mean {w.mean():.3f}, 90th pct {np.quantile(w,.9):.3f}, max {w.max():.3f}.")
P(f"Share of flagged fits with range > 0.05: {pct(w > 0.05)}; > 0.10: {pct(w > 0.10)}.\n")
P("### Accuracy: DecoDiPhy's point vs the canonical class representative (same loss, same d)\n")
P("| subset | n | EMD DecoDiPhy | EMD canonical | Wilcoxon p | canonical better / worse | EMD range over class vertices (mean max-min) |")
P("|---|---|---|---|---|---|---|")
def emdrow(name, g):
    if len(g) < 5: return
    a = np.array([r["emd_fit"] for r in g]); b = np.array([r["emd_canon"] for r in g])
    rng = np.array([r["emd_class_max_vertex"] - r["emd_class_min_vertex"] for r in g])
    d = b - a
    try: pv = wilcoxon(d[np.abs(d) > 1e-12]).pvalue
    except ValueError: pv = float("nan")
    P(f"| {name} | {len(g)} | {a.mean():.5f} | {b.mean():.5f} | {pv:.2g} | {(d < -1e-9).sum()}/{(d > 1e-9).sum()} | {rng.mean():.5f} |")
emdrow("all flagged", fl)
for nz in ("noise0", "noise1", "noise2"):
    emdrow(f"flagged, {nz}", [r for r in fl if r["noise"] == nz])
emdrow("flagged, adjacent fits", [r for r in fl if r["fit_adjacent"]])
emdrow("flagged, non-adjacent (kernel) fits", [r for r in fl if not r["fit_adjacent"]])
P("")
se = [r for r in fl if r.get("same_edges")]
P(f"### Flagged fits with exactly the true edge set: {len(se)}\n")
if se:
    a = np.array([r["pL1_fit"] for r in se]); b = np.array([r["pL1_canon"] for r in se]); c = np.array([r["pL1_class_lb"] for r in se])
    P(f"Abundance L1 error: DecoDiPhy {a.mean():.4f}, canonical {b.mean():.4f}, best point in class (lower bound) {c.mean():.4f}.")
    d = b - a
    try: pv = wilcoxon(d[np.abs(d) > 1e-12]).pvalue
    except ValueError: pv = float('nan')
    P(f"canonical better/worse: {(d < -1e-9).sum()}/{(d > 1e-9).sum()}, Wilcoxon p = {pv:.2g}.")
    P(f"True p inside the reported per-placement intervals [lo, hi] (+-0.001): {pct([r['truth_in_p_intervals'] for r in se])}; "
      f"DecoDiPhy's point equal to the true p (+-0.001): {pct([r['fit_p_equal_truth'] for r in se])}.")
    for nz in ("noise0", "noise1", "noise2"):
        g = [r for r in se if r["noise"] == nz]
        if g:
            P(f"- {nz}: n={len(g)}, L1 DecoDiPhy {np.mean([r['pL1_fit'] for r in g]):.4f} / canonical {np.mean([r['pL1_canon'] for r in g]):.4f} / class lb {np.mean([r['pL1_class_lb'] for r in g]):.4f}; "
              f"truth in intervals {pct([r['truth_in_p_intervals'] for r in g])}, point exact {pct([r['fit_p_equal_truth'] for r in g])}")
P("")
P("## 3. Whole-benchmark effect of reporting the canonical representative instead of DecoDiPhy's point\n")
for nz in ("noise0", "noise1", "noise2", None):
    g = [r for r in rows if nz is None or r["noise"] == nz]
    a = np.mean([r["emd_fit"] for r in g]); b = np.mean([r["emd_canon"] for r in g])
    P(f"- {nz or 'all'}: n={len(g)}, mean EMD {a:.5f} -> {b:.5f} ({100*(b-a)/a:+.2f}%)")
open(sys.argv[2], "w").write("\n".join(out) + "\n")
print("\n".join(out))
P2 = []
P2.append("\n## 4. The flag predicts abundance error (runs with exactly the true edge set)\n")
P2.append("| noise | class trivial: n | abundance L1 mean (median) | class non-trivial: n | abundance L1 mean (median) | share of total L1 error in flagged runs |")
P2.append("|---|---|---|---|---|---|")
for nz in ("noise0", "noise1", "noise2"):
    a = [r["pL1_fit"] for r in rows if r["noise"] == nz and r.get("same_edges") and not r["nonid_class"]]
    b = [r["pL1_fit"] for r in rows if r["noise"] == nz and r.get("same_edges") and r["nonid_class"]]
    P2.append(f"| {nz} | {len(a)} | {np.mean(a):.4f} ({np.median(a):.4f}) | {len(b)} | {np.mean(b):.4f} ({np.median(b):.4f}) | {100*sum(b)/(sum(a)+sum(b)):.0f}% |")
open(sys.argv[2], "a").write("\n".join(P2) + "\n")
print("\n".join(P2))
