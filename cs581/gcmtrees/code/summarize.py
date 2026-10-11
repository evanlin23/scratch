"""Tables for REPORT.md from results/{aln,trees,iqtrees,magus}.jsonl.

    python3 summarize.py > ../results/summary.md
"""
import json
import os
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import wilcoxon, spearmanr

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
VAR = {"linsi": "magus", "wsoft0.03:linsi&fftns2#es4": "recipe", "linsi#es3": "es3", "linsi&fftns2-op3": "hard"}
METH = ["recipe", "es3", "hard"]
TIE = 0.1


DIRS = [R] + sorted(os.path.join(R, "..", d) for d in os.listdir(os.path.join(R, "..")) if d.startswith("results_h"))


def rows(f):
    """Rows of results/F plus the helper machines' results_h*/F."""
    out = []
    for d in DIRS:
        p = os.path.join(d, f)
        if os.path.exists(p) and os.path.getsize(p):
            out += [json.loads(l) for l in open(p) if l.strip()]
    return out


aln = defaultdict(dict)
for r in rows("aln.jsonl"):
    if r["variant"] in VAR:
        aln[r["rep"]][VAR[r["variant"]]] = r
tre = defaultdict(dict)
for r in rows("trees.jsonl"):
    tre[r["dataset"]][r["method"]] = r
iq = defaultdict(dict)
for r in rows("iqtrees.jsonl"):
    iq[r["dataset"]][r["method"]] = r


def key(d):
    lvl, rest = d.split("_R")
    return ({"SIMHIGH": 0, "SIMMOD": 1}.get(lvl, 2), lvl, int(rest.split("_")[0]), rest)


def pval(d):
    d = np.asarray(d)
    if len(d) < 2 or np.all(d == 0):
        return float("nan")
    return wilcoxon(d).pvalue


def wtl(d):
    d = np.round(np.asarray(d), 2)  # RF deltas are multiples of ~0.1 point; avoid float edge effects at the band
    return "{}/{}/{}".format((d < -TIE).sum(), (abs(d) <= TIE).sum(), (d > TIE).sum())


done = sorted([d for d in tre if all(m in tre[d] for m in ["true", "magus"] + METH) and "_d" not in d], key=key)
print("## Per dataset: FastTree nRF (%) to the true tree and SP error\n")
print("| dataset | RF true | RF MAGUS | room | RF recipe | RF es3 | RF hard | Δ recipe | Δ es3 | Δ hard | SP err MAGUS | Δ SP recipe | Δ SPFN / ΔSPFP recipe | Δ SP es3 | Δ SP hard |")
print("|---|" + "---|" * 14)
for d in done:
    t, a = tre[d], aln[d]
    rf = {m: 100 * t[m]["RF"] for m in t}
    er = {m: 100 * a[m]["avgErr"] if a[m]["avgErr"] <= 1 else a[m]["avgErr"] for m in a}
    def sp(m, k):
        return 100 * (a[m][k] - a["magus"][k])
    print("| {} | {:.2f} | {:.2f} | {:+.2f} | {:.2f} | {:.2f} | {:.2f} | {:+.2f} | {:+.2f} | {:+.2f} | {:.2f} | {:+.2f} | {:+.2f} / {:+.2f} | {:+.2f} | {:+.2f} |".format(
        d, rf["true"], rf["magus"], rf["magus"] - rf["true"], rf["recipe"], rf["es3"], rf["hard"],
        rf["recipe"] - rf["magus"], rf["es3"] - rf["magus"], rf["hard"] - rf["magus"],
        50 * (a["magus"]["SPFN"] + a["magus"]["SPFP"]), 50 * (a["recipe"]["SPFN"] + a["recipe"]["SPFP"] - a["magus"]["SPFN"] - a["magus"]["SPFP"]),
        sp("recipe", "SPFN"), sp("recipe", "SPFP"),
        50 * (a["es3"]["SPFN"] + a["es3"]["SPFP"] - a["magus"]["SPFN"] - a["magus"]["SPFP"]),
        50 * (a["hard"]["SPFN"] + a["hard"]["SPFP"] - a["magus"]["SPFN"] - a["magus"]["SPFP"])))

print("\n## Paired summaries (method − MAGUS; RF in nRF points, SP in error points; W/T/L tie band 0.1 RF point)\n")
print("| subset | method | n | mean Δ RF | median Δ RF | W/T/L | Wilcoxon p | mean Δ FN | mean Δ FP | mean Δ SP err | mean ΔSPFN | mean ΔSPFP | mean room |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
subsets = [("all", done), ("SIMHIGH", [d for d in done if d.startswith("SIMHIGH")]),
           ("SIMMOD", [d for d in done if d.startswith("SIMMOD")])]
for lab, ds in subsets:
    if not ds:
        continue
    for m in METH:
        drf = [100 * (tre[d][m]["RF"] - tre[d]["magus"]["RF"]) for d in ds]
        dfn = [100 * (tre[d][m]["FN"] - tre[d]["magus"]["FN"]) for d in ds]
        dfp = [100 * (tre[d][m]["FP"] - tre[d]["magus"]["FP"]) for d in ds]
        dsfn = [100 * (aln[d][m]["SPFN"] - aln[d]["magus"]["SPFN"]) for d in ds]
        dsfp = [100 * (aln[d][m]["SPFP"] - aln[d]["magus"]["SPFP"]) for d in ds]
        room = [100 * (tre[d]["magus"]["RF"] - tre[d]["true"]["RF"]) for d in ds]
        print("| {} | {} | {} | {:+.2f} | {:+.2f} | {} | {:.3g} | {:+.2f} | {:+.2f} | {:+.2f} | {:+.2f} | {:+.2f} | {:+.2f} |".format(
            lab, m, len(ds), np.mean(drf), np.median(drf), wtl(drf), pval(drf), np.mean(dfn), np.mean(dfp),
            (np.mean(dsfn) + np.mean(dsfp)) / 2, np.mean(dsfn), np.mean(dsfp), np.mean(room)))
    # true alignment vs MAGUS (how much room)
    drf = [100 * (tre[d]["true"]["RF"] - tre[d]["magus"]["RF"]) for d in ds]
    print("| {} | true aln | {} | {:+.2f} | {:+.2f} | {} | {:.3g} | | | | | | |".format(lab, len(ds), np.mean(drf), np.median(drf), wtl(drf), pval(drf)))

# does the tree change track SP change?
print("\n## Does Δ RF track Δ SPFN or Δ SPFP? (Spearman over dataset × method, n = {})\n".format(3 * len(done)))
x = defaultdict(list)
for d in done:
    for m in METH:
        x["drf"].append(tre[d][m]["RF"] - tre[d]["magus"]["RF"])
        x["dfn"].append(aln[d][m]["SPFN"] - aln[d]["magus"]["SPFN"])
        x["dfp"].append(aln[d][m]["SPFP"] - aln[d]["magus"]["SPFP"])
        x["derr"].append(aln[d][m]["SPFN"] + aln[d][m]["SPFP"] - aln[d]["magus"]["SPFN"] - aln[d]["magus"]["SPFP"])
if len(x["drf"]) > 3:
    for k, lab in (("dfn", "Δ SPFN"), ("dfp", "Δ SPFP"), ("derr", "Δ SP error")):
        rho, p = spearmanr(x["drf"], x[k])
        print("- Δ RF vs {}: rho = {:+.2f} (p = {:.3g})".format(lab, rho, p))
    # across datasets: RF(MAGUS) - RF(true) vs MAGUS SP error
    room = [tre[d]["magus"]["RF"] - tre[d]["true"]["RF"] for d in done]
    for k, lab in (("SPFN", "MAGUS SPFN"), ("SPFP", "MAGUS SPFP")):
        rho, p = spearmanr(room, [aln[d]["magus"][k] for d in done])
        print("- room (RF MAGUS − RF true) vs {}: rho = {:+.2f} (p = {:.3g}, n = {})".format(lab, rho, p, len(done)))

if iq:
    print("\n## IQ-TREE 3 (-m LG+G4 --fast) nRF (%)\n")
    ms = ["true", "magus", "recipe", "es3", "hard"]
    print("| dataset | " + " | ".join(ms) + " | Δ recipe | Δ es3 | Δ hard |")
    print("|---|" + "---|" * (len(ms) + 3))
    dd = defaultdict(list)
    for d in sorted(iq, key=key):
        r = iq[d]
        cells = ["{:.2f}".format(100 * r[m]["RF"]) if m in r else "" for m in ms]
        deltas = []
        for m in ("recipe", "es3", "hard"):
            if m in r and "magus" in r:
                v = 100 * (r[m]["RF"] - r["magus"]["RF"])
                dd[m].append(v)
                deltas.append("{:+.2f}".format(v))
            else:
                deltas.append("")
        print("| {} | {} | {} |".format(d, " | ".join(cells), " | ".join(deltas)))
    for m, v in dd.items():
        print("\n- IQ-TREE {} − MAGUS: n = {}, mean {:+.2f}, W/T/L {}, p = {:.3g}".format(m, len(v), np.mean(v), wtl(v), pval(v)), end="")
    print()

sp = [d for d in sorted(tre, key=key) if "split_magus" in tre[d] and "split_recipe" in tre[d]]
if sp:
    print("\n## Why-not diagnostic: zero-SPFP refinements (FastTree nRF %)\n")
    print("split(X) = common refinement of the true alignment and X: X's true-positive pairs only (X's SPFN, SPFP = 0).\n")
    print("| dataset | true | split(MAGUS) | MAGUS | split(recipe) | recipe | FP cost MAGUS | FP cost recipe | split(recipe) − split(MAGUS) |")
    print("|---|---|---|---|---|---|---|---|---|")
    a, b, c = [], [], []
    for d in sp:
        r = {m: 100 * tre[d][m]["RF"] for m in tre[d]}
        a.append(r["magus"] - r["split_magus"]); b.append(r["recipe"] - r["split_recipe"]); c.append(r["split_recipe"] - r["split_magus"])
        print("| {} | {:.2f} | {:.2f} | {:.2f} | {:.2f} | {:.2f} | {:+.2f} | {:+.2f} | {:+.2f} |".format(
            d, r["true"], r["split_magus"], r["magus"], r["split_recipe"], r["recipe"], a[-1], b[-1], c[-1]))
    print("| mean | | | | | | {:+.2f} | {:+.2f} | {:+.2f} |".format(np.mean(a), np.mean(b), np.mean(c)))

ca = rows("colana.jsonl")
if ca:
    print("\n## Where the alignment changes land (true-alignment columns)\n")
    print("Gains = recovered true pairs as % of all true pairs. 'dense' = true columns holding >= 50% of taxa; 'informative' = parsimony-informative columns. Residue level: misplaced = not in the largest true-column group of its estimated column; split = outside the largest estimated fragment of its true column (% of residues).\n")
    print("| dataset | pairs in dense cols | MAGUS recall dense / gappy | recipe gain dense / gappy | es3 gain dense / gappy | gain in informative cols (recipe) | misplaced res MAGUS / recipe / es3 | of which in blocks >= 10 (MAGUS / recipe) | split res MAGUS / recipe / es3 |")
    print("|---|---|---|---|---|---|---|---|---|")
    for r in sorted(ca, key=lambda r: key(r["rep"])):
        g, dn = r["gappy"], r["dense"]
        print("| {} | {:.1%} | {:.1%} / {:.1%} | {:+.2f} / {:+.2f} | {:+.2f} / {:+.2f} | {:+.2f} | {:.1%} / {:.1%} / {:.1%} | {:.1%} / {:.1%} | {:.1%} / {:.1%} / {:.1%} |".format(
            r["rep"], dn["share_true_pairs"], dn["recall_magus"], g["recall_magus"], dn["gain_recipe"], g["gain_recipe"],
            dn["gain_es3"], g["gain_es3"], r["bins"]["informative"]["gain_recipe"],
            r["magus_misplaced_res"], r["recipe_misplaced_res"], r["es3_misplaced_res"],
            r["magus_misplaced_block_res"], r["recipe_misplaced_block_res"],
            r["magus_split_res"], r["recipe_split_res"], r["es3_split_res"]))

mg = rows("magus.jsonl")
if mg:
    print("\n## MAGUS runtime\n")
    for r in mg:
        print("- {} draw {}: MAGUS wall {} s".format(r.get("dataset", r.get("name")), r.get("draw"), (r.get("magus") or {}).get("wall", r.get("magus_wall"))))
