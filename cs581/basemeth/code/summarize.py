"""Tables for REPORT.md from results/merge.jsonl (+ baselines.jsonl, pasta.jsonl, trees.jsonl if present).

    python3 cs581/basemeth/code/summarize.py cs581/basemeth/results > cs581/basemeth/results/tables.md
"""
import collections
import json
import os
import sys

from scipy.stats import wilcoxon

R = sys.argv[1]
TIE = 0.05  # SP-error points
PROT_CANDIDATES = ["muscle5", "famsa", "probcons", "clustalo"]


def load(name):
    p = os.path.join(R, name)
    return [json.loads(l) for l in open(p)] if os.path.exists(p) else []


def group(ds):
    if ds.startswith("BBA"):
        return "BAliBASE"
    if ds.startswith("HF_"):
        return "HomFam"
    if ds.startswith("SIMMOD"):
        return "SIMMOD"
    if ds.startswith("SIMHIGH"):
        return "SIMHIGH"
    return "nucleotide"


def holm(ps):
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    adj, run = [0] * len(ps), 0
    for rank, i in enumerate(order):
        run = max(run, min(1, (len(ps) - rank) * ps[i]))
        adj[i] = run
    return adj


def pval(d):
    d = [x for x in d if x != 0]
    if len(d) < 2:
        return float("nan")
    return wilcoxon(d).pvalue


def wtl(d):
    return "{}/{}/{}".format(sum(x <= -TIE for x in d), sum(abs(x) < TIE for x in d), sum(x >= TIE for x in d))


rows = load("merge.jsonl")
E = collections.defaultdict(dict)  # dataset -> method -> row
for r in rows:
    E[r["dataset"]][r["method"]] = r
methods = ["linsi", "ginsi", "muscle5", "famsa", "probcons", "clustalo", "kalign", "prank"]
prot = [d for d in E if group(d) != "nucleotide"]
nuc = [d for d in E if group(d) == "nucleotide"]
order = ["BAliBASE", "HomFam", "SIMMOD", "SIMHIGH", "nucleotide"]
ds_sorted = sorted(E, key=lambda d: (order.index(group(d)), d))


def pct(x):
    return "{:.2f}".format(100 * x)


print("## Final MAGUS alignment SP error (%) by subset aligner (merge-only, same subsets & backbones)\n")
print("| dataset | " + " | ".join(methods) + " |")
print("|---|" + "---|" * len(methods))
for d in ds_sorted:
    cells = []
    best = min(E[d][m]["avgErr"] for m in E[d])
    for m in methods:
        if m in E[d]:
            v = pct(E[d][m]["avgErr"])
            cells.append("**" + v + "**" if E[d][m]["avgErr"] == best else v)
        else:
            cells.append("")
    print("| {} | {} |".format(d, " | ".join(cells)))

print("\n## Paired differences vs L-INS-i subsets (SP-error points, negative = better than MAGUS default)\n")
print("| method | set | n | mean d | median d | W/T/L | p (Wilcoxon) | d SPFN | d SPFP | d TC |")
print("|---|---|---|---|---|---|---|---|---|---|")
primary = {}
for m in methods[1:]:
    for setname, dsl in [("all protein", prot)] + [(g, [d for d in E if group(d) == g]) for g in order]:
        pairs = [d for d in dsl if m in E[d] and "linsi" in E[d]]
        if not pairs:
            continue
        diff = [100 * (E[d][m]["avgErr"] - E[d]["linsi"]["avgErr"]) for d in pairs]
        fn = [100 * (E[d][m]["SPFN"] - E[d]["linsi"]["SPFN"]) for d in pairs]
        fp = [100 * (E[d][m]["SPFP"] - E[d]["linsi"]["SPFP"]) for d in pairs]
        tc = [100 * (E[d][m]["TC"] - E[d]["linsi"]["TC"]) for d in pairs]
        p = pval(diff)
        if setname == "all protein" and m in PROT_CANDIDATES:
            primary[m] = (sum(diff) / len(diff), p, len(diff), wtl(diff))
        med = sorted(diff)[len(diff) // 2] if len(diff) % 2 else sum(sorted(diff)[len(diff) // 2 - 1:len(diff) // 2 + 1]) / 2
        print("| {} | {} | {} | {:+.2f} | {:+.2f} | {} | {:.3g} | {:+.2f} | {:+.2f} | {:+.1f} |".format(
            m, setname, len(diff), sum(diff) / len(diff), med, wtl(diff), p,
            sum(fn) / len(fn), sum(fp) / len(fp), sum(tc) / len(tc)))

if primary:
    ms = list(primary)
    adj = holm([primary[m][1] for m in ms])
    print("\n## Pre-registered primary test (protein sets, Holm over {} candidates)\n".format(len(ms)))
    print("| candidate | n | mean d vs L-INS-i | W/T/L | raw p | Holm p |")
    print("|---|---|---|---|---|---|")
    for m, a in sorted(zip(ms, adj), key=lambda x: primary[x[0]][0]):
        md, p, n, w = primary[m]
        print("| {} | {} | {:+.2f} | {} | {:.3g} | {:.3g} |".format(m, n, md, w, p, a))

print("\n## Subset-level accuracy (pooled over subsets; HomFam only subsets with >= 2 seed sequences)\n")
print("| set | " + " | ".join(methods) + " |")
print("|---|" + "---|" * len(methods))
for g in order:
    gd = [d for d in E if group(d) == g]
    gm = {m for d in gd for m in E[d]}
    common = [d for d in gd if gm <= set(E[d])]  # datasets with every method of this group
    cells = []
    for m in methods:
        rs = [E[d][m] for d in common if m in E[d] and E[d][m].get("sub_SPFN") is not None]
        cells.append("{:.2f}".format(100 * sum((r["sub_SPFN"] + r["sub_SPFP"]) / 2 for r in rs) / len(rs)) if rs else "")
    print("| {} (n={}) | {} |".format(g, len(common), " | ".join(cells)))

print("\n## Subset-alignment cost: total CPU seconds over the 25 subsets (1 thread each); ratio to L-INS-i\n")
print("| set | " + " | ".join(methods) + " |")
print("|---|" + "---|" * len(methods))
for g in order:
    cells = []
    for m in methods:
        dl = [d for d in E if group(d) == g and m in E[d] and "linsi" in E[d]]
        if not dl:
            cells.append("")
            continue
        c = sum(E[d][m]["sub_cpu"] for d in dl) / len(dl)
        cl = sum(E[d]["linsi"]["sub_cpu"] for d in dl) / len(dl)
        cells.append("{:.0f} ({:.2g}x)".format(c, c / cl))
    print("| {} | {} |".format(g, " | ".join(cells)))
mw = [E[d]["linsi"]["merge_wall"] for d in E if "linsi" in E[d]]
print("\nGCM merge wall time (4 threads), L-INS-i subsets: median {:.0f} s, range {:.0f}-{:.0f} s.".format(
    sorted(mw)[len(mw) // 2], min(mw), max(mw)))

# control: merge-linsi vs the MAGUS run that produced the subsets
print("\n## Control: GCM merge on re-aligned L-INS-i subsets vs the original MAGUS run (SP error %)\n")
print("| dataset | MAGUS (original) | merge-linsi |")
print("|---|---|---|")
for d in ds_sorted:
    pj = os.path.join("/opt/work/basemeth", d, "prep.json")
    if "linsi" in E[d] and os.path.exists(pj):
        pr = json.load(open(pj))
        print("| {} | {} | {} |".format(d, pct(pr["magus_avgErr"]), pct(E[d]["linsi"]["avgErr"])))

B = load("baselines.jsonl") + load("pasta.jsonl")
if B:
    print("\n## Whole-dataset baselines and PASTA variants (SP error %, wall s, 4 threads)\n")
    bm = sorted({r["method"] for r in B})
    print("| dataset | MAGUS (L-INS-i subsets) | best MAGUS variant | " + " | ".join(bm) + " |")
    print("|---|---|---|" + "---|" * len(bm))
    BB = collections.defaultdict(dict)
    for r in B:
        BB[r["dataset"]][r["method"]] = r
    for d in sorted(BB, key=lambda d: (order.index(group(d)), d)):
        cells = []
        for m in bm:
            r = BB[d].get(m)
            if r is None:
                cells.append("")
            elif "error" in r:
                cells.append("fail/timeout")
            else:
                cells.append("{} ({:.0f} s)".format(pct(r["avgErr"]), r["wall"]))
        lin = pct(E[d]["linsi"]["avgErr"]) if "linsi" in E.get(d, {}) else ""
        bestm = min(E[d], key=lambda m: E[d][m]["avgErr"]) if d in E else None
        best = "{} ({})".format(pct(E[d][bestm]["avgErr"]), bestm) if bestm else ""
        print("| {} | {} | {} | {} |".format(d, lin, best, " | ".join(cells)))

T = load("trees.jsonl")
if T:
    print("\n## Trees: FastTree -lg -gamma, normalized RF (%) vs the true tree\n")
    tm = sorted({r["method"] for r in T})
    print("| dataset | " + " | ".join(tm) + " |")
    print("|---|" + "---|" * len(tm))
    TT = collections.defaultdict(dict)
    for r in T:
        TT[r["dataset"]][r["method"]] = r["RF"]
    for d in sorted(TT):
        print("| {} | {} |".format(d, " | ".join("{:.2f}".format(100 * TT[d][m]) if m in TT[d] else "" for m in tm)))
