"""Tables for REPORT.md from bbtool_bench rows (+ baselines and trees).

    python summarize.py results/bbtool.jsonl [results/baselines.jsonl] [results/trees.jsonl] > results/tables.md

Differences are in error points (x100), MAGUS+Clustal-backbones minus MAGUS; negative = Clustal backbones better.
W/T/L counts Clustal-backbone wins/ties/losses with |d| < 0.05 points a tie; p = two-sided Wilcoxon signed-rank.
Categories come from the dataset-name prefix (10AA_, HF_, SIMMOD_, SIMHIGH_, ...).
"""
import json, sys
from collections import defaultdict
from scipy.stats import wilcoxon

rows = [json.loads(l) for l in open(sys.argv[1])]
base = [json.loads(l) for l in open(sys.argv[2])] if len(sys.argv) > 2 and sys.argv[2] else []
trees = [json.loads(l) for l in open(sys.argv[3])] if len(sys.argv) > 3 and sys.argv[3] else []
cat = lambda n: n.split("_")[0]
pt = lambda x: 100 * x


def stats(ds):
    w = sum(d < -0.05 for d in ds); l = sum(d > 0.05 for d in ds); t = len(ds) - w - l
    p = wilcoxon(ds).pvalue if len(ds) >= 2 and any(abs(d) > 1e-9 for d in ds) else float("nan")
    return "{:+.2f}".format(sum(ds) / len(ds)), "{}/{}/{}".format(w, t, l), "{:.3g}".format(p)


print("## Per dataset (error = (SPFN+SPFP)/2, %)\n")
print("| dataset | n | ref | MAGUS | e2e-Clustal-bb | d e2e | merge L-INS-i bb | merge Clustal bb | d merge | merge mafft-auto bb | merge union | MAGUS wall s | e2e wall s | ratio | bb SPFP L-INS-i / Clustal |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in rows:
    m, e = r["magus"], r.get("e2e-clustalo")
    mm, mc = r.get("merge-mafft"), r.get("merge-clustalo")
    ma, mu = r.get("merge-mafft-auto"), r.get("merge-union-clustalo")
    f = lambda x: "{:.2f}".format(pt(x["avgErr"])) if x else "-"
    print("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {:.0f} | {} | {} | {} |".format(
        r["dataset"], r["nseq"], r.get("nref", r["nseq"]), f(m), f(e),
        "{:+.2f}".format(pt(e["avgErr"] - m["avgErr"])) if e else "-", f(mm), f(mc),
        "{:+.2f}".format(pt(mc["avgErr"] - mm["avgErr"])) if mc and mm else "-", f(ma), f(mu), m["wall"],
        "{:.0f}".format(e["wall"]) if e else "-", "{:.2f}".format(e["wall"] / m["wall"]) if e else "-",
        "{:.1f} / {:.1f}".format(pt(m["bb_SPFP"]), pt(mc["bb_SPFP"])) if mc and m.get("bb_SPFP") is not None and mc.get("bb_SPFP") is not None else "-"))

print("\n## Summary by category (d in error points, Clustal-bb minus L-INS-i-bb)\n")
print("| category | comparison | n | mean d err | W/T/L | p | mean d SPFN | mean d SPFP | mean d TC | wall ratio |")
print("|---|---|---|---|---|---|---|---|---|---|")
groups = defaultdict(list)
for r in rows:
    groups[cat(r["dataset"])].append(r)
groups["ALL"] = rows
for g, rs in groups.items():
    for label, a, b in (("end-to-end", "magus", "e2e-clustalo"), ("merge-only (paired)", "merge-mafft", "merge-clustalo"),
                        ("merge mafft-auto bb", "merge-mafft", "merge-mafft-auto"), ("merge union", "merge-mafft", "merge-union-clustalo")):
        rr = [r for r in rs if a in r and b in r]
        if not rr:
            continue
        d = [pt(r[b]["avgErr"] - r[a]["avgErr"]) for r in rr]
        mean, wtl, p = stats(d)
        dfn = sum(pt(r[b]["SPFN"] - r[a]["SPFN"]) for r in rr) / len(rr)
        dfp = sum(pt(r[b]["SPFP"] - r[a]["SPFP"]) for r in rr) / len(rr)
        tcs = [pt(r[b]["TC"] - r[a]["TC"]) for r in rr if "TC" in r[a] and "TC" in r[b]]
        tc = "{:+.2f}".format(sum(tcs) / len(tcs)) if tcs else "-"
        if b == "e2e-clustalo":
            ratio = "{:.2f}".format(sum(r[b]["wall"] for r in rr) / sum(r[a]["wall"] for r in rr))
        elif "backbone_wall" in rr[0][b]:
            ratio = "bb realign {:.0f} s avg".format(sum(r[b]["backbone_wall"] for r in rr) / len(rr))
        else:
            ratio = "-"
        print("| {} | {} | {} | {} | {} | {} | {:+.2f} | {:+.2f} | {} | {} |".format(g, label, len(rr), mean, wtl, p, dfn, dfp, tc, ratio))

if base:
    print("\n## Baselines (error %, wall s; 4 threads, 30-min cap)\n")
    print("| dataset | MAGUS | L-INS-i | L-INS-i wall | Clustal Omega | Clustal wall |")
    print("|---|---|---|---|---|---|")
    bd = {(b["dataset"], b["method"]): b for b in base}
    for r in rows:
        cell = lambda m: ("{:.2f}".format(pt(bd[(r["dataset"], m)]["avgErr"])) if bd[(r["dataset"], m)].get("wall") else "timeout") if (r["dataset"], m) in bd else "-"
        wall = lambda m: "{:.0f}".format(bd[(r["dataset"], m)]["wall"]) if bd.get((r["dataset"], m), {}).get("wall") else "-"
        print("| {} | {:.2f} | {} | {} | {} | {} |".format(r["dataset"], pt(r["magus"]["avgErr"]), cell("linsi"), wall("linsi"), cell("clustalo"), wall("clustalo")))

if trees:
    print("\n## Trees (FastTree -lg -gamma; RF = normalized Robinson-Foulds vs true tree, %)\n")
    td = defaultdict(dict)
    for t in trees:
        td[t["dataset"]][t["method"]] = t
    ms = ["true", "magus", "e2e-clustalo", "merge-clustalo"]
    print("| dataset | " + " | ".join(ms) + " | d RF e2e | d RF merge |")
    print("|---|" + "---|" * (len(ms) + 2))
    de, dm = defaultdict(list), defaultdict(list)
    for n, t in sorted(td.items()):
        g = lambda m: "{:.2f}".format(pt(t[m]["RF"])) if m in t else "-"
        x = pt(t["e2e-clustalo"]["RF"] - t["magus"]["RF"]) if "e2e-clustalo" in t and "magus" in t else None
        # merge-mafft has exactly magus's SP scores and length (columns only reordered), so magus's tree stands in
        y = pt(t["merge-clustalo"]["RF"] - t["magus"]["RF"]) if "merge-clustalo" in t and "magus" in t else None
        if x is not None: de[cat(n)].append(x)
        if y is not None: dm[cat(n)].append(y)
        print("| {} | {} | {} | {} |".format(n, " | ".join(g(m) for m in ms), "{:+.2f}".format(x) if x is not None else "-", "{:+.2f}".format(y) if y is not None else "-"))
    print("\n| category | comparison | n | mean d RF | W/T/L | p |\n|---|---|---|---|---|---|")
    for g in de:
        print("| {} | end-to-end | {} | {} | {} | {} |".format(g, len(de[g]), *stats(de[g])))
        if dm[g]:
            print("| {} | merge-only | {} | {} | {} | {} |".format(g, len(dm[g]), *stats(dm[g])))
