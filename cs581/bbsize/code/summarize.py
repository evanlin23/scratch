# AI-assisted (Claude), exploration code for CS581 project
"""Tables for REPORT.md from results/results.jsonl and backbone timing (python3 summarize.py > results/tables.md)."""
import glob, json, os, collections

R = "/home/user/scratch/cs581/bbsize/results"
rows = [json.loads(l) for l in open(os.path.join(R, "results.jsonl"))]
DS = ["BBA0101", "BBA0067", "SIMHIGH", "SIMMOD", "1000M2", "16S.M"]
CONDS = [(100, 10), (200, 10), (400, 10), (200, 5), (200, 20), (200, 40), (400, 5)]
VAR = ["magus", "frac0.2", "frac0.3", "frac0.4", "frac0.5", "hard-bb"]
by = {(r["dataset"].split("_R")[0], r["size"], r["B"], r["variant"]): r for r in rows}

print("## SP error (%) = (SPFN+SPFP)/2, merge-only on MAGUS's 25 subsets\n")
print("best frac = lowest of frac0.2-0.5 (oracle choice); Δ columns vs magus at the same (size, B).\n")
print("| dataset | size | B | magus | F0.2 | F0.3 | F0.4 | F0.5 | hard-bb | best frac Δ | hard-bb Δ | hard-bb cutoff (k at full n) |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
for d in DS:
    for s, b in CONDS:
        g = {v: by.get((d, s, b, v)) for v in VAR}
        if not g["magus"]:
            continue
        e = {v: (100 * g[v]["avgErr"] if g[v] else None) for v in VAR}
        f = lambda x: "%.2f" % x if x is not None else "–"
        fr = [e[v] for v in VAR[1:5] if e[v] is not None]
        bf = min(fr) - e["magus"] if fr else None
        hb = e["hard-bb"] - e["magus"] if e["hard-bb"] is not None else None
        cut = g["hard-bb"]["cutoff_by_n"].get(str(b)) if g["hard-bb"] else None
        print("| %s | %d | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            d, s, b, f(e["magus"]), f(e["frac0.2"]), f(e["frac0.3"]), f(e["frac0.4"]), f(e["frac0.5"]),
            f(e["hard-bb"]), ("%+.2f" % bf) if bf is not None else "–", ("%+.2f" % hb) if hb is not None else "–", cut))
print("\n## Evidence agreement (MAGUS graph, all cross-subset edges) and beta-binomial fit\n")
print("| dataset | size | B | edges | mean k/n | weighted k/n | unanimous | k=1 | mean n | p1 | p0 | pi |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
for d in DS:
    for s, b in CONDS:
        r = by.get((d, s, b, "magus"))
        if not r or "mean_kn" not in r:
            continue
        fb = r["fit_bb"]
        print("| %s | %d | %d | %d | %.3f | %.3f | %.3f | %.3f | %.2f | %.3f | %.3f | %.3f |" % (
            d, s, b, r["edges"], r["mean_kn"], r["wmean_kn"], r["frac_unanimous"], r["frac_k1"], r["mean_n"],
            fb["p1"], fb["p0"], fb["pi"]))
print("\n## Backbone cost (new backbones, L-INS-i --thread 1, up to 4 at once; CPU = user+sys)\n")
print("| dataset | size | n backbones | seqs/backbone | mean CPU s | mean wall s | max RSS MB |")
print("|---|---|---|---|---|---|---|")
for p in sorted(glob.glob("/opt/work/bbsize/bb/*/s*/timing.jsonl")):
    t = [json.loads(l) for l in open(p)]
    d = p.split("/")[-3].split("_b20")[0].split("_R0_B20")[0].split("_R1")[0]
    s = int(p.split("/")[-2][1:])
    m = lambda k: sum(x[k] for x in t) / len(t)
    print("| %s | %d | %d | %.0f | %.0f | %.0f | %.0f |" % (d, s, len(t), m("nseq"), m("cpu"), m("wall"), max(x["maxrss_mb"] for x in t)))
