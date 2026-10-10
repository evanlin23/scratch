"""Tables for REPORT.md from results/meth_*.jsonl and results/est.jsonl."""
import glob, json, os
from collections import defaultdict
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
ORDER = ["true", "ovl", "own", "apro_bin", "recon_true", "recon_first", "disco_astral", "wqfm_gdl",
         "astral_multi", "gtp_dup", "gtp_dl", "duploss2"]
print("## wrong / blocks (K families per block)\n")
pools = sorted(glob.glob(os.path.join(R, "meth_*.jsonl")))
tab = defaultdict(dict)
Ks = set()
for p in pools:
    name = os.path.basename(p)[5:-6]
    acc = defaultdict(lambda: [0, 0])
    for r in map(json.loads, open(p)):
        acc[(r["method"], r["K"])][0] += r["wrong"]; acc[(r["method"], r["K"])][1] += 1
        Ks.add(r["K"])
    for (m, K), (w, n) in acc.items():
        tab[(m, K)][name] = "%d/%d" % (w, n)
names = [os.path.basename(p)[5:-6] for p in pools]
print("| method | K | " + " | ".join(names) + " |")
print("|---|---|" + "---|" * len(names))
for m in ORDER:
    for K in sorted(Ks):
        if (m, K) in tab:
            print("| %s | %d | " % (m, K) + " | ".join(tab[(m, K)].get(n, "-") for n in names) + " |")
e = os.path.join(R, "est.jsonl")
if os.path.exists(e):
    print("\n## estimated gene trees: mean FN (bipartitions missed) over reps\n")
    acc = defaultdict(list)
    for r in map(json.loads, open(e)):
        acc[(r["setting"], r["kind"], r["method"])].append(r["fn"])
    for k in sorted(acc):
        print("| %s | %s | %s | %s | n=%d |" % (k[0], k[1], k[2], " ".join(map(str, acc[k])), len(acc[k])))
