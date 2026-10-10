"""Pool high-precision batches; limiting matrix D = sum(per-gene avg)/sum(present); delete-one-batch
jackknife SE for every four-point margin (wrong-or-right split minus the next).
Usage: python analyze_hiprec.py TRUE_TREE out.md batches1.jsonl [batches2.jsonl ...]"""
import itertools
import json
import sys

import numpy as np

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from phylo import parse_newick, bipartitions  # noqa: E402

true = parse_newick(sys.argv[1])
recs = [json.loads(l) for f in sys.argv[3:] for l in open(f)]
sp = recs[0]["species"]
idx = {s: i for i, s in enumerate(sp)}
tb = bipartitions(true, idx)
modes = list(recs[0]["acc"].keys())
L = ["# High-precision limiting matrices: %s\n" % sys.argv[1],
     "Batches: %d, families: %d (two independent seeds pooled). Margin = (smallest wrong split sum) - (true split sum): "
     "> 0 means the limiting matrix resolves the quartet correctly. SE: delete-one-batch jackknife.\n" % (
         len(recs), sum(r["nfam"] for r in recs))]
full = (1 << len(sp)) - 1
for m in modes:
    S = np.array([r["acc"][m][0] for r in recs]); C = np.array([r["acc"][m][1] for r in recs])
    def margins(Ssum, Csum):
        D = Ssum / np.maximum(Csum, 1)
        out = {}
        for q in itertools.combinations(range(len(sp)), 4):
            a, b, c, d = q
            splits = {(a, b): D[a, b] + D[c, d], (a, c): D[a, c] + D[b, d], (a, d): D[a, d] + D[b, c]}
            # true split of this quartet
            tq = None
            for (x, y), _ in splits.items():
                m_ = (1 << x) | (1 << y)
                for bp in tb:
                    if (bp & ((1 << a) | (1 << b) | (1 << c) | (1 << d))) in (m_, ((1 << a) | (1 << b) | (1 << c) | (1 << d)) ^ m_):
                        tq = (x, y)
            if tq is None:
                continue
            wrong = min(v for k, v in splits.items() if k != tq)
            out["".join(sp[i] for i in q)] = wrong - splits[tq]
        return out
    tot = margins(S.sum(0), C.sum(0))
    B = len(recs)
    jk = [margins(S.sum(0) - S[i], C.sum(0) - C[i]) for i in range(B)]
    L.append("\n## %s\n" % m)
    L.append("| quartet | margin | jackknife SE | z |")
    L.append("|---|---|---|---|")
    for q, v in tot.items():
        vals = np.array([j[q] for j in jk])
        se = np.sqrt((B - 1) / B * ((vals - vals.mean()) ** 2).sum())
        L.append("| %s | %+.4f | %.4f | %+.1f |" % (q, v, se, v / se if se > 0 else float("nan")))
open(sys.argv[2], "w").write("\n".join(L) + "\n")
print("\n".join(L))
