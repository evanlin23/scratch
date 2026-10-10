"""Summarize the identifiability enumeration (output of run_enum.py).

Per instance we classify the true edge set S:
  matching : no two edges of S share an endpoint (stronger than Assumption 1)
  claw     : some internal node v has all three neighbours in V(S) (endpoints of S)
and report: true solution rank-deficient (continuum of (p,x) on the true edges), alternative
edge sets with k' < k, k' = k, k' = k+1.
"""
import json, glob, os, sys, collections
sys.path.insert(0, os.path.dirname(__file__))
from pdd import Tree

D = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else None
rows = collections.defaultdict(collections.Counter)
examples = []
pred = collections.Counter()
for fn in sorted(glob.glob(os.path.join(D, "*.jsonl"))):
    L = [json.loads(l) for l in open(fn)]
    meta = L[0]["meta"]
    t = Tree(meta["n"], [tuple(e) for e in meta["edges"]])
    nbr = {v: set() for v in range(t.N)}
    for a, b, _ in t.edges:
        nbr[a].add(b); nbr[b].add(a)
    for r in L[1:]:
        S = r["S"]; k = r["k"]
        VS = [x for e in S for x in t.edges[e][:2]]
        matching = len(set(VS)) == len(VS)
        claw = any(nbr[v] <= set(VS) for v in range(t.n, t.N))
        key = (r["regime"], k, "matching" if matching else "adjacent", "claw" if claw else "noclaw")
        c = rows[key]
        c["N"] += 1
        c["true_found"] += r["true_found"]
        c["continuum"] += r["true_full_rank"] is False
        alts = r["alts"]
        lt = any(a["kp"] < k for a in alts); eq = any(a["kp"] == k for a in alts); p1 = any(a["kp"] == k + 1 for a in alts)
        c["alt_k'<k"] += lt; c["alt_k'=k"] += eq; c["alt_k'=k+1"] += p1
        c["any_alt_k'<=k"] += (lt or eq)
        if matching and k == 3:
            pred[(r["regime"], claw, lt or eq)] += 1
        if (lt or eq) and len(examples) < 6 and r["regime"] == "generic" and matching:
            examples.append(dict(n=r["n"], newick=meta["newick_unit"], S=[t.edges[e][:2] for e in S],
                                 alts=[[t.edges[e][:2] for e in a["S"]] for a in alts if a["kp"] <= k][:3]))
cols = ["N", "true_found", "continuum", "alt_k'<k", "alt_k'=k", "any_alt_k'<=k", "alt_k'=k+1"]
lines = ["| regime | k | S structure | claw | " + " | ".join(cols) + " |", "|" + "---|" * (4 + len(cols))]
for key in sorted(rows):
    lines.append("| " + " | ".join(map(str, key)) + " | " + " | ".join(str(rows[key][c]) for c in cols) + " |")
lines.append("")
lines.append("k=3, S a matching: (regime, claw, has alternative with k'<=k) -> count")
for k_ in sorted(pred):
    lines.append(f"  {k_}: {pred[k_]}")
lines.append("")
lines.append("examples (generic regime, S a matching, alternative with k'<=k):")
for e in examples:
    lines.append("  " + json.dumps(e))
txt = "\n".join(lines)
print(txt)
if out:
    open(out, "w").write(txt + "\n")
