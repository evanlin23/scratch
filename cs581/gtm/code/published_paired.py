"""Paired comparisons among the published DTM trees (rescore_published.tsv)."""
import csv, sys, collections
from stats import paired, fmt
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
v = collections.defaultdict(dict)
for r in rows:
    v[(r["condition"], r["method"])][r["replicate"]] = float(r["FN"])
conds = ["RNASim1000", "Cox1-HET", "1000M1-HF"]
for c in conds:
    for g in ("FT", "IQ"):
        for other in ("TreeMerge", "CINC"):
            A, B = v[(c, f"GTM/{g}")], v[(c, f"{other}/{g}")]
            reps = sorted(A)
            print(f"{c:10s} {g}  A=GTM B={other:9s} " + fmt(paired([A[r] for r in reps], [B[r] for r in reps])))
# pooled over all 40 replicate x guide pairs
for other in ("TreeMerge", "CINC"):
    a, b = [], []
    for c in conds:
        for g in ("FT", "IQ"):
            A, B = v[(c, f"GTM/{g}")], v[(c, f"{other}/{g}")]
            for r in sorted(A):
                a.append(A[r]); b.append(B[r])
    print(f"POOLED     A=GTM B={other:9s} " + fmt(paired(a, b)))
