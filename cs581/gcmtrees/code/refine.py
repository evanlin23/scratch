"""Diagnostic alignments that separate the two kinds of alignment error (why-not analysis).

    python3 refine.py TRUE EST OUT_SPLIT [OUT_MERGE]

OUT_SPLIT = common refinement of the true and the estimated alignment: every residue goes to column
(true column, estimated column). Its homology pairs are exactly the estimate's true-positive pairs, so it has
the estimate's SPFN (it splits true columns the way the estimate does) and SPFP = 0. Columns are ordered by
true column, then estimated column, which keeps every sequence's residue order.
"""
import sys

import numpy as np

sys.path.insert(0, "/home/user/scratch/cs581/code")
from gcmx import fasta  # noqa: E402

true = fasta.upper(fasta.read(sys.argv[1]))
est = fasta.upper(fasta.read(sys.argv[2]))
taxa = list(true)
cols = {}
keyed = {}
for t in taxa:
    ts, es = true[t], est[t]
    tc = [i for i, c in enumerate(ts) if c != "-"]
    ec = [i for i, c in enumerate(es) if c != "-"]
    res = [c for c in ts if c != "-"]
    assert len(tc) == len(ec), t
    keyed[t] = [((a, b), r) for a, b, r in zip(tc, ec, res)]
    for a, b in zip(tc, ec):
        cols[(a, b)] = 1
order = {k: i for i, k in enumerate(sorted(cols))}
out = {}
for t in taxa:
    row = np.full(len(order), "-", dtype="<U1")
    for k, r in keyed[t]:
        row[order[k]] = r
    out[t] = "".join(row)
fasta.write(out, sys.argv[3])
print(len(order), "columns")
