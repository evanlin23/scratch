"""Mechanism check: precision of MAGUS's GCM evidence by edge support (how many of the 10 L-INS-i backbones
contribute to a cross-subset edge), with the reference (diagnostic only).

    python3 edge_support.py REP_DIR [...] > ../results/edge_support.md
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gg  # noqa: E402,F401  (patches bbe)
bbe = gg.bbe

print("| rep | edges | share of edges with support < 4 | share of evidence units on them | precision of units, support 1 / 2 / 3 / ≥ 4 |")
print("|---|---|---|---|---|")
for rep in sys.argv[1:]:
    R = bbe.Rep(rep)
    files, _, _ = bbe.parse_variant(rep, "linsi")
    keys, w, wt, nbb = bbe.graph_edges(R, files)
    lo = nbb < 4
    prec = []
    for sel in (nbb == 1, nbb == 2, nbb == 3, nbb >= 4):
        prec.append("{:.2f}".format(wt[sel].sum() / max(w[sel].sum(), 1)))
    print("| {} | {} | {:.2f} | {:.3f} | {} |".format(os.path.basename(rep), len(keys), lo.mean(), w[lo].sum() / w.sum(),
                                                      " / ".join(prec)), flush=True)
