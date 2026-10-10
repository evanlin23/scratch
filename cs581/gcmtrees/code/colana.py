"""Why-not analysis: where (in which true columns) do the alignment changes land?

    python3 colana.py REP_DIR [...] >> ../results/colana.jsonl

For each true-alignment column c (n_c residues): true pairs C(n_c, 2), and for each estimate the number of
those pairs it recovers (sum over estimated columns of C(k, 2)). Columns are binned by how much they can
inform a tree: 'small' (n_c < 4 residues), 'const' (one amino acid), 'singleton' (variable but not
parsimony-informative), 'informative' (>= 2 states each in >= 2 sequences). Reports per bin: share of true
pairs, recall of MAGUS, and the gain of each method in recovered true pairs (as % of all true pairs), plus
SPFP-side counts (false pairs per estimate). Also the gain restricted to the gappy region (n_c < 50% of taxa).
"""
import json
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, "/home/user/scratch/cs581/code")
from gcmx import fasta  # noqa: E402

V = {"magus": "linsi", "recipe": "wsoft0.03_c_linsi_i_fftns2_es_4", "es3": "linsi_es_3", "hard": "linsi_i_fftns2-op3"}


def colids(aln, taxa):
    """Per taxon: array of column indices of its residues."""
    out = {}
    for t in taxa:
        s = np.frombuffer(aln[t].encode(), dtype=np.uint8)
        out[t] = np.nonzero(s != ord("-"))[0]
    return out


def c2(x):
    x = np.asarray(x, dtype=np.int64)
    return x * (x - 1) // 2


for rep in sys.argv[1:]:
    true = fasta.upper(fasta.read(os.path.join(rep, "true.fasta")))
    taxa = sorted(true)
    T = colids(true, taxa)
    tc = np.concatenate([T[t] for t in taxa])
    L = len(true[taxa[0]])
    # column classes from the true alignment
    cls = np.empty(L, dtype=object)
    ncol = np.bincount(tc, minlength=L)
    chars = np.array([list(true[t]) for t in taxa])
    for c in range(L):
        col = [x for x in chars[:, c] if x != "-"]
        cnt = Counter(col)
        if len(col) < 4:
            cls[c] = "small"
        elif len(cnt) == 1:
            cls[c] = "const"
        elif sum(1 for v in cnt.values() if v >= 2) >= 2:
            cls[c] = "informative"
        else:
            cls[c] = "singleton"
    gappy = ncol < 0.5 * len(taxa)
    tpairs = c2(ncol)
    tot = tpairs.sum()
    row = {"rep": os.path.basename(rep.rstrip("/")), "true_pairs": int(tot)}
    rec = {}
    for m, v in V.items():
        p = os.path.join(rep, "variants", v, "out.fasta")
        if not os.path.exists(p):
            continue
        est = fasta.upper(fasta.read(p))
        E = colids(est, taxa)
        ec = np.concatenate([E[t] for t in taxa])
        key = tc.astype(np.int64) * (ec.max() + 1) + ec
        u, k = np.unique(key, return_counts=True)
        tp_per_tc = np.bincount(u // (ec.max() + 1), weights=c2(k), minlength=L)
        est_pairs = c2(np.bincount(ec)).sum()
        rec[m] = tp_per_tc
        row[m + "_fp_pairs"] = int(est_pairs - tp_per_tc.sum())
        row[m + "_tp_pairs"] = int(tp_per_tc.sum())
    bins = {}
    for b in ("small", "const", "singleton", "informative"):
        sel = cls == b
        d = {"cols": int(sel.sum()), "share_true_pairs": round(tpairs[sel].sum() / tot, 4)}
        if "magus" in rec:
            d["recall_magus"] = round(rec["magus"][sel].sum() / max(tpairs[sel].sum(), 1), 4)
            for m in rec:
                if m != "magus":
                    d["gain_" + m] = round(100 * (rec[m][sel].sum() - rec["magus"][sel].sum()) / tot, 3)
        bins[b] = d
    row["bins"] = bins
    for lab, sel in (("gappy", gappy), ("dense", ~gappy)):
        d = {"cols": int(sel.sum()), "share_true_pairs": round(tpairs[sel].sum() / tot, 4)}
        if "magus" in rec:
            d["recall_magus"] = round(rec["magus"][sel].sum() / max(tpairs[sel].sum(), 1), 4)
            for m in rec:
                if m != "magus":
                    d["gain_" + m] = round(100 * (rec[m][sel].sum() - rec["magus"][sel].sum()) / tot, 3)
        row[lab] = d
    print(json.dumps(row), flush=True)
