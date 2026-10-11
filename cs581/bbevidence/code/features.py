"""Per-replicate dataset features (reference-based and reference-free) for predicting when
alternative GCM evidence helps.  python3 features.py REP_DIR [...]  -> JSON lines"""
import json
import os
import random
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bbe  # noqa: E402
from gcmx import fasta  # noqa: E402

for rep in sys.argv[1:]:
    ref = fasta.read(os.path.join(rep, "true.fasta"))
    rows = list(ref.values())
    arr = np.frombuffer("".join(rows).encode(), dtype=np.uint8).reshape(len(rows), -1)
    res = arr != ord("-")
    lens = res.sum(1)
    occ = res.mean(0)
    rng = random.Random(0)
    ids = []
    for _ in range(300):
        i, j = rng.sample(range(len(rows)), 2)
        both = res[i] & res[j]
        if both.sum():
            ids.append((arr[i][both] == arr[j][both]).mean())
    # share of residues in sparse reference columns (< 50% occupancy): unalignable / insertion regions
    sparse = (res * (occ < 0.5)).sum() / res.sum()
    subs = [list(a) for _, a in bbe.subsets(rep)]
    print(json.dumps({"rep": os.path.basename(rep.rstrip("/")), "nseq": len(rows), "mean_len": round(lens.mean(), 1),
                      "len_cv": round(lens.std() / lens.mean(), 3), "ref_cols": arr.shape[1],
                      "pid": round(float(np.mean(ids)), 3), "res_in_sparse_cols": round(float(sparse), 3),
                      "subset_size": round(len(rows) / len(subs), 1)}))
