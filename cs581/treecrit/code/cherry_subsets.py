# AI-assisted (Claude), exploration code for CS581 project
"""Share of true-tree cherries (and of random taxon pairs) whose two taxa sit in the same MAGUS subset.
Pairs inside one subset are aligned by the subset alignment, so the GCM merge step cannot change them.

    python3 cherry_subsets.py  -> results/cherry_subsets.json
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from measures import read_fasta, W  # noqa: E402
from measures2 import pairs_for  # noqa: E402

out = {}
for key in sorted(os.listdir(W)):
    rep = os.path.join(W, key)
    sd = os.path.join(rep, "inputs", "subalignments")
    if not os.path.isdir(sd):
        continue
    sub = {}
    for i, f in enumerate(sorted(os.listdir(sd))):
        for t in read_fasta(os.path.join(sd, f)):
            sub[t] = i
    taxa = sorted(read_fasta(os.path.join(rep, "true.fasta")))
    cher, rnd = pairs_for(key, taxa)
    same = lambda ps: float(np.mean([sub[taxa[a]] == sub[taxa[b]] for a, b in ps]))  # noqa: E731
    out[key] = {"cherries": len(cher), "cherry_same_subset": same(cher), "random_same_subset": same(rnd)}
    print(key, out[key], flush=True)
json.dump(out, open(os.path.join(HERE, "..", "results", "cherry_subsets.json"), "w"), indent=1)
v = [o["cherry_same_subset"] for o in out.values()]
print("mean cherry same-subset", np.mean(v), "random", np.mean([o["random_same_subset"] for o in out.values()]))
