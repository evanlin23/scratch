"""Pick the general recipe on the training replicates only (../SPLIT.md), then score it on held-out.

    python3 select_recipe.py [RESULTS_DIR] >> ../results/tables.md

Candidates: every variant run on all 10 training replicates. Criterion: lowest mean Δ error over the
training replicates (protein and DNA/RNA pooled); the worst single-replicate Δ is shown as a safety check.
"""
import collections
import glob
import json
import os
import sys

import numpy as np
from scipy.stats import wilcoxon

R = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
TRAIN = ["BBA0101", "BBA0134", "BBA0067", "BBA0039", "SIMMOD_R1", "SIMHIGH_R1", "1000M2", "1000L1", "1000L2", "16S.M"]
HELD = ["BBA0154", "BBA0190", "SIMMOD_R2", "SIMHIGH_R2", "1000L3", "1000M3", "1000S1", "1000S2", "RNASim"]
prot = lambda r: r.startswith(("BBA", "SIM"))

rows = collections.defaultdict(dict)
for f in glob.glob(os.path.join(R, "*.results.jsonl")):
    for l in open(f):
        r = json.loads(l)
        rows[r["rep"]][r["variant"]] = 100 * r["avgErr"]


def d(rep, v):
    return rows[rep][v] - rows[rep]["linsi"]


def fmt(ds):
    ds = np.array(ds)
    if not len(ds):
        return ""
    p = wilcoxon(ds).pvalue if len(ds) >= 2 and np.any(ds != 0) else float("nan")
    return "{:+.2f} ({}/{}/{}, p={:.3f})".format(ds.mean(), int((ds < -0.05).sum()), int((abs(ds) <= 0.05).sum()),
                                                int((ds > 0.05).sum()), p)


cands = [v for v in set().union(*(rows[r] for r in TRAIN)) if v != "linsi" and all(v in rows[r] for r in TRAIN)]
cands.sort(key=lambda v: np.mean([d(r, v) for r in TRAIN]))
print("\n## Recipe selection on the 10 training replicates (no data-type switch)\n")
print("| variant | train all | train protein (6) | train DNA/RNA (4) | worst train Δ |")
print("|---|---|---|---|---|")
for v in cands:
    print("| `{}` | {} | {} | {} | {:+.2f} |".format(v, fmt([d(r, v) for r in TRAIN]),
                                                    fmt([d(r, v) for r in TRAIN if prot(r)]),
                                                    fmt([d(r, v) for r in TRAIN if not prot(r)]),
                                                    max(d(r, v) for r in TRAIN)))
best = cands[0] if cands else None
print("\nSelected: `{}`\n".format(best))
print("## Held-out replicates (scored after selection)\n")
print("| variant | held-out all | held-out protein | held-out DNA/RNA | worst | per replicate |")
print("|---|---|---|---|---|---|")
for v in cands[:6] + ["linsi&fftns2", "linsi#es3", "wsoft0.03:linsi&fftns2"]:
    hs = [r for r in HELD if r in rows and v in rows[r] and "linsi" in rows[r]]
    if not hs:
        continue
    print("| `{}` | {} | {} | {} | {:+.2f} | {} |".format(
        v, fmt([d(r, v) for r in hs]), fmt([d(r, v) for r in hs if prot(r)]), fmt([d(r, v) for r in hs if not prot(r)]),
        max(d(r, v) for r in hs), ", ".join("{} {:+.2f}".format(r, d(r, v)) for r in hs)))
