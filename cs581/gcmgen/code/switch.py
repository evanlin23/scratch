"""Variant 4: adaptive switch from the reference-free agreement (see ../SPLIT.md).

    python3 switch.py [RESULTS_DIR] >> ../results/tables.md

For a second aligner X and a filtered variant V, the switched recipe uses V when overlap(linsi, X) >= theta
and ELSE otherwise (ELSE = MAGUS's own evidence 'linsi', or a soft variant). theta is chosen on the
training sets only: every midpoint between adjacent training overlaps (plus 'always' / 'never'), maximising
the mean Δ error on training; the held-out sets are then scored once with that theta.
"""
import collections
import glob
import json
import os
import sys

import numpy as np
from scipy.stats import wilcoxon

R = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
TRAIN = {"BBA0101", "BBA0134", "BBA0067", "BBA0039", "SIMMOD_R1", "SIMHIGH_R1", "1000M2", "1000L1", "1000L2", "16S.M"}
HELD = {"BBA0154", "BBA0190", "BBA0081", "BBA0117", "SIMMOD_R2", "SIMHIGH_R2", "1000L3", "1000M3", "1000S1",
        "1000S2", "RNASim", "1000M4", "1000S3", "1000M2_R1", "1000L1_R1", "RNASim_R1"}
RECIPES = [  # (second aligner, variant if overlap >= theta, variant otherwise)
    ("fftns2", "linsi&fftns2", "linsi"),
    ("fftns2", "wsoft0.01:linsi&fftns2", "linsi"),
    ("fftns2", "linsi&fftns2", "wsoft0.03:linsi&fftns2"),
    ("fftns2", "linsi&fftns2", "wsoft0.03:linsi&fftns2#es4"),
    ("fftns2", "linsi&fftns2#es2", "linsi#es4"),
    ("fftns2-op3", "linsi&fftns2-op3", "linsi"),
    ("fftns2", "linsi|cons0.7", "linsi"),
    ("linsi-rt", "linsi&linsi-rt", "linsi"),
    ("linsi-rt", "wsoft0.01:linsi&linsi-rt", "linsi"),
    ("linsi-sh", "linsi&linsi-sh", "linsi"),
]

rows = collections.defaultdict(dict)
for f in glob.glob(os.path.join(R, "*.results.jsonl")):
    for l in open(f):
        r = json.loads(l)
        rows[r["rep"]][r["variant"]] = 100 * r["avgErr"]
agree = collections.defaultdict(dict)
for f in glob.glob(os.path.join(R, "*.agree.jsonl")):
    for l in open(f):
        r = json.loads(l)
        agree[r["rep"]][r["tool"]] = r["overlap"]


def delta(rep, tool, v_hi, v_lo, theta):
    v = v_hi if agree[rep][tool] >= theta else v_lo
    return rows[rep][v] - rows[rep]["linsi"]


def fmt(ds):
    ds = np.array(ds)
    if not len(ds):
        return "n=0"
    p = wilcoxon(ds).pvalue if len(ds) >= 2 and np.any(ds != 0) else float("nan")
    return "n={} {:+.2f} ({}/{}/{}, p={:.3f})".format(len(ds), ds.mean(), int((ds < -0.05).sum()),
                                                     int((abs(ds) <= 0.05).sum()), int((ds > 0.05).sum()), p)


print("\n## Adaptive switch (variant 4): θ fixed on training, evaluated on held-out\n")
print("| X | if overlap ≥ θ | else | θ (train) | train | held-out protein | held-out DNA/RNA | held-out all |")
print("|---|---|---|---|---|---|---|---|")
for tool, hi, lo in RECIPES:
    ok = lambda r: tool in agree[r] and all(v in rows[r] for v in ("linsi", hi, lo))
    tr = sorted((r for r in TRAIN if ok(r)), key=lambda r: agree[r][tool])
    he = sorted(r for r in HELD if ok(r))
    if len(tr) < 3:
        continue
    ovs = [agree[r][tool] for r in tr]
    cands = [-1.0] + [(a + b) / 2 for a, b in zip(ovs, ovs[1:])] + [2.0]
    best = min(cands, key=lambda t: (np.mean([delta(r, tool, hi, lo, t) for r in tr]), t))
    prot = [delta(r, tool, hi, lo, best) for r in he if r.startswith(("BBA", "SIM"))]
    dna = [delta(r, tool, hi, lo, best) for r in he if not r.startswith(("BBA", "SIM"))]
    print("| {} | `{}` | `{}` | {} | {} | {} | {} | {} |".format(
        tool, hi, lo, "never" if best > 1 else "always" if best < 0 else "{:.3f}".format(best),
        fmt([delta(r, tool, hi, lo, best) for r in tr]), fmt(prot), fmt(dna), fmt(prot + dna)))
