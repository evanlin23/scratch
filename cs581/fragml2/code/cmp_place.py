"""Compare patched vs stock EPA-ng placements (best-LWR edge per query) for each replicate.

    python cmp_place.py > ../results/place_cmp.md
"""
import glob
import json
import os

import numpy as np

MLDATA = os.environ.get("MLDATA", "/opt/data/mlcache")


def best(path):
    j = json.load(open(path))
    f = j["fields"]
    ie, il, ill = f.index("edge_num"), f.index("like_weight_ratio"), f.index("likelihood")
    out = {}
    for pq in j["placements"]:
        p = max(pq["p"], key=lambda x: x[il])
        for n in pq.get("n", []) + [x[0] for x in pq.get("nm", [])]:
            out[n] = (p[ie], p[il], p[ill])
    return out


print("| dataset | rep | backbone | queries | same best edge | mean best-LWR patched / stock | "
      "share LWR > 0.99 patched / stock | mean (stock − patched) logL of best placement |")
print("|---|---|---|---|---|---|---|---|")
for fx in sorted(glob.glob(os.path.join(MLDATA, "*", "R*", "trees", "cache", "*.place_fix.jplace"))):
    st = fx.replace("place_fix", "place_stock")
    if not os.path.exists(st):
        continue
    a, b = best(fx), best(st)
    q = sorted(set(a) & set(b))
    same = sum(a[n][0] == b[n][0] for n in q)
    parts = fx.split(os.sep)
    print("| %s | %s | %s | %d | %d (%.1f%%) | %.3f / %.3f | %.1f%% / %.1f%% | %+.1f |" % (
        parts[-5], parts[-4], os.path.basename(fx).split(".place")[0].split(".", 1)[1], len(q), same, 100.0 * same / len(q),
        np.mean([a[n][1] for n in q]), np.mean([b[n][1] for n in q]),
        100 * np.mean([a[n][1] > 0.99 for n in q]), 100 * np.mean([b[n][1] > 0.99 for n in q]),
        np.mean([b[n][2] - a[n][2] for n in q])))
