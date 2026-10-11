# AI-assisted (Claude), exploration code for CS581 project
"""Fit every vote2 model on a cached feature file (no reference, no scoring): kept edges, time, fit summary."""
import json, sys, time
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__file__))
import vote2
F = vote2.load(sys.argv[1])
print("edges", len(F["w"]), "B", F["X"].shape[1], "n hist", np.bincount(F["nexp"]).tolist())
for v in sys.argv[2:]:
    t = time.time()
    w, info = vote2.weights(F, v)
    info = {k: v_ for k, v_ in info.items() if k not in ("post_hist",)}
    print(v, "kept", int((w > 0).sum()), "wfrac %.3f" % (w.sum() / F["w"].sum()), "%.1fs" % (time.time() - t),
          json.dumps(info, default=float)[:400], flush=True)
