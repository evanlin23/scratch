# AI-assisted (Claude), exploration code for CS581 project
"""M5 dataset-level gate: reference-free per-replicate signals (PREREG addendum 1) -> REP/vote2/gate.json.

    python3 gate.py REP [REP ...]

Signals: hardbb_pi, hardbb_sep (p1 - p0), lowk_wshare (weight share on k < 4 edges), med_neff, aos, disp.
(removed_wshare of the selected variant is read from its model.json by the selection script.)
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vote2  # noqa: E402


def signals(F):
    S, E, w = F["S"], F["E"], F["w"].astype(float)
    pbb, fit = vote2.m_hardbb(F)
    out = {"hardbb_pi": fit["pi"], "hardbb_sep": fit["p1"] - fit["p0"],
           "lowk_wshare": float(w[F["k"] < 4].sum() / w.sum()), "med_neff": float(np.median(vote2.neff(F)))}
    B = S.shape[1]
    ov = []
    for i in range(B):
        for j in range(i + 1, B):
            both = E[:, i] & E[:, j]
            mx = np.maximum(S[both, i], S[both, j]).sum()
            if mx > 0:
                ov.append(np.minimum(S[both, i], S[both, j]).sum() / mx)
    out["aos"] = float(np.mean(ov))
    n = F["nexp"]
    m = S.sum(1) / np.maximum(n, 1)
    var = (np.where(E, (S - m[:, None]) ** 2, 0).sum(1)) / np.maximum(n, 1)
    sel = n >= 2
    out["disp"] = float((w[sel] * var[sel]).sum() / w[sel].sum())
    return out


if __name__ == "__main__":
    for rep in sys.argv[1:]:
        d = os.path.join(rep, "vote2")
        feats = [f for f in os.listdir(d) if f.startswith("features_") and f.endswith(".npz") and ".tmp" not in f]
        F = vote2.load(os.path.join(d, feats[0]))
        s = signals(F)
        json.dump(s, open(os.path.join(d, "gate.json"), "w"), indent=1)
        print(os.path.basename(rep), json.dumps(s), flush=True)
