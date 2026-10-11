# AI-assisted (Claude), exploration code for CS581 project
"""M5 gate selection on training (PREREG addendum 1): threshold rule per signal/direction, chosen by LOO-CV.

    python3 select_gate.py SELECTED REP [REP ...]    -> prints table + JSON rule (fit on all training reps)
"""
import json
import os
import sys

import numpy as np

SIG = ["hardbb_pi", "hardbb_sep", "lowk_wshare", "removed_wshare", "med_neff", "aos", "disp"]


def rows(sel, reps):
    out = []
    for rep in reps:
        r = {json.loads(l)["variant"]: json.loads(l) for l in open(os.path.join(rep, "results2.jsonl"))}
        g = json.load(open(os.path.join(rep, "vote2", "gate.json")))
        g["removed_wshare"] = 1 - json.load(open(os.path.join(rep, "vote2", sel, "model.json")))["kept_weight_frac"]
        out.append((os.path.basename(rep), 100 * (r[sel]["avgErr"] - r["magus"]["avgErr"]), g))
    return out


def fit(d, x):
    """Best (direction, theta): filter iff dir * x >= dir * theta; theta = midpoints (and +-inf). Max mean gain."""
    xs = np.sort(np.unique(x))
    cands = [xs[0] - 1] + list((xs[:-1] + xs[1:]) / 2) + [xs[-1] + 1]
    best = None
    for dr in (1, -1):
        for t in cands:
            on = dr * x >= dr * t
            m = np.where(on, d, 0).mean()
            if best is None or m < best[0] - 1e-12:
                best = (m, dr, t)
    return best


def main():
    sel, reps = sys.argv[1], sys.argv[2:]
    R = rows(sel, reps)
    d = np.array([r[1] for r in R])
    res = {}
    for s in SIG:
        x = np.array([r[2][s] for r in R])
        loo = []
        for i in range(len(R)):
            m = np.arange(len(R)) != i
            _, dr, t = fit(d[m], x[m])
            loo.append(d[i] if dr * x[i] >= dr * t else 0.0)
        full = fit(d, x)
        res[s] = {"loo_mean": float(np.mean(loo)), "train_mean": float(full[0]), "dir": int(full[1]),
                  "theta": float(full[2]), "loo": [round(v, 2) for v in loo]}
    print("| rep | Δ {} vs magus | ".format(sel) + " | ".join(SIG) + " |")
    for n, dd, g in R:
        print("| {} | {:+.2f} | ".format(n, dd) + " | ".join("{:.3f}".format(g[s]) for s in SIG) + " |")
    print("ungated mean {:+.3f}".format(d.mean()))
    for s in SIG:
        print(s, json.dumps(res[s]))
    order = sorted(SIG, key=lambda s: res[s]["loo_mean"])
    best = res[order[0]]["loo_mean"]
    tied = [s for s in order if res[s]["loo_mean"] <= best + 0.05]
    pick = "aos" if "aos" in tied else "disp" if "disp" in tied else order[0]
    print("PICK", pick, json.dumps(res[pick]))


if __name__ == "__main__":
    main()
