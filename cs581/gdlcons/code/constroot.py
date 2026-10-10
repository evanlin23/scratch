"""Constant-rate GDL, stock ASTRAL-Pro3 (own min-duplication rooting + overlap tags) on
4-taxon caterpillars from true gene trees: single-sample quartet margins at NFAM families,
next to correct-root overlap tags and the exact prediction for the latter.
Configurations: random (lam, mu = r*lam, branch lengths) with feasible copy numbers,
over-sampling high turnover (lam*t large), the regime of the branch-specific failure cand4.
usage: python constroot.py N NFAM OUT.jsonl [NPROC]"""
import json, os, random, sys
import multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rootscan import job
from predict_margin import predict

LAM = [0.5, 1, 2, 4, 8]
RATIO = [0.5, 0.8, 1.0, 1.25, 2]
TT = [0.05, 0.2, 0.5, 1, 2]
TI = [0.01, 0.05, 0.2, 0.5, 1]

if __name__ == "__main__":
    n, nfam, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    nproc = int(sys.argv[4]) if len(sys.argv) > 4 else 2
    rng = random.Random(77)
    sel = []
    while len(sel) < n:
        lam = rng.choice(LAM); mu = lam * rng.choice(RATIO)
        t = {k: rng.choice(TT) for k in "ABC"}; t.update({k: rng.choice(TI) for k in "xy"})
        g = lam - mu
        if g * (t["y"] + t["x"] + max(t["A"], t["B"], t["C"])) > 3.5 or lam * max(t.values()) > 16:
            continue
        br = {k: (lam, mu, t[k]) for k in "ABCxy"}
        O, H = predict(br, n=80)
        c, w = O + H["AB"], max(H["AC"], H["BC"]); tot = c + H["AC"] + H["BC"]
        if O < 1e-3 or tot <= 0:
            continue
        tD = t["y"] + t["x"] + t["A"]  # D gets the same rates (constant-rate model)
        br["D"] = (lam, mu, tD)
        sel.append({"br": br, "ovl_margin": (c - w) / tot})
    with mp.Pool(nproc) as p:
        for r in p.imap_unordered(job, [(i, r, nfam) for i, r in enumerate(sel)]):
            with open(out, "a") as f:
                f.write(json.dumps(r) + "\n")
            print(r.get("i"), round(r.get("pred_ovl_margin", 0), 3), r.get("ovl"), r.get("own"), r.get("error", ""), flush=True)
