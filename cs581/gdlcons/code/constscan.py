"""Constant-rate GDL (one lam, one mu on every branch): does correct-root species-overlap
tagging ever make ASTRAL-Pro inconsistent on a 4-taxon caterpillar (((A,B)x,C)y,D)?
Exact limiting scores from predict_margin.predict (D's branch does not matter: only
P(D present) multiplies every score).  Balanced 4-taxon trees with no root branch have
no wrong SQs, so only caterpillars are scanned.
usage: python constscan.py OUT.jsonl"""
import itertools, json, sys
from predict_margin import predict

LAM = [0.1, 0.25, 0.5, 1, 2, 4, 8]
RATIO = [0, 0.25, 0.5, 0.8, 1.0, 1.25, 2, 4]          # mu = ratio * lam  (ratio < 1: supercritical)
T_TIP = [0.05, 0.2, 0.5, 1, 2, 4]
T_INT = [0.01, 0.05, 0.2, 0.5, 1, 2]
out = sys.argv[1]
if len(sys.argv) > 2:
    LAM = [float(sys.argv[2])]
n = bad = 0
with open(out, "w") as f:
    for lam, r, tA, tB, tC, tx, ty in itertools.product(LAM, RATIO, T_TIP, T_TIP, T_TIP, T_INT, T_INT):
        if tB < tA:            # A/B symmetric
            continue
        mu = r * lam
        if (lam - mu) * max(tA + tx, tB + tx, tC) + (lam - mu) * ty > 40:
            continue
        br = {"A": (lam, mu, tA), "B": (lam, mu, tB), "C": (lam, mu, tC), "x": (lam, mu, tx), "y": (lam, mu, ty)}
        try:
            O, H = predict(br, n=60)
        except (OverflowError, ZeroDivisionError, ValueError):
            continue
        c, w = O + H["AB"], max(H["AC"], H["BC"])
        tot = c + H["AC"] + H["BC"]
        if tot <= 0 or O < 1e-12:
            continue
        n += 1
        m = (c - w) / tot
        rec = {"lam": lam, "mu": mu, "tA": tA, "tB": tB, "tC": tC, "tx": tx, "ty": ty, "O": O, "H": H, "margin": m,
               "qstar": O / (w - H["AB"]) if w > H["AB"] else None}
        if m < 0.02 or rec["qstar"] is not None:
            f.write(json.dumps(rec) + "\n")
        if m < 0:
            bad += 1
print("configs", n, "inconsistent (ovl)", bad)
