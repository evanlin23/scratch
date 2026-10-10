"""How often is ASTRAL-Pro with species-overlap tags (correct roots) inconsistent?
Uses the exact limiting formula of predict_margin.py on random rate configurations of
the 4-taxon caterpillar (((A,B)x,C)y,D) with a single-copy outgroup D.
usage: python prevalence.py N OUT.jsonl"""
import json
import random
import sys

from predict_margin import predict

LAM = [0, 0.5, 1, 2, 4, 8]
MU = [0, 0.5, 1, 2, 4, 8]
TS = [0.05, 0.2, 0.5, 1.0, 2.0, 4.0]

if __name__ == "__main__":
    n, out = int(sys.argv[1]), sys.argv[2]
    rng = random.Random(2026)
    with open(out, "w") as f:
        for i in range(n):
            br = {k: (rng.choice(LAM), rng.choice(MU), rng.choice(TS)) for k in "ABCxy"}
            try:
                O, H = predict(br, n=120)
            except (OverflowError, ZeroDivisionError):
                continue
            c = O + H["AB"]
            w = max(H["AC"], H["BC"])
            tot = c + H["AC"] + H["BC"]
            d = w - H["AB"]
            f.write(json.dumps({"br": br, "O": O, "H": H, "ovl_margin": (c - w) / tot if tot > 0 else None,
                                "q_star": O / d if d > 0 else None}) + "\n")
