"""Family-R phase map: exact correct-root overlap margin vs simulated ASTRAL-Pro3 own-rooting
margin (apro.py, tie-averaged). Margins are (correct - best wrong)/total per family.
Restartable jsonl."""
import json, math, os, random, sys
from collections import defaultdict
from famR_sim import family
from apro import family_scores
from famR import margin as exact_margin

out = sys.argv[1]; N = int(sys.argv[2])
done = set()
if os.path.exists(out):
    done = {(r["a"], r["c"], r["lamT"]) for r in map(json.loads, open(out))}
for lamT in [1.0, 2.0, 3.0]:
    for a in [0.02, 0.05, 0.1, 0.2]:
        for c in [0.2, 0.4, 0.6]:
            if (a, c, lamT) in done or c <= a: continue
            m, r = exact_margin(a, a, c, lamT)
            tot = r["O"] + r["AB"] + r["AC"] + r["BC"]
            rng = random.Random(int(1000 * a + 10 * c + lamT))
            s = {"true": defaultdict(float), "own": defaultdict(float)}
            blocks = {"own": [], "true": []}
            for blk in range(10):
                bs = {"true": defaultdict(float), "own": defaultdict(float)}
                for i in range(N // 10):
                    g = family(lamT, 1.0, a, a, c, rng, [0])
                    if g is None or g.kind == "L": continue
                    for md in ("true", "own"):
                        for k, v in family_scores(g, md).items():
                            if "|" in k: bs[md][k] += v
                for md in bs:
                    d = bs[md]["AB|CD"] - max(bs[md]["AC|BD"], bs[md]["AD|BC"])
                    blocks[md].append(d / (N // 10))
                    for k, v in bs[md].items(): s[md][k] += v
            res = {"a": a, "c": c, "lamT": lamT, "exact_rel_margin": m / tot}
            for md in s:
                T = sum(s[md].values())
                d = s[md]["AB|CD"] - max(s[md]["AC|BD"], s[md]["AD|BC"])
                bl = blocks[md]; mu = sum(bl) / len(bl)
                se = (sum((x - mu) ** 2 for x in bl) / (len(bl) - 1) / len(bl)) ** 0.5
                res[md + "_rel_margin"] = d / T if T else None
                res[md + "_z"] = mu / se if se > 0 else None
            open(out, "a").write(json.dumps(res) + "\n")
            print(json.dumps(res), flush=True)
