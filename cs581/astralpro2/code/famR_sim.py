"""Simulate family R gene trees directly (Yule on y-branch, per-copy patterns) and
compare: exact formula vs simulation (true root + overlap tags) vs ASTRAL-Pro3's own
rooting (re-implementation, averaged over tied roots).  Restartable: appends jsonl."""
import json, math, random, sys
from collections import defaultdict
from gdlsim import GNode
from apro import family_scores


def yule(T, lam, rng):
    """returns nested structure of the y-branch dup tree: 'TIP' or ('D', l, r)."""
    def rec(rem):
        dt = rng.expovariate(lam) if lam > 0 else math.inf
        if dt >= rem:
            return "TIP"
        return ("D", rec(rem - dt), rec(rem - dt))
    return rec(T)


def family(T, lam, a, b, c, rng, ctr):
    def ycopy():
        i = ctr[0]; ctr[0] += 1
        xa = GNode("L", label="A_%d" % i) if rng.random() < a else None
        xb = GNode("L", label="B_%d" % i) if rng.random() < b else None
        xc = GNode("L", label="C_%d" % i) if rng.random() < c else None
        x = GNode("S", [xa, xb]) if xa and xb else (xa or xb)
        return GNode("S", [x, xc]) if x and xc else (x or xc)
    def fill(s):
        if s == "TIP":
            return ycopy()
        l, r = fill(s[1]), fill(s[2])
        if l is None: return r
        if r is None: return l
        return GNode("D", [l, r])
    top = fill(yule(T, lam, rng))
    if top is None:
        return None
    return GNode("S", [top, GNode("L", label="D_0")])


if __name__ == "__main__":
    a, b, c, lamT, n, seed = map(float, sys.argv[1:7])
    out = sys.argv[7]
    rng = random.Random(int(seed))
    tot = {m: defaultdict(float) for m in ("true", "own")}
    nacc, ntie, nopt = 0, 0, 0
    for it in range(int(n)):
        g = family(lamT, 1.0, a, b, c, rng, [0])
        if g is None or g.kind == "L":
            continue
        for m in ("true", "own"):
            r = family_scores(g, m)
            for k, v in r.items():
                if k in ("nroots", "true_in_opt"):
                    continue
                tot[m][k] += v
            if m == "own":
                ntie += r["nroots"] > 1
                nopt += bool(r["true_in_opt"])
        nacc += 1
    rec = {"a": a, "b": b, "c": c, "lamT": lamT, "n": int(n), "seed": int(seed),
           "true": {k: v / n for k, v in tot["true"].items()}, "own": {k: v / n for k, v in tot["own"].items()},
           "frac_true_root_optimal": nopt / max(nacc, 1), "frac_ties": ntie / max(nacc, 1)}
    with open(out, "a") as f:
        f.write(json.dumps(rec) + "\n")
    print(json.dumps(rec))
