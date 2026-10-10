"""Prevalence of ASTRAL-Pro inconsistency (correct root + overlap tags; exact formula of
predict_margin.py) for induced caterpillar quartets of S25-like species trees
(26-taxon ultrametric Yule, height 1; D taken from the other side of the root, so the
family starts at the quartet root as in SimPhy).  Rates: lam*H on a grid, mu = r*lam,
optional per-branch log-normal multipliers (sigma) on lam and mu independently.
ASTRAL-Pro / DISCO S25 default: lam*H = 4.9e-10 * 1.9e9 ~ 0.93, r = 1, sigma = 0.
usage: python realscan.py OUT.jsonl   (restartable: skips finished grid cells)"""
import json, math, os, random, sys
from phylo import parse_newick
from gdlsim import yule_tree
from predict_margin import predict

out = sys.argv[1]
done = set()
if os.path.exists(out):
    for l in open(out):
        r = json.loads(l); done.add((r["lamH"], r["r"], r["sigma"]))

trees = []
for s in range(40):
    t = parse_newick(yule_tree(26, seed=1000 + s, height=1.0))
    trees.append(t)


def height(t):
    h = {}
    for v in t.postorder():
        h[v] = 0.0 if not t.children[v] else max(h[c] + t.length[c] for c in t.children[v])
    return h


def caterpillars(t, rng, k):
    h = height(t)
    leaves_below = {}
    for v in t.postorder():
        leaves_below[v] = [v] if not t.children[v] else sum((leaves_below[c] for c in t.children[v]), [])
    r = t.root
    sides = t.children[r]
    res = []
    for _ in range(k * 5):
        if len(res) >= k: break
        i = rng.randrange(2)
        ing, outg = leaves_below[sides[i]], leaves_below[sides[1 - i]]
        if len(ing) < 3: continue
        A, B, C = rng.sample(ing, 3)
        D = rng.choice(outg)
        # LCAs
        def anc(v):
            s = [v]
            while t.parent[v] is not None and t.parent[v] != -1:
                v = t.parent[v]; s.append(v)
            return s
        def lca(u, v):
            a = set(anc(u))
            for w in anc(v):
                if w in a: return w
        lab, lac, lbc = lca(A, B), lca(A, C), lca(B, C)
        # make (A,B) the cherry
        if lac != lab and h[lac] < h[lab]: B, C = C, B; lab = lac
        elif lbc != lab and h[lbc] < h[lab]: A, C = C, A; lab = lbc
        y = lca(lab, C)
        tx, ty = h[lab], h[y]
        hr = h[r]
        res.append({"A": tx, "B": tx, "C": ty, "x": ty - tx, "y": hr - ty})
    return res


LAMH = [0.25, 0.5, 1, 2, 4, 8, 16]
RS = [0.0, 0.5, 1.0, 2.0]
SIG = [0.0, 0.5, 1.0]
with open(out, "a") as f:
    for lamH in LAMH:
        for rr in RS:
            for sig in SIG:
                if (lamH, rr, sig) in done: continue
                rng = random.Random(hash((lamH, rr, sig)) & 0xffff)
                n = nfail = ninf = 0
                minm, worst = 9, None
                for t in trees:
                    for L in caterpillars(t, rng, 10):
                        br = {}
                        for k in "ABCxy":
                            ml = math.exp(rng.gauss(0, sig)) if sig else 1
                            mm = math.exp(rng.gauss(0, sig)) if sig else 1
                            br[k] = (lamH * ml, rr * lamH * mm, L[k])
                        try:
                            O, H = predict(br, n=200)
                        except (OverflowError, ZeroDivisionError):
                            continue
                        tot = O + sum(H.values())
                        if tot < 1e-12: continue
                        n += 1
                        m = (O + H["AB"] - max(H["AC"], H["BC"])) / tot
                        if m < minm: minm, worst = m, br
                        nfail += m < 0
                rec = {"lamH": lamH, "r": rr, "sigma": sig, "n": n, "nfail": nfail, "min_margin": minm, "worst": worst}
                f.write(json.dumps(rec) + "\n"); f.flush()
                print(lamH, rr, sig, n, nfail, round(minm, 4), flush=True)
