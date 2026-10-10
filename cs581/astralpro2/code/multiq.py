"""Limiting quartet frequencies for ASTRAL-ONE (one uniform copy per species, one quartet per
family) and ASTRAL-multi (all copy 4-tuples counted) on a 4-taxon pool; block SEs.
usage: python multiq.py POOL [NBLOCK]"""
import itertools, sys, statistics as stt
from collections import defaultdict
from runmeth import parse_g
import apro

pool = sys.argv[1]; NB = int(sys.argv[2]) if len(sys.argv) > 2 else 20
lines = [l.strip() for l in open(pool) if ";" in l]
bs = len(lines) // NB
blocks = []
for b in range(NB):
    one, multi = defaultdict(float), defaultdict(float)
    for nw in lines[b * bs:(b + 1) * bs]:
        g = parse_g(nw)
        adj, sp, te = apro.to_unrooted(g)
        par, ch, order = apro.rooted(adj, sp, te)
        depth = {-1: 0}
        for v in order[1:]: depth[v] = depth[par[v]] + 1
        def lca(u, v):
            while u != v:
                if depth[u] < depth[v]: u, v = v, u
                u = par[u]
            return u
        L = defaultdict(list)
        for v in order:
            if v != -1 and not ch[v]: L[sp[v]].append(v)
        cnt = defaultdict(int); n = 0
        for q in itertools.product(L["A"], L["B"], L["C"], L["D"]):
            # unrooted topology via four-point on path lengths (edge counts)
            dist = lambda u, v: depth[u] + depth[v] - 2 * depth[lca(u, v)]
            s = {"AB|CD": dist(q[0], q[1]) + dist(q[2], q[3]), "AC|BD": dist(q[0], q[2]) + dist(q[1], q[3]),
                 "AD|BC": dist(q[0], q[3]) + dist(q[1], q[2])}
            k = min(s, key=s.get)
            cnt[k] += 1; n += 1
        for k, v in cnt.items():
            multi[k] += v; one[k] += v / n
    blocks.append((dict(one), dict(multi)))
for i, name in enumerate(["ONE", "multi"]):
    tops = ["AB|CD", "AC|BD", "AD|BC"]
    per = [[bl[i].get(t, 0) / bs for t in tops] for bl in blocks]
    means = [stt.mean(p[j] for p in per) for j in range(3)]
    d1 = [p[0] - max(p[1], p[2]) for p in per]
    print(name, "per-family", ["%s=%.5f" % (t, m) for t, m in zip(tops, means)],
          "correct-bestwrong %.5f +- %.5f" % (stt.mean(d1), stt.stdev(d1) / NB ** 0.5))
