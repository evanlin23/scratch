"""Quartet margins at S25-like rates on 26-taxon trees: true tags vs correct-root overlap vs
ASTRAL-Pro3 own rooting (apro.py, on the FULL gene tree), for random quartets whose D is on
the other side of the root.  Margin = (correct - best wrong)/total, summed over families.
usage: python realown.py SETTING NFAM NQ OUT.jsonl"""
import json, math, random, sys
from collections import defaultdict
from phylo import parse_newick
from gdlsim import yule_tree
from gdl_len import sim_family
from runmeth import parse_g
import apro

setting, NF, NQ, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
for rep in range(3):
    rng = random.Random(900 + rep)
    nwk = yule_tree(26, seed=600 + rep, height=1.0)
    st = parse_newick(nwk)
    lamH, sig = {"S25": (0.93, 0), "S25x4": (4, 0), "S25het": (2, 1.0), "S25het4": (4, 1.0)}[setting]
    lam = {v: lamH * (math.exp(rng.gauss(0, sig)) if sig else 1) for v in range(len(st.parent))}
    mu = {v: lamH * (math.exp(rng.gauss(0, sig)) if sig else 1) for v in range(len(st.parent))}
    below = {}
    for v in st.postorder():
        below[v] = [st.label[v]] if not st.children[v] else sum((below[c] for c in st.children[v]), [])
    r0, r1 = st.children[st.root]
    side, other = (below[r0], below[r1]) if len(below[r0]) >= 3 else (below[r1], below[r0])
    quartets = []
    tr = {}
    # species-tree topology of A,B,C to name the true cherry
    def lca_sp(xs):
        best = None
        for v in st.postorder():
            if set(xs) <= set(below[v]) and (best is None or len(below[v]) < len(below[best])):
                best = v
        return best
    for _ in range(NQ):
        A, B, C = rng.sample(side, 3); D = rng.choice(other)
        for x, y, z in [(A, B, C), (A, C, B), (B, C, A)]:
            if lca_sp([x, y]) != lca_sp([x, y, z]):
                quartets.append((x, y, z, D)); break
    fams = []
    while len(fams) < NF:
        r = sim_family(st, lam, mu, rng, tags=True)
        if r is None: continue
        fams.append(parse_g(r[0]))
    acc = {m: [defaultdict(float) for _ in quartets] for m in ("true", "ovl", "own")}
    for g in fams:
        adj, sp, te = apro.to_unrooted(g)
        present = set(s for s in sp if s)
        qs = [i for i, q in enumerate(quartets) if set(q) <= present]
        if not qs: continue
        kinds = []
        def rec(x):
            kinds.append(x.kind)
            for c in x.ch: rec(c)
        rec(g)
        par, ch, order = apro.rooted(adj, sp, te)
        tdup = {v: kinds[v] == "D" for v in range(len(kinds)) if sp[v] is None}
        odup, _, _ = apro.tag_and_score(sp, ch, order)
        edges = [(u, v) for u in range(len(adj)) for v in adj[u] if u < v]
        best, opt = None, []
        for e in edges:
            p2, c2, o2 = apro.rooted(adj, sp, e)
            d2, s, _ = apro.tag_and_score(sp, c2, o2)
            if best is None or s < best: best, opt = s, [(p2, c2, o2, d2)]
            elif s == best: opt.append((p2, c2, o2, d2))
        for i in qs:
            q = quartets[i]
            for k, v in apro.quartet_scores(sp, par, ch, order, tdup, q).items(): acc["true"][i][k] += v
            for k, v in apro.quartet_scores(sp, par, ch, order, odup, q).items(): acc["ovl"][i][k] += v
            for p2, c2, o2, d2 in opt:
                for k, v in apro.quartet_scores(sp, p2, c2, o2, d2, q).items(): acc["own"][i][k] += v / len(opt)
    for m in acc:
        ms = []
        for i, q in enumerate(quartets):
            d = acc[m][i]; tot = sum(d.values())
            if tot == 0: continue
            ms.append((d.get("AB|CD", 0) - max(d.get("AC|BD", 0), d.get("AD|BC", 0))) / tot)
        ms.sort()
        rec_ = {"setting": setting, "rep": rep, "mode": m, "nq": len(ms), "min": ms[0] if ms else None,
                "median": ms[len(ms) // 2] if ms else None, "nneg": sum(x < 0 for x in ms)}
        open(out, "a").write(json.dumps(rec_) + "\n"); print(json.dumps(rec_), flush=True)
