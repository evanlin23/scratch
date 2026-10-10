"""Independent re-implementation of ASTRAL-Pro's per-family quartet scores for one set of
four species, with (i) a given root + species-overlap tags, or (ii) ASTRAL-Pro3's own
rooting: minimise sum over D-nodes of 1 + [L != L u R] + [R != L u R] (scoreSubtree in
ASTER src/astral-pro.cpp), ties broken uniformly at random (we return the exact average
over tied optimal roots, i.e. the expectation of the tie-break).
SQ = four leaves of four species whose every triple LCA is S; anchor LCA as in Zhang et
al. 2020 Def. 3; per-family score of a topology = #distinct anchor LCAs (Def. 4/5)."""
import itertools
from collections import defaultdict

from gdlsim import GNode


def to_unrooted(g, sp_of=lambda l: l.rsplit("_", 1)[0]):
    """GNode -> (adj list, leaf species list or None for internal). Root degree-2 node suppressed."""
    adj, sp = [], []

    def new(s):
        adj.append([]); sp.append(s); return len(adj) - 1

    def rec(x):
        if x.kind == "L":
            return new(sp_of(x.label))
        v = new(None)
        for c in x.ch:
            u = rec(c)
            adj[v].append(u); adj[u].append(v)
        return v
    r = rec(g)
    true_edge = None
    if sp[r] is None and len(adj[r]) == 2:          # suppress root
        a, b = adj[r]
        adj[a][adj[a].index(r)] = b; adj[b][adj[b].index(r)] = a
        adj[r] = []
        true_edge = (a, b)
    return adj, sp, true_edge


def rooted(adj, sp, edge):
    """Root on edge (a,b). Returns parent dict, children dict with a virtual root id -1."""
    a, b = edge
    par, ch = {a: -1, b: -1}, {-1: [a, b]}
    order, st = [-1], [a, b]
    other = {a: b, b: a}
    while st:
        v = st.pop(); order.append(v)
        ch[v] = [u for u in adj[v] if u != par[v] and u != other.get(v)]
        for u in ch[v]:
            par[u] = v; st.append(u)
    return par, ch, order


def tag_and_score(sp, ch, order):
    """species-overlap tags + ASTRAL-Pro DL score."""
    sets, dup, score = {}, {}, 0
    for v in reversed(order):
        if not ch[v]:
            sets[v] = frozenset([sp[v]]); continue
        L, R = sets[ch[v][0]], sets[ch[v][1]]
        U = L | R
        sets[v] = U
        if L & R:
            dup[v] = True
            score += 1 + (L != U) + (R != U)
        else:
            dup[v] = False
    return dup, score, sets


def quartet_scores(sp, par, ch, order, dup, species=("A", "B", "C", "D")):
    """{topology: #classes}; topology key 'AB|CD' etc. Assumes binary tree."""
    depth = {-1: 0}
    for v in order[1:]:
        depth[v] = depth[par[v]] + 1
    def lca(u, v):
        while u != v:
            if depth[u] < depth[v]: u, v = v, u
            u = par[u]
        return u
    leaves = defaultdict(list)
    for v in order:
        if v != -1 and not ch[v]:
            leaves[sp[v]].append(v)
    A, B, C, D = species
    anchors = defaultdict(set)
    for q in itertools.product(leaves[A], leaves[B], leaves[C], leaves[D]):
        ok = True
        for tri in itertools.combinations(q, 3):
            w = lca(lca(tri[0], tri[1]), tri[2])
            if dup.get(w, False):
                ok = False; break
        if not ok:
            continue
        # rooted restriction: deepest pairwise LCA is a cherry (i,j); topology ij|km
        l = {(i, j): lca(q[i], q[j]) for i, j in itertools.combinations(range(4), 2)}
        i, j = max(l, key=lambda p: depth[l[p]])
        k, m = [x for x in range(4) if x not in (i, j)]
        r = lca(l[(i, j)], l[(k, m)])
        if l[(k, m)] != r and lca(q[i], q[k]) == r:
            anc = r                                   # balanced ((ij),(km)): anchors u, LCA(km)
        else:
            wk, wm = lca(l[(i, j)], q[k]), lca(l[(i, j)], q[m])
            anc = wk if depth[wk] > depth[wm] else wm  # caterpillar: anchor LCA = 2nd internal node
        top = frozenset(["ABCD"[i], "ABCD"[j]])   # positional names: (A,B,C,D) = species order
        key = "".join(sorted(top)) + "|" + "".join(sorted(set("ABCD") - top))
        if "A" not in key.split("|")[0]:
            key = key.split("|")[1] + "|" + key.split("|")[0]
        anchors[key].add(anc)
    return {k: len(v) for k, v in anchors.items()}


def family_scores(g, mode="true", species=("A", "B", "C", "D")):
    """mode 'true': true root + overlap tags. 'own': ASTRAL-Pro3 rooting, average over ties.
    Returns dict topology -> expected #classes, plus 'nroots' (ties) for own."""
    adj, sp, te = to_unrooted(g)
    if mode == "true":
        par, ch, order = rooted(adj, sp, te)
        dup, _, _ = tag_and_score(sp, ch, order)
        return quartet_scores(sp, par, ch, order, dup, species)
    edges = [(u, v) for u in range(len(adj)) for v in adj[u] if u < v]
    best, res = None, []
    for e in edges:
        par, ch, order = rooted(adj, sp, e)
        dup, s, _ = tag_and_score(sp, ch, order)
        if best is None or s < best:
            best, res = s, [(par, ch, order, dup, e)]
        elif s == best:
            res.append((par, ch, order, dup, e))
    tot = defaultdict(float)
    for par, ch, order, dup, e in res:
        for k, v in quartet_scores(sp, par, ch, order, dup, species).items():
            tot[k] += v / len(res)
    out = dict(tot)
    out["nroots"] = len(res)
    out["true_in_opt"] = any(set(e) == set(te) for *_, e in res) if te else None
    return out


# ---------------------------------------------------------------- reconciliation (LCA)
SPECIES_TREES = {   # rooted 4-taxon species trees as nested tuples
    "AB|CD": ((("A", "B"), "C"), "D"),
}


def _st_index(st):
    """nested tuple species tree -> (clusters list, map species-set -> smallest cluster id)."""
    cl = []
    def rec(x):
        if isinstance(x, str):
            s = frozenset([x])
        else:
            s = frozenset().union(*[rec(c) for c in x])
        cl.append(s)
        return s
    rec(st)
    return cl


def recon_root_tag(adj, sp, st):
    """DupTree-style: root minimising #duplications under LCA mapping to rooted species tree st
    (ties: fewest losses, then uniform). Returns list of (par, ch, order, dup) optimal roots."""
    cl = _st_index(st)
    def M(s):  # smallest cluster containing s
        return min((c for c in cl if s <= c), key=len)
    def depthc(c):
        return sum(1 for d in cl if c < d)
    edges = [(u, v) for u in range(len(adj)) for v in adj[u] if u < v]
    best, res = None, []
    for e in edges:
        par, ch, order = rooted(adj, sp, e)
        sets, mp, dup, nd, nl = {}, {}, {}, 0, 0
        for v in reversed(order):
            if not ch[v]:
                sets[v] = frozenset([sp[v]]); mp[v] = M(sets[v]); continue
            sets[v] = sets[ch[v][0]] | sets[ch[v][1]]
            mp[v] = M(sets[v])
            d = mp[v] == mp[ch[v][0]] or mp[v] == mp[ch[v][1]]
            dup[v] = d
            nd += d
            for c in ch[v]:   # losses on edge v->c (standard LCA loss count)
                nl += depthc(mp[c]) - depthc(mp[v]) - (0 if d else 1)
        key = (nd, nl)
        if best is None or key < best:
            best, res = key, [(par, ch, order, dup)]
        elif key == best:
            res.append((par, ch, order, dup))
    return res, best


def truetag_scores(g, species=("A", "B", "C", "D")):
    adj, sp, te = to_unrooted(g)
    # tags from GNode kinds: rebuild by walking g in the same order as to_unrooted
    kinds = []
    def rec(x):
        kinds.append(x.kind)
        for c in x.ch:
            rec(c)
    rec(g)
    par, ch, order = rooted(adj, sp, te)
    dup = {v: kinds[v] == "D" for v in range(len(kinds)) if sp[v] is None}
    return quartet_scores(sp, par, ch, order, dup, species)


def recon_scores(g, st, species=("A", "B", "C", "D")):
    adj, sp, te = to_unrooted(g)
    res, _ = recon_root_tag(adj, sp, st)
    tot = defaultdict(float)
    for par, ch, order, dup in res:
        for k, v in quartet_scores(sp, par, ch, order, dup, species).items():
            tot[k] += v / len(res)
    return dict(tot)


def gtp_cost(g, st, loss_w=1.0):
    """min over gene-tree roots of dups + loss_w * losses w.r.t. rooted species tree st."""
    adj, sp, te = to_unrooted(g)
    cl = _st_index(st)
    M = lambda s: min((c for c in cl if s <= c), key=len)
    depthc = lambda c: sum(1 for d in cl if c < d)
    best = None
    for e in [(u, v) for u in range(len(adj)) for v in adj[u] if u < v]:
        par, ch, order = rooted(adj, sp, e)
        sets, mp, cost = {}, {}, 0.0
        for v in reversed(order):
            if not ch[v]:
                sets[v] = frozenset([sp[v]]); mp[v] = M(sets[v]); continue
            sets[v] = sets[ch[v][0]] | sets[ch[v][1]]; mp[v] = M(sets[v])
            d = mp[v] == mp[ch[v][0]] or mp[v] == mp[ch[v][1]]
            cost += d
            for c in ch[v]:
                cost += loss_w * (depthc(mp[c]) - depthc(mp[v]) - (0 if d else 1))
        if best is None or cost < best:
            best = cost
    return best
