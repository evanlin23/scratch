"""Fast re-implementation of the Daskalakis-Mossel-Roch (2011) Forest algorithm.

Follows the logic of lanl/distphylo (forest_algorithm.py + prune_deep_helper_functions.py):
  1. clustering graph: edge (i,j) iff d(i,j) < m; each connected component is solved separately
  2. MiniContractor for every leaf pair (u,v) of a component: ball B = {w : max(d_uw,d_vw) < M},
     sort ball leaves by phi(w) = (d_uv + d_uw - d_vw)/2, cut wherever consecutive phi gap >= 2*tau
  3. Extender: if the ball is not the whole component, remove edges between the two sides of the
     mini-bipartition from the component graph and attach each connected piece to the side it touches
  4. tree popping: the set of unique bipartitions must be pairwise compatible, else no forest

Splits are python-int bitmasks over global leaf indices. Differences from distphylo:
  * exact pairwise compatibility check instead of 10 random tree-popping orders
  * all tau values for one (m, M) are evaluated from one sort per leaf pair
"""
import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.sparse import csr_matrix


def popcount(x):
    return bin(x).count("1")


def mask_of(idx):
    m = 0
    for i in idx:
        m |= 1 << int(i)
    return m


def components(D, m):
    A = (D < m)
    np.fill_diagonal(A, False)
    nc, lab = connected_components(csr_matrix(A), directed=False)
    return [np.flatnonzero(lab == c) for c in range(nc)]


def _bits(boolarr):
    return int.from_bytes(np.packbits(boolarr, bitorder="little").tobytes(), "little")


def _extend(adj, side_u, side_v):
    """Extender with int bitmasks (local indices): remove edges between the two sides of the
    mini-bipartition, return the union of connected pieces touching side_u (distphylo assigns a
    piece touching both sides to u, pieces touching neither are dropped)."""
    visited = side_u
    frontier = side_u
    while frontier:
        low = frontier & -frontier
        x = low.bit_length() - 1
        frontier ^= low
        nb = adj[x]
        if side_u & low:
            nb &= ~side_v
        elif side_v & low:
            nb &= ~side_u
        new = nb & ~visited
        if new:
            visited |= new
            frontier |= new
    return visited


def run_forest(D, m, M, taus):
    """Forest for one (m, M) and a list of taus (one sort per leaf pair shared by all taus).
    Returns {tau: dict(valid=bool, comps=[global idx arrays], splits=[set per comp])}."""
    comps = components(D, m)
    taus = list(taus)
    res = {t: dict(valid=True, comps=comps, splits=[set() for _ in comps]) for t in taus}
    A_all = (D < m)
    np.fill_diagonal(A_all, False)
    tmin2 = 2 * min(taus)
    for ci, comp in enumerate(comps):
        c = len(comp)
        if c < 4:
            continue
        Dc = D[np.ix_(comp, comp)]
        Ac = A_all[np.ix_(comp, comp)]
        adj = [_bits(Ac[i]) for i in range(c)]
        full = (1 << c) - 1
        cache = {}
        found = {}  # local split mask -> largest gap that produced it
        for u in range(c):
            du = Dc[u]
            for v in range(u + 1, c):
                dv = Dc[v]
                inball = np.maximum(du, dv) < M
                inball[u] = True
                inball[v] = True
                S = np.flatnonzero(inball)
                S = S[S != u]
                phi = 0.5 * (du[v] + du[S] - dv[S])
                order = np.argsort(phi, kind="stable")
                gaps = np.diff(np.concatenate(([0.0], phi[order])))
                cuts = np.flatnonzero(gaps >= tmin2)
                if len(cuts) == 0:
                    continue
                whole = bool(inball.all())
                ballmask = full if whole else _bits(inball)
                rank = np.full(c, -1)
                rank[S[order]] = np.arange(len(S))
                for i in cuts:
                    side_v = _bits(rank >= i)
                    side_u = ballmask & ~side_v
                    if whole:
                        ext = side_u
                    else:
                        key = (ballmask, side_u)
                        ext = cache.get(key)
                        if ext is None:
                            ext = _extend(adj, side_u, side_v)
                            cache[key] = ext
                    k = popcount(ext)
                    if k < 2 or k > c - 2:
                        continue
                    if ext & 1:
                        ext = full & ~ext
                    g = gaps[i]
                    if found.get(ext, -1.0) < g:
                        found[ext] = g
        gm = [1 << int(x) for x in comp]
        glob = {}
        for s_loc, g in found.items():
            mm = 0
            x = s_loc
            while x:
                low = x & -x
                mm |= gm[low.bit_length() - 1]
                x ^= low
            glob[mm] = g
        cmask = mask_of(comp)
        for t in taus:
            sp = {s for s, g in glob.items() if g >= 2 * t}
            res[t]["splits"][ci] = sp
            if res[t]["valid"] and not compatible_set(sp, cmask):
                res[t]["valid"] = False
    return res


def compatible(a, b, full):
    ac, bc = full & ~a, full & ~b
    return (a & b) == 0 or (a & bc) == 0 or (ac & b) == 0 or (ac & bc) == 0


def compatible_set(splits, full):
    s = list(splits)
    for i in range(len(s)):
        for j in range(i + 1, len(s)):
            if not compatible(s[i], s[j], full):
                return False
    return True


def splits_to_newick(splits, leaves, names):
    """Build an (unrooted, possibly multifurcating) newick tree on `leaves` (global idx) from a
    compatible set of splits, each a bitmask over global indices restricted to `leaves`."""
    leaves = [int(x) for x in leaves]
    if len(leaves) == 1:
        return names[leaves[0]] + ";"
    full = mask_of(leaves)
    r = leaves[0]
    clades = []
    for s in splits:
        cl = s if not (s >> r) & 1 else full & ~s
        clades.append(cl)
    clades = sorted(set(clades), key=popcount)
    # children lists: each clade's parent is the smallest strictly-containing clade (or root)
    children = {full: []}
    for cl in clades:
        children[cl] = []
    allc = clades + [full]
    for cl in clades:
        for p in allc:
            if p != cl and (cl & p) == cl and popcount(p) > popcount(cl):
                children[p].append(cl)
                break
    def render(cl):
        covered = 0
        parts = []
        for ch in children[cl]:
            parts.append(render(ch))
            covered |= ch
        rest = cl & ~covered
        for i in leaves:
            if (rest >> i) & 1:
                parts.append(names[i])
        return "(" + ",".join(parts) + ")" if len(parts) > 1 else parts[0]
    return render(full) + ";"
