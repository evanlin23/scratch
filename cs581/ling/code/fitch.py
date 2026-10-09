"""Numba Fitch parsimony + a generic 'homoplasy-penalty' objective and SPR hill-climbing.

Each character c has leaf state sets (uint64 bitmasks; polymorphic leaves = several bits),
a weight w_c and a lower bound lb_c = (#states - 1). With extra_c = fitch_c - lb_c the
objective minimised is

    sum_c w_c * g(extra_c)  +  eps * sum_c w_c * fitch_c

with g given as a lookup table. g(e) = e is (weighted) maximum parsimony; g(e) = [e > 0] is
(weighted) maximum compatibility; g(e) = min(e, cap) is "capped parsimony", which tolerates a
bounded amount of homoplasy per character and interpolates between the two (cap=1 -> MC).
eps > 0 breaks the many ties of compatibility scores by parsimony.
"""
import numpy as np
from numba import njit

from trees import ArrTree, random_tree


@njit(cache=True)
def _postorder_internal(left, right, root, n, out):
    # iterative postorder of internal nodes; returns count
    stack = np.empty(2 * len(left), np.int64)
    flag = np.zeros(2 * len(left), np.int8)
    sp = 0
    stack[0] = root; flag[0] = 0; sp = 1
    k = 0
    while sp > 0:
        sp -= 1
        v = stack[sp]; f = flag[sp]
        if v < n:
            continue
        if f == 1:
            out[k] = v; k += 1
        else:
            stack[sp] = v; flag[sp] = 1; sp += 1
            stack[sp] = right[v]; flag[sp] = 0; sp += 1
            stack[sp] = left[v]; flag[sp] = 0; sp += 1
    return k


@njit(cache=True)
def fitch_costs(left, right, root, n, leafsets, sets, order, costs):
    C = leafsets.shape[1]
    for v in range(n):
        for c in range(C):
            sets[v, c] = leafsets[v, c]
    k = _postorder_internal(left, right, root, n, order)
    for c in range(C):
        costs[c] = 0
    for i in range(k):
        v = order[i]
        a = left[v]; b = right[v]
        for c in range(C):
            x = sets[a, c] & sets[b, c]
            if x == 0:
                sets[v, c] = sets[a, c] | sets[b, c]
                costs[c] += 1
            else:
                sets[v, c] = x


@njit(cache=True)
def _score(left, right, root, n, leafsets, w, lb, gtab, eps, sets, order, costs):
    fitch_costs(left, right, root, n, leafsets, sets, order, costs)
    G = len(gtab)
    s = 0.0
    p = 0.0
    for c in range(len(costs)):
        e = costs[c] - lb[c]
        if e < 0:
            e = 0
        if e >= G:
            e = G - 1
        s += w[c] * gtab[e]
        p += w[c] * costs[c]
    return s + eps * p


@njit(cache=True)
def _in_subtree(x, p, parent):
    while x != -1:
        if x == p:
            return True
        x = parent[x]
    return False


@njit(cache=True)
def _apply_spr(L, R, P, root, p, x):
    q = P[p]
    s = R[q] if L[q] == p else L[q]
    g = P[q]
    if g == -1:
        root = s; P[s] = -1
    else:
        if L[g] == q:
            L[g] = s
        else:
            R[g] = s
        P[s] = g
    y = P[x]
    L[q] = x; R[q] = p; P[x] = q; P[p] = q; P[q] = y
    if y == -1:
        root = q
    else:
        if L[y] == x:
            L[y] = q
        else:
            R[y] = q
    return root


@njit(cache=True)
def spr_climb(left, right, parent, root, n, leafsets, w, lb, gtab, eps, seed, max_rounds):
    np.random.seed(seed)
    N = len(left)
    C = leafsets.shape[1]
    sets = np.zeros((N, C), np.uint64)
    order = np.zeros(N, np.int64)
    costs = np.zeros(C, np.int64)
    best = _score(left, right, root, n, leafsets, w, lb, gtab, eps, sets, order, costs)
    L = left.copy(); R = right.copy(); P = parent.copy()
    nevals = 0
    for rnd in range(max_rounds):
        improved = False
        prune = np.random.permutation(N)
        for pi in range(N):
            p = prune[pi]
            if p == root:
                continue
            q = parent[p]
            s = right[q] if left[q] == p else left[q]
            targets = np.random.permutation(N)
            for xi in range(N):
                x = targets[xi]
                if x == q or x == s or _in_subtree(x, p, parent):
                    continue
                L[:] = left; R[:] = right; P[:] = parent
                r2 = _apply_spr(L, R, P, root, p, x)
                sc = _score(L, R, r2, n, leafsets, w, lb, gtab, eps, sets, order, costs)
                nevals += 1
                if sc < best - 1e-9:
                    best = sc
                    left[:] = L; right[:] = R; parent[:] = P; root = r2
                    improved = True
                    break
            # continue scanning other prune nodes on the (possibly) new tree
        if not improved:
            break
    return root, best, nevals


def encode_sets(states_list, n):
    """states_list: list over characters of length-n lists of state sets (iterables of ints >=0).
    Returns leafsets (n x C uint64) and lb (C,) where lb = (#distinct states - 1), a valid
    lower bound only for monomorphic data (polymorphic leaves can lower the true minimum,
    in which case extra is clipped at 0)."""
    C = len(states_list)
    ls = np.zeros((n, C), np.uint64)
    lb = np.zeros(C, np.int64)
    for c, col in enumerate(states_list):
        labels = {}
        for i, st in enumerate(col):
            m = 0
            for x in st:
                if x not in labels:
                    labels[x] = len(labels)
                m |= 1 << labels[x]
            ls[i, c] = np.uint64(m)
        assert len(labels) <= 64
        lb[c] = len(labels) - 1
    return ls, lb


class Objective:
    def __init__(self, kind, cap=None, eps=None):
        self.kind = kind
        G = 400
        if kind == 'mp':
            self.gtab = np.arange(G, dtype=float); self.eps = 0.0 if eps is None else eps
        elif kind == 'mc':
            self.gtab = np.minimum(np.arange(G), 1).astype(float); self.eps = 1e-4 if eps is None else eps
        elif kind == 'cap':
            self.gtab = np.minimum(np.arange(G), cap).astype(float); self.eps = 1e-4 if eps is None else eps
        else:
            raise ValueError(kind)


def score_tree(t, leafsets, w, lb, obj):
    N = len(t.left); C = leafsets.shape[1]
    sets = np.zeros((N, C), np.uint64); order = np.zeros(N, np.int64); costs = np.zeros(C, np.int64)
    s = _score(t.left, t.right, t.root, t.n, leafsets, w, lb, obj.gtab, obj.eps, sets, order, costs)
    return s, costs.copy()


def char_costs(t, leafsets):
    N = len(t.left); C = leafsets.shape[1]
    sets = np.zeros((N, C), np.uint64); order = np.zeros(N, np.int64); costs = np.zeros(C, np.int64)
    fitch_costs(t.left, t.right, t.root, t.n, leafsets, sets, order, costs)
    return costs


def search(leafsets, lb, w, obj, starts, rng, max_rounds=10000):
    """SPR hill-climb from each start tree (list of ArrTree, plus 'random' strings). Best wins."""
    n = leafsets.shape[0]
    best = None
    for st in starts:
        t = random_tree(n, rng) if isinstance(st, str) else st.copy()
        root, sc, ne = spr_climb(t.left, t.right, t.parent, t.root, n, leafsets, w.astype(float),
                                 lb, obj.gtab, obj.eps, int(rng.integers(1 << 30)), max_rounds)
        t.root = root
        if best is None or sc < best[1] - 1e-9:
            best = (t, sc)
    return best
