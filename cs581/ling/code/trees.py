"""Small tree utilities: rooted-binary array trees, Newick I/O, bipartitions, FN/RF, NJ.

A tree on n leaves is stored as arrays left/right/parent of length 2n-1 plus a root index.
Leaves are 0..n-1; internal nodes n..2n-2. The rooting is arbitrary (all methods here are
unrooted); bipartitions are compared unrooted.
"""
import re
import numpy as np


class ArrTree:
    def __init__(self, n, left, right, parent, root, blen=None):
        self.n, self.left, self.right, self.parent, self.root = n, left, right, parent, root
        self.blen = blen  # optional per-node branch length (edge above node)

    def copy(self):
        return ArrTree(self.n, self.left.copy(), self.right.copy(), self.parent.copy(), self.root,
                       None if self.blen is None else self.blen.copy())


def from_nested(nested, names):
    """nested: tuple tree of leaf names (binary). Returns ArrTree with leaves indexed by names."""
    n = len(names)
    idx = {nm: i for i, nm in enumerate(names)}
    left = -np.ones(2 * n - 1, np.int64); right = left.copy(); parent = left.copy()
    nxt = [n]

    def rec(t):
        if isinstance(t, str):
            return idx[t]
        assert len(t) == 2, t
        a, b = rec(t[0]), rec(t[1])
        v = nxt[0]; nxt[0] += 1
        left[v], right[v] = a, b
        parent[a] = parent[b] = v
        return v
    root = rec(nested)
    assert nxt[0] == 2 * n - 1, "tree must be binary and contain all leaves"
    return ArrTree(n, left, right, parent, root)


def parse_newick(s):
    """Parse Newick into nested tuples (names only, lengths ignored). Multifurcations are
    resolved as left-combs (only used for binary trees here); a trifurcating root is fine."""
    s = s.strip().rstrip(';')
    s = re.sub(r':[-0-9.eE+]+', '', s)
    s = re.sub(r'\)[^,():;]+', ')', s)  # drop internal labels / support values
    pos = [0]

    def rec():
        if s[pos[0]] == '(':
            pos[0] += 1
            kids = [rec()]
            while s[pos[0]] == ',':
                pos[0] += 1
                kids.append(rec())
            assert s[pos[0]] == ')'
            pos[0] += 1
            t = kids[0]
            for k in kids[1:]:
                t = (t, k)
            return t
        m = re.match(r'[^,():;]+', s[pos[0]:])
        pos[0] += m.end()
        return m.group(0).strip()
    return rec()


def to_newick(t, names):
    def rec(v):
        if v < t.n:
            return names[v]
        return '(' + rec(t.left[v]) + ',' + rec(t.right[v]) + ')'
    return rec(t.root) + ';'


def bipartitions(t):
    """Set of nontrivial splits as frozen int bitmasks normalized to exclude leaf 0."""
    n = t.n
    full = (1 << n) - 1
    mask = {}
    order = postorder(t)
    for v in order:
        if v < n:
            mask[v] = 1 << v
        else:
            mask[v] = mask[t.left[v]] | mask[t.right[v]]
    out = set()
    for v in order:
        if v == t.root:
            continue
        m = mask[v]
        if m & 1:
            m = full ^ m
        c = bin(m).count('1')
        if 2 <= c <= n - 2:
            out.add(m)
    return out


def postorder(t):
    out, stack = [], [(t.root, False)]
    while stack:
        v, done = stack.pop()
        if v < t.n or done:
            out.append(v)
        else:
            stack.append((v, True)); stack.append((t.right[v], False)); stack.append((t.left[v], False))
    return out


def fn_fp(true_t, est_t):
    a, b = bipartitions(true_t), bipartitions(est_t)
    fn = len(a - b) / max(1, len(a))
    fp = len(b - a) / max(1, len(b)) if b else 0.0
    return fn, fp, len(a - b) + len(b - a)


def splits_from_clades(clades, names):
    """Bitmask splits for named clades (list of lists of names)."""
    idx = {nm: i for i, nm in enumerate(names)}
    full = (1 << len(names)) - 1
    out = []
    for c in clades:
        m = 0
        for nm in c:
            m |= 1 << idx[nm]
        if m & 1:
            m = full ^ m
        out.append(m)
    return out


def nj(D):
    """Neighbor joining on distance matrix D (n x n). Returns ArrTree (rooted at last join)."""
    n = D.shape[0]
    D = D.astype(float).copy()
    left = -np.ones(2 * n - 1, np.int64); right = left.copy(); parent = left.copy()
    nxt = n
    M = D.copy()
    ids = list(range(n))  # tree node for each row of M
    while len(ids) > 3:
        m = len(ids)
        r = M.sum(1)
        Q = (m - 2) * M - r[:, None] - r[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = np.unravel_index(np.argmin(Q), Q.shape)
        if i > j:
            i, j = j, i
        v = nxt; nxt += 1
        left[v], right[v] = ids[i], ids[j]
        parent[ids[i]] = parent[ids[j]] = v
        dnew = 0.5 * (M[i] + M[j] - M[i, j])
        keep = [k for k in range(m) if k not in (i, j)]
        M2 = np.zeros((m - 1, m - 1))
        M2[:m - 2, :m - 2] = M[np.ix_(keep, keep)]
        M2[m - 2, :m - 2] = dnew[keep]; M2[:m - 2, m - 2] = dnew[keep]
        M = M2
        ids = [ids[k] for k in keep] + [v]
    # join last three: ((a,b),c)
    a, b, c = ids
    v = nxt; nxt += 1
    left[v], right[v] = a, b; parent[a] = parent[b] = v
    r = nxt; nxt += 1
    left[r], right[r] = v, c; parent[v] = parent[c] = r
    return ArrTree(n, left, right, parent, r)


def random_tree(n, rng):
    """Random binary topology by random sequential addition (rooted representation)."""
    left = -np.ones(2 * n - 1, np.int64); right = left.copy(); parent = left.copy()
    left[n], right[n] = 0, 1; parent[0] = parent[1] = n
    root = n
    nodes = [0, 1, n]
    for k in range(2, n):
        x = nodes[rng.integers(len(nodes))]
        q = n + k - 1
        y = parent[x]
        left[q], right[q] = x, k
        parent[x] = q; parent[k] = q; parent[q] = y
        if y == -1:
            root = q
        elif left[y] == x:
            left[y] = q
        else:
            right[y] = q
        nodes += [k, q]
    return ArrTree(n, left, right, parent, root)
