"""Shared helpers: alignment IO, distances, split extraction and FN/FP scoring."""
import numpy as np
import dendropy

ENC = {c: i for i, c in enumerate("ACGT")}


def read_fasta(path):
    names, seqs, cur = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
                names.append(line[1:].split()[0])
            else:
                cur.append(line)
    if cur:
        seqs.append("".join(cur))
    lut = np.full(256, 4, np.uint8)
    for c, i in ENC.items():
        lut[ord(c)] = i
        lut[ord(c.lower())] = i
    X = np.stack([lut[np.frombuffer(s.encode(), np.uint8)] for s in seqs])
    return names, X


def _pair_counts(X, chunk=20000):
    """same[i,j] = #identical sites, ts[i,j] = #transitions (A<->G, C<->T)."""
    n, k = X.shape
    same = np.zeros((n, n))
    ts = np.zeros((n, n))
    for s in range(0, k, chunk):
        Y = X[:, s:s + chunk]
        oh = [(Y == b).astype(np.float32) for b in range(4)]
        for b in range(4):
            same += oh[b] @ oh[b].T
        ts += oh[0] @ oh[2].T + oh[2] @ oh[0].T + oh[1] @ oh[3].T + oh[3] @ oh[1].T
    return same, ts


def jc_distance(X):
    n, k = X.shape
    same, _ = _pair_counts(X)
    p = 1.0 - same / k
    with np.errstate(divide="ignore", invalid="ignore"):
        arg = 1 - 4.0 * p / 3.0
        D = np.where(arg > 0, -0.75 * np.log(np.where(arg > 0, arg, 1)), np.inf)
    np.fill_diagonal(D, 0.0)
    return D


def k2p_distance(X):
    n, k = X.shape
    same, ts = _pair_counts(X)
    P = ts / k
    Q = 1.0 - same / k - P
    a1 = 1 - 2 * P - Q
    a2 = 1 - 2 * Q
    ok = (a1 > 0) & (a2 > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        D = np.where(ok, -0.5 * np.log(np.where(ok, a1, 1)) - 0.25 * np.log(np.where(ok, a2, 1)), np.inf)
    np.fill_diagonal(D, 0.0)
    return D


def cap_distances(D, factor=2.0):
    """Replace undefined (saturated) distances by factor * max finite distance (same matrix is
    given to every distance method)."""
    D = D.copy()
    fin = np.isfinite(D)
    mx = D[fin].max() if fin.any() else 1.0
    D[~fin] = factor * mx
    return D


def tree_splits(tree, names):
    """non-trivial splits of a dendropy tree as global bitmasks (index = position in names),
    normalized to the side NOT containing leaf 0."""
    idx = {nm: i for i, nm in enumerate(names)}
    n = len(names)
    full = (1 << n) - 1
    res = set()
    def rec(node):
        if node.is_leaf():
            return 1 << idx[node.taxon.label]
        m = 0
        for ch in node.child_nodes():
            m |= rec(ch)
        node._mask = m
        return m
    import sys
    sys.setrecursionlimit(100000)
    rec(tree.seed_node)
    for node in tree.preorder_node_iter():
        if node.is_leaf() or node is tree.seed_node:
            continue
        m = node._mask
        if m & 1:
            m = full & ~m
        k = bin(m).count("1")
        if 2 <= k <= n - 2:
            res.add(m)
    return res


true_splits = tree_splits


def splits_from_newick(nwk, names):
    tns = dendropy.TaxonNamespace(names)
    t = dendropy.Tree.get(data=nwk, schema="newick", taxon_namespace=tns, rooting="force-unrooted",
                          preserve_underscores=True)
    return tree_splits(t, names)


def score_tree(est, true, n):
    """FN rate = missed true splits / (n-3); FP rate = false splits / #estimated splits."""
    tp = len(est & true)
    fn = (len(true) - tp) / max(1, len(true))
    fp = (len(est) - tp) / max(1, len(est))
    return fn, fp, len(est) - tp


def restrict(splits, sub, n):
    """restrict global splits to leaf subset mask `sub`, normalized away from min leaf of sub."""
    low = (sub & -sub)
    out = set()
    k_sub = bin(sub).count("1")
    for s in splits:
        a = s & sub
        k = bin(a).count("1")
        if 2 <= k <= k_sub - 2:
            if a & low:
                a = sub & ~a
            out.add(a)
    return out


def score_forest(comps, comp_splits, true, n):
    """Forest scored as a partial tree.
    correct = component splits present in true tree restricted to the component;
    FN = 1 - correct/(n-3)  (resolution relative to a full binary tree)
    FP = false/(#forest splits)."""
    correct = false = 0
    for comp, S in zip(comps, comp_splits):
        sub = 0
        for i in comp:
            sub |= 1 << int(i)
        Tr = restrict(true, sub, n)
        low = sub & -sub
        for s in S:
            a = s if not (s & low) else sub & ~s
            if a in Tr:
                correct += 1
            else:
                false += 1
    tot = correct + false
    return 1 - correct / (n - 3), (false / tot if tot else 0.0), false, tot
