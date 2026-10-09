"""Internode-distance species-tree methods (NJst / ASTRID) and a missing-data correction.

Distances
---------
For a gene tree g on leaf set S_g and leaves i, j in S_g, the internode distance is the
number of internal nodes on the i-j path. NJst/ASTRID average it over the genes that
contain both i and j, then run NJ (NJst) or FastME (ASTRID).

Missing-data correction ("ASTRID-HT")
-------------------------------------
Under i.i.d. taxon deletion with rate p, an internal node v of the full gene tree on the
i-j path survives (is still a node of the restricted tree) iff at least one of the K_v
leaves on its off-path side survives; P = 1 - p^K_v. We observe k_v ~ Bin(K_v, 1-p) and
the node iff k_v >= 1. The unique unbiased estimator of the indicator "node exists" that
depends only on k_v is

    f(k) = 1 - (-p/(1-p))^k      (f(0) = 0),

since sum_k C(K,k) q^k p^(K-k) (1 - (-p/q)^k) = (q+p)^K - (p-p)^K = 1 for K >= 1.
So sum over observed nodes of f(k_v) is an unbiased estimator of the full-tree internode
distance, conditional on i, j present. Averaging over genes gives a consistent estimator
of the complete-data NJst distance (which is additive on the species tree), hence a
statistically consistent method under i.i.d. deletion. The variance explodes for p >= 1/2
(|p/q| >= 1), so we also provide a clipped version.
"""
import subprocess, tempfile, os, itertools
import numpy as np
import treeswift as ts


# ---------------------------------------------------------------- tree utilities
def read_trees(path, limit=None):
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            out.append(line)
            if limit and len(out) >= limit:
                break
    return out


def parse(nwk):
    t = ts.read_tree_newick(nwk)
    t.suppress_unifurcations()
    return t


def to_adj(t):
    """Return (labels of leaves, adjacency list over node ids, leaf ids)."""
    nodes = list(t.traverse_preorder())
    idx = {id(n): k for k, n in enumerate(nodes)}
    adj = [[] for _ in nodes]
    for n in nodes:
        for c in n.children:
            adj[idx[id(n)]].append(idx[id(c)])
            adj[idx[id(c)]].append(idx[id(n)])
    leaves = [idx[id(n)] for n in nodes if n.is_leaf()]
    labels = [nodes[k].label for k in leaves]
    # unroot: a degree-2 root is not a node of the unrooted tree
    return labels, adj, leaves


def gene_contrib(nwk, taxa_index, mode="plain", p=0.0, clip=None, count_root=False):
    """Return (D_add, present_mask) for one gene tree.

    D_add[a,b] = (corrected) number of internal nodes on the path between taxa a and b.
    mode: 'plain' -> each internal node counts 1
          'ht'    -> each internal node counts f(k) with k = #observed leaves off the path
    """
    labels, adj, leaves = to_adj(parse(nwk))
    m = len(labels)
    if m < 3:
        tix = np.array([taxa_index[l] for l in labels])
        return tix, np.zeros((m, m))
    N = len(adj)
    leafpos = {lf: r for r, lf in enumerate(leaves)}
    # For each internal node of degree >= 3, the leaf sets of its incident components.
    # Root the tree at leaf 0 for subtree sizes.
    root = leaves[0]
    parent = [-1] * N
    order = []
    stack = [root]
    seen = [False] * N
    seen[root] = True
    while stack:
        u = stack.pop()
        order.append(u)
        for w in adj[u]:
            if not seen[w]:
                seen[w] = True
                parent[w] = u
                stack.append(w)
    below = [None] * N  # list of leaf positions below node (rooted at leaf 0)
    for u in reversed(order):
        if u in leafpos and u != root:
            below[u] = [leafpos[u]]
        else:
            acc = []
            for w in adj[u]:
                if w != parent[u]:
                    acc.extend(below[w])
            below[u] = acc
    allpos = np.arange(m)
    D = np.zeros((m, m))
    q = 1.0 - p
    for u in range(N):
        if u in leafpos:
            continue
        comps = [below[w] for w in adj[u] if w != parent[u]]
        if parent[u] != -1:
            mask = np.ones(m, bool)
            for c in comps:
                mask[c] = False
            comps.append(list(allpos[mask]))
        if len(comps) < 3 and not count_root:
            continue  # degree-2 node (root of a rooted input) is not an internal node
        sizes = [len(c) for c in comps]
        for a, b in itertools.combinations(range(len(comps)), 2):
            k = m - sizes[a] - sizes[b]
            if mode == "plain":
                w = 1.0
            elif mode == "expected":
                # exact expectation of the plain observed count under iid deletion (full trees in)
                w = 1.0 - p ** k
            elif mode == "ht1":
                # first-order HT: exact for K=1 nodes, 1 otherwise (bounded by 1/q)
                w = 1.0 / q if k == 1 else 1.0
            elif mode == "plugin":
                # plug-in: K_hat = k/q, weight 1/(1 - p^K_hat)  (bounded, biased, low variance)
                w = 1.0 / (1.0 - p ** (k / q)) if p > 0 else 1.0
            else:
                w = 1.0 - (-p / q) ** k if p > 0 else 1.0
                if clip is not None:
                    w = min(max(w, 0.0), clip)
            A, B = comps[a], comps[b]
            D[np.ix_(A, B)] += w
            D[np.ix_(B, A)] += w
    tix = np.array([taxa_index[l] for l in labels])
    return tix, D


def distance_matrix(trees, taxa, mode="plain", p=None, clip=None, per_gene_p=False, count_root=False):
    """Average (corrected) internode distance over genes containing both taxa."""
    n = len(taxa)
    ti = {t: k for k, t in enumerate(taxa)}
    S = np.zeros((n, n))
    C = np.zeros((n, n))
    for nwk in trees:
        if mode != "plain" and per_gene_p:
            m = nwk.count(",") + 1
            pg = 1.0 - m / n
        else:
            pg = p or 0.0
        tix, D = gene_contrib(nwk, ti, mode=mode, p=pg, clip=clip, count_root=count_root)
        S[np.ix_(tix, tix)] += D
        C[np.ix_(tix, tix)] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        M = S / C
    np.fill_diagonal(M, 0.0)
    return M, C


def estimate_p(trees, n):
    sizes = [t.count(",") + 1 for t in trees]
    return 1.0 - np.mean(sizes) / n


# ---------------------------------------------------------------- tree from distances
def nj(M, taxa):
    """Plain neighbor joining (NJst). Returns newick. Missing entries not allowed."""
    D = M.copy().astype(float)
    nodes = list(taxa)
    while len(nodes) > 3:
        r = len(nodes)
        s = D.sum(1)
        Q = (r - 2) * D - s[:, None] - s[None, :]
        np.fill_diagonal(Q, np.inf)
        i, j = np.unravel_index(np.argmin(Q), Q.shape)
        if i > j:
            i, j = j, i
        new = f"({nodes[i]},{nodes[j]})"
        dnew = 0.5 * (D[i] + D[j] - D[i, j])
        keep = [k for k in range(r) if k not in (i, j)]
        D2 = np.zeros((r - 1, r - 1))
        D2[:-1, :-1] = D[np.ix_(keep, keep)]
        D2[-1, :-1] = dnew[keep]
        D2[:-1, -1] = dnew[keep]
        D = D2
        nodes = [nodes[k] for k in keep] + [new]
    return "(" + ",".join(nodes) + ");"


FASTME = os.environ.get("FASTME", "fastme")


def fastme(M, taxa, mode="astrid"):
    """ASTRID's tree step: FastME BioNJ start + NNI/SPR on balanced minimum evolution."""
    n = len(taxa)
    M = M.copy()
    if np.isnan(M).any():  # ASTRID fills missing entries; here: replace by row means (rare)
        fill = np.nanmax(M)
        M[np.isnan(M)] = fill
    off = ~np.eye(n, dtype=bool)
    if M[off].min() < 0:  # BME/NJ topology is invariant to a constant shift of off-diagonal entries
        M[off] += 1.0 - M[off].min()
    with tempfile.TemporaryDirectory() as d:
        fi = os.path.join(d, "in.phy")
        fo = os.path.join(d, "out.tre")
        with open(fi, "w") as fh:
            fh.write(f"{n}\n")
            for a in range(n):
                fh.write(f"T{a} " + " ".join(f"{x:.8f}" for x in M[a]) + "\n")
        subprocess.run([FASTME, "-i", fi, "-o", fo, "-m", "B", "-n", "-s", "-T", "1"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True, cwd=d)
        nwk = open(fo).read().strip()
    t = ts.read_tree_newick(nwk)
    for lf in t.traverse_leaves():
        lf.label = taxa[int(lf.label[1:])]
    for nd in t.traverse_preorder():
        nd.edge_length = None
    return t.newick()


# ---------------------------------------------------------------- comparison
def bipartitions(nwk, taxa_subset=None):
    t = parse(nwk)
    leaves = [l.label for l in t.traverse_leaves()]
    if taxa_subset is not None:
        keep = set(taxa_subset)
        if set(leaves) != keep:
            t = t.extract_tree_with(keep)
            t.suppress_unifurcations()
            leaves = [l.label for l in t.traverse_leaves()]
    allset = frozenset(leaves)
    ref = min(allset)
    out = set()
    desc = {}
    for nd in t.traverse_postorder():
        if nd.is_leaf():
            desc[nd] = frozenset([nd.label])
        else:
            s = frozenset().union(*(desc[c] for c in nd.children))
            desc[nd] = s
            if 1 < len(s) < len(allset) - 1:
                out.add(s if ref not in s else allset - s)
    return out, allset


def fn_fp(true_nwk, est_nwk):
    bt, L = bipartitions(true_nwk)
    be, L2 = bipartitions(est_nwk, L)
    nint = len(L) - 3
    fn = len(bt - be) / nint
    fp = len(be - bt) / max(1, nint)
    return fn, fp, bt - be


ASTRAL = os.environ.get("ASTRAL", "astral")


def astral(trees, extra=()):
    """ASTRAL (ASTER implementation, quartet score, no branch lengths needed)."""
    with tempfile.TemporaryDirectory() as d:
        fi = os.path.join(d, "g.tre")
        fo = os.path.join(d, "s.tre")
        with open(fi, "w") as fh:
            fh.write("\n".join(trees) + "\n")
        subprocess.run([ASTRAL, "-i", fi, "-o", fo, "-t", "1", "-u", "0", *extra],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        return open(fo).read().strip()


class Accumulator:
    """Running sums for nested gene sets (so every k in a grid costs one pass)."""

    def __init__(self, taxa):
        self.taxa = list(taxa)
        self.ti = {t: k for k, t in enumerate(taxa)}
        n = len(taxa)
        self.S = np.zeros((n, n))
        self.C = np.zeros((n, n))

    def add(self, nwk, **kw):
        tix, D = gene_contrib(nwk, self.ti, **kw)
        self.S[np.ix_(tix, tix)] += D
        self.C[np.ix_(tix, tix)] += 1

    def matrix(self):
        with np.errstate(invalid="ignore", divide="ignore"):
            M = self.S / self.C
        np.fill_diagonal(M, 0.0)
        return M


ASTEROID = os.environ.get("ASTEROID", "/opt/tools/Asteroid/build/bin/asteroid")


def asteroid(trees, correction=True):
    """Asteroid (Morel, Williams & Stamatakis 2023): missing-data-aware internode distance."""
    with tempfile.TemporaryDirectory() as d:
        fi = os.path.join(d, "g.tre")
        with open(fi, "w") as fh:
            fh.write("\n".join(trees) + "\n")
        cmd = [ASTEROID, "-i", fi, "-p", os.path.join(d, "o"), "-r", "0"]
        if not correction:
            cmd.append("-n")
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True, cwd=d)
        return open(os.path.join(d, "o.bestTree.newick")).read().strip()


def strip(nwk):
    t = parse(nwk)
    for nd in t.traverse_preorder():
        nd.edge_length = None
        if not nd.is_leaf():
            nd.label = None
    out = t.newick()
    return out.split("] ", 1)[1] if out.startswith("[&") else out


def _adj_newick(adj, lab, start):
    def rec(u, p):
        ch = [w for w in adj[u] if w != p]
        if not ch:
            return lab[u]
        return "(" + ",".join(rec(w, u) for w in ch) + ")"
    return rec(start, -1) + ";"


def complete_gene(nwk, ref, taxa, root_leaf):
    """Greedy insertion of the taxa missing from a gene tree, guided by a reference species tree
    `ref` (treeswift, rooted at leaf `root_leaf`). Each missing x is attached as sibling of the
    gene-tree LCA (rooted at an outside anchor) of the present taxa in the smallest reference
    clade around x that contains present taxa ("OCTAL-lite")."""
    labels, adj0, leaves = to_adj(parse(nwk))
    adj = {u: list(v) for u, v in enumerate(adj0)}
    lab = {lf: l for lf, l in zip(leaves, labels)}
    leafof = {l: lf for lf, l in lab.items()}
    # remove degree-2 nodes (rooted input)
    for u in list(adj):
        if len(adj[u]) == 2 and u not in lab:
            a, b = adj[u]
            adj[a] = [b if w == u else w for w in adj[a]]
            adj[b] = [a if w == u else w for w in adj[b]]
            del adj[u]
    refleaf = {l.label: l for l in ref.traverse_leaves()}
    refclade = {}
    for nd in ref.traverse_postorder():
        refclade[nd] = {nd.label} if nd.is_leaf() else set().union(*(refclade[c] for c in nd.children))
    nxt = max(adj) + 1
    for x in taxa:
        present = set(leafof)
        if x in present:
            continue
        node = refleaf[x]
        P = set()
        while node.parent is not None:
            node = node.parent
            P = (refclade[node] & present) - {root_leaf}
            if P:
                break
        out = [y for y in present if y not in P]
        if P and not out:
            # every present taxon is inside the clade: split the x-free part of it instead
            nd = node
            while True:
                kids = [c for c in nd.children if (refclade[c] & present) and x not in refclade[c]]
                xk = [c for c in nd.children if x in refclade[c]]
                if len(kids) >= 2 or (len(kids) == 1 and xk and (refclade[xk[0]] & present)):
                    break
                if len(kids) == 1 and (not xk or not (refclade[xk[0]] & present)):
                    nd = kids[0]
                    continue
                break
            if len(kids) >= 2:
                P = refclade[kids[0]] & present
            out = [y for y in present if y not in P]
        if not P or not out:
            a = leafof[sorted(present)[0]]
            m, par = a, adj[a][0]
        else:
            anchor = leafof[root_leaf] if root_leaf in out else leafof[sorted(out)[0]]
            parent = {anchor: -1}
            order = [anchor]
            for u in order:
                for w in adj[u]:
                    if w not in parent:
                        parent[w] = u
                        order.append(w)
            cnt = {}
            for u in reversed(order):
                c = 1 if (u in lab and lab[u] in P) else 0
                cnt[u] = c + sum(cnt[w] for w in adj[u] if parent.get(w) == u)
            # LCA: deepest node with cnt == |P|
            m = next(u for u in reversed(order) if cnt[u] == len(P))
            par = parent[m]
        u, v = nxt, nxt + 1
        nxt += 2
        adj[par] = [u if w == m else w for w in adj[par]]
        adj[m] = [u if w == par else w for w in adj[m]]
        adj[u] = [par, m, v]
        adj[v] = [u]
        lab[v] = x
        leafof[x] = v
    start = next(u for u in adj if u not in lab)
    return _adj_newick(adj, lab, start)
