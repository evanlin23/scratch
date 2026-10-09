"""Pure gene duplication-and-loss simulator (no ILS): gene-family tree = locus tree.

Species tree: rooted, with branch lengths. Each branch e may carry its own
duplication rate lam[e] and loss rate mu[e] (birth-death on copies).
A family starts as one copy at the top of a root branch of length root_len.
Output: Newick with leaves 'species_copyid' and internal labels 'D'/'S'
(true duplication / speciation). Unifurcations created by loss are suppressed.
"""
import random

from phylo import Tree, parse_newick


class GNode:
    __slots__ = ("kind", "ch", "label")

    def __init__(self, kind, ch=None, label=None):
        self.kind, self.ch, self.label = kind, ch or [], label


def _evolve_edge(t, lam, mu, rng, cap):
    """One copy along an edge of length t. Returns list of (tip placeholder) via
    a nested structure: a tree of 'D' nodes whose open tips are None (to be filled)."""
    # returns a builder: function that maps each surviving tip to the subtree below
    tot = lam + mu
    if tot <= 0:
        return "TIP"
    # Gillespie on one lineage, recursive
    def rec(rem):
        while True:
            dt = rng.expovariate(tot)
            if dt >= rem:
                return "TIP"
            rem -= dt
            if rng.random() < lam / tot:
                if cap[0] <= 0:
                    raise OverflowError
                cap[0] -= 1
                return ("D", rec(rem), rec(rem))
            return None  # lost
    return rec(t)


def simulate_family(stree, lam, mu, rng, root_len=0.0, cap=5000):
    """stree: phylo.Tree with lengths; lam/mu: lists indexed by species node id
    (rate on the branch above that node; the root entry is used for the root branch)."""
    capbox = [cap]

    def fill(struct, v):
        # struct from _evolve_edge; replace TIPs with the speciation/leaf at node v
        if struct is None:
            return None
        if struct == "TIP":
            return at_node(v)
        _, a, b = struct
        a, b = fill(a, v), fill(b, v)
        if a is None:
            return b
        if b is None:
            return a
        return GNode("D", [a, b])

    counter = {}

    def at_node(v):
        ch = stree.children[v]
        if not ch:
            sp = stree.label[v]
            k = counter.get(sp, 0)
            counter[sp] = k + 1
            return GNode("L", label="%s_%d" % (sp, k))
        kids = []
        for c in ch:
            s = _evolve_edge(stree.length[c] or 0.0, lam[c], mu[c], rng, capbox)
            g = fill(s, c)
            if g is not None:
                kids.append(g)
        if not kids:
            return None
        if len(kids) == 1:
            return kids[0]
        return GNode("S", kids)

    r = stree.root
    s = _evolve_edge(root_len, lam[r], mu[r], rng, capbox) if root_len > 0 else "TIP"
    try:
        return fill(s, r)
    except OverflowError:
        return "OVERFLOW"


def to_newick(g):
    def rec(x):
        if x.kind == "L":
            return x.label
        return "(" + ",".join(rec(c) for c in x.ch) + ")" + x.kind
    return rec(g) + ";"


def leaf_species(g):
    out, st = set(), [g]
    while st:
        x = st.pop()
        if x.kind == "L":
            out.add(x.label.split("_")[0])
        else:
            st.extend(x.ch)
    return out


def simulate(stree, lam, mu, nfam, seed, min_species=4, root_len=0.0, cap=5000):
    rng = random.Random(seed)
    fams, tries, over = [], 0, 0
    while len(fams) < nfam:
        tries += 1
        g = simulate_family(stree, lam, mu, rng, root_len=root_len, cap=cap)
        if g == "OVERFLOW":
            over += 1
            continue
        if g is None or g.kind == "L":
            continue
        if len(leaf_species(g)) < min_species:
            continue
        fams.append(to_newick(g))
    return fams, tries, over


def yule_tree(n, seed, height=1.0):
    """Ultrametric Yule tree with n leaves 'T0'..'T{n-1}', scaled to the given height."""
    rng = random.Random(seed)
    lineages = [[None, [], "T%d" % i, 0.0] for i in range(n)]  # node: [parent, children, label, height]
    nodes = list(lineages)
    t = 0.0
    while len(lineages) > 1:
        k = len(lineages)
        t += rng.expovariate(k)
        i, j = rng.sample(range(k), 2)
        a, b = lineages[i], lineages[j]
        p = [None, [a, b], None, t]
        a[0] = b[0] = p
        nodes.append(p)
        lineages = [x for idx, x in enumerate(lineages) if idx not in (i, j)] + [p]
    root = lineages[0]
    scale = height / root[3]
    def rec(x):
        if not x[1]:
            s = x[2]
        else:
            s = "(" + ",".join(rec(c) for c in x[1]) + ")"
        if x[0] is not None:
            s += ":%.6f" % ((x[0][3] - x[3]) * scale)
        return s
    return rec(root) + ";"
