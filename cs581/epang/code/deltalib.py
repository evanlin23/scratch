"""Delta error and placement mapping on a fixed backbone tree.

Delta error of query q placed on backbone edge e (as in SCAMPP/BSCAMPP):
    FN(T + q@e, T*|B+q) - FN(T, T*|B)
where T is the estimated backbone tree, T* the true tree, B the backbone leaves,
FN = number of non-trivial true bipartitions missing from the estimate.

Bipartitions are hashed (xor of random 64-bit leaf keys) with every tree rooted at a
fixed backbone leaf r, so each bipartition is represented by its side without r.
"""
import json, random, re
import treeswift


def leaf_keys(names, seed=12345):
    rng = random.Random(seed)
    return {n: rng.getrandbits(64) for n in sorted(names)}


def postorder_hash(tree, key, only=None):
    """For each node, (hash, count) of leaves below restricted to `only` (None = all keyed)."""
    H, C = {}, {}
    for v in tree.traverse_postorder():
        if v.is_leaf():
            lab = v.get_label()
            if lab in key and (only is None or lab in only):
                H[v], C[v] = key[lab], 1
            else:
                H[v], C[v] = 0, 0
        else:
            h, c = 0, 0
            for ch in v.children:
                h ^= H[ch]; c += C[ch]
            H[v], C[v] = h, c
    return H, C


class Node:
    __slots__ = ('label', 'parent', 'children', 'length')

    def __init__(self, label=None):
        self.label, self.parent, self.children, self.length = label, None, [], 0.0

    def get_label(self):
        return self.label

    def get_edge_length(self):
        return self.length

    def is_leaf(self):
        return not self.children


class RTree:
    """Tree re-rooted at leaf r (r is the root; its single child is r's neighbour).
    Degree-2 nodes of the input are suppressed. Every non-root node v stands for the
    edge above it, so the edges of the unrooted tree are exactly the non-root nodes."""

    def __init__(self, newick_path, r):
        t = treeswift.read_tree_newick(newick_path)
        adj = {}
        def add(a, b, l):
            adj.setdefault(a, {})[b] = l; adj.setdefault(b, {})[a] = l
        for v in t.traverse_preorder():
            adj.setdefault(v, {})
            if v.parent is not None:
                add(v, v.parent, v.edge_length or 0.0)
        # suppress degree-2 nodes
        for v in list(adj):
            if len(adj[v]) == 2 and not v.is_leaf():
                (a, la), (b, lb) = adj[v].items()
                del adj[a][v], adj[b][v]; del adj[v]
                add(a, b, la + lb)
        start = [v for v in adj if v.is_leaf() and v.label == r][0]
        mk = {start: Node(r)}
        self.root = mk[start]
        stack = [start]
        while stack:
            u = stack.pop()
            for w, l in adj[u].items():
                if w in mk:
                    continue
                n = Node(w.label if w.is_leaf() else None)
                n.parent, n.length = mk[u], l
                mk[u].children.append(n); mk[w] = n
                stack.append(w)

    def traverse_postorder(self):
        out, stack = [], [(self.root, False)]
        while stack:
            v, done = stack.pop()
            if done:
                out.append(v)
            else:
                stack.append((v, True))
                for c in v.children:
                    stack.append((c, False))
        return out

    def traverse_leaves(self):
        return [v for v in self.traverse_postorder() if v.is_leaf()]


class Backbone:
    """Estimated backbone tree T rooted at leaf r, with true tree for scoring."""

    def __init__(self, est_tree_path, true_tree_path, backbone_names, root_leaf=None):
        self.B = set(backbone_names)
        self.r = root_leaf or sorted(self.B)[0]
        self.T = RTree(est_tree_path, self.r)
        self.Ts = RTree(true_tree_path, self.r)
        alln = [v.get_label() for v in self.Ts.traverse_leaves()] + [self.r]
        self.key = leaf_keys(alln)
        self.nB = len(self.B)
        # estimated tree
        self.H, self.C = postorder_hash(self.T, self.key, self.B)
        self.nodes = [v for v in self.T.traverse_postorder()]
        self.by_hash = {self.H[v]: v for v in self.nodes if v.parent is not None}
        # true tree hashes restricted to backbone leaves
        self.Hs, self.Cs = postorder_hash(self.Ts, self.key, self.B)
        self.leaf_s = {v.get_label(): v for v in self.Ts.traverse_leaves()}
        TB = {self.Hs[v] for v in self.Ts.traverse_postorder() if 2 <= self.Cs[v] <= self.nB - 2}
        est = {self.H[v] for v in self.nodes if v.parent is not None and 2 <= self.C[v] <= self.nB - 2}
        self.TB = TB
        self.FN_T = len(TB - est)
        self.RF_note = (len(TB), len(est), self.FN_T)
        self._q = {}

    def true_set_q(self, q):
        if q in self._q:
            return self._q[q]
        hq = self.key[q]
        anc = set()
        v = self.leaf_s[q].parent
        while v is not None:
            anc.add(v); v = v.parent
        N = self.nB + 1
        S = set()
        for v in self.Ts.traverse_postorder():
            if v in anc:
                h, c = self.Hs[v] ^ hq, self.Cs[v] + 1
            else:
                h, c = self.Hs[v], self.Cs[v]
            if 2 <= c <= N - 2:
                S.add(h)
        # membership arrays on the estimated tree
        hq = self.key[q]
        in0 = {v: (self.H[v] in S) for v in self.nodes}
        in1 = {v: ((self.H[v] ^ hq) in S) for v in self.nodes}
        base = sum(in0[v] for v in self.nodes if v.parent is not None and 2 <= self.C[v])
        res = (S, in0, in1, base)
        self._q[q] = res
        return res

    def delta(self, q, c):
        """delta error for q inserted on edge above node c of T (rooted at r)."""
        S, in0, in1, base = self.true_set_q(q)
        N = self.nB + 1
        m = base
        v = c.parent
        while v.parent is not None:  # the root (leaf r) carries no edge
            if 2 <= self.C[v]:
                m -= in0[v]
            if 2 <= self.C[v] + 1 <= N - 2:
                m += in1[v]
            v = v.parent
        if 2 <= self.C[c] + 1 <= N - 2:
            m += in1[c]
        fn = len(S) - m
        return fn - self.FN_T

    def node_dist(self, a, b):
        """number of edges between edges above nodes a and b (0 if same)."""
        if a is b:
            return 0
        da, db = {}, {}
        v, d = a, 0
        while v is not None:
            da[v] = d; v = v.parent; d += 1
        v, d = b, 0
        while v is not None:
            if v in da:
                return da[v] + d
            v = v.parent; d += 1
        return None

    # ---- mapping placements on (sub)trees onto T ----
    def map_jplace(self, jplace_path):
        """Return {query: (T node c, lwr)} for the best placement of each query.
        Works for jplace trees whose leaf set is any subset of B (edges of the
        subtree correspond to paths in T; the position along the path follows
        distal_length proportionally)."""
        J = json.load(open(jplace_path))
        tstr = J['tree']
        t, emap = parse_edge_tree(tstr)
        leaves = {v.get_label() for v in t.traverse_leaves()}
        S = leaves & self.B
        hS = 0
        for l in S:
            hS ^= self.key[l]
        Hj, Cj = postorder_hash(t, self.key, S)
        # T node hashes restricted to S
        HT, CT = postorder_hash(self.T, self.key, S)
        chains = {}
        for v in self.nodes:
            if v.parent is None:
                continue
            chains.setdefault(HT[v], []).append(v)
        fields = J['fields']
        fi = {f: i for i, f in enumerate(fields)}
        out = {}
        for pl in J['placements']:
            ps = pl['p']
            best = max(ps, key=lambda p: p[fi['like_weight_ratio']])
            en = int(best[fi['edge_num']])
            node = emap[en]
            hA = Hj[node]
            dist = best[fi['distal_length']]
            blen = node.get_edge_length() or 0.0
            frac = 0.0 if blen <= 0 else min(max(dist / blen, 0.0), 1.0)
            if hA in chains:
                ch = chains[hA]; from_bottom = True
            else:
                ch = chains.get(hS ^ hA); from_bottom = False
            if not ch:
                raise RuntimeError('cannot map edge')
            # order chain bottom -> top (postorder list already bottom-up)
            ch = sorted(ch, key=lambda v: self.C[v])
            if not from_bottom:
                ch = ch[::-1]
            lens = [max(v.get_edge_length() or 0.0, 0.0) for v in ch]
            tot = sum(lens)
            pos = frac * tot
            c = ch[-1]
            acc = 0.0
            for v, l in zip(ch, lens):
                acc += l
                if pos <= acc:
                    c = v; break
            names = pl.get('n') or [x[0] for x in pl.get('nm', [])]
            for qn in names:
                out[qn] = (c, best[fi['like_weight_ratio']])
        return out


def parse_edge_tree(s):
    """Parse a jplace newick string with {edge} tags. Returns (treeswift tree, {edge_num: node})."""
    # tag form: label:len{n}  -> we move {n} into label position by inserting marker before ':'
    s2 = re.sub(r'([^(),:;]*):([0-9eE.+-]+)\{(\d+)\}', lambda m: f"{m.group(1)}@@{m.group(3)}:{m.group(2)}", s)
    t = treeswift.read_tree_newick(s2)
    emap = {}
    for v in t.traverse_postorder():
        lab = v.get_label()
        if lab and '@@' in lab:
            base, en = lab.rsplit('@@', 1)
            v.set_label(base if base else None)
            emap[int(en)] = v
    return t, emap
