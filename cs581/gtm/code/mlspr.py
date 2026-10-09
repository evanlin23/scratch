"""GTM-Blend-ML: constrained SPR hill-climbing scored by maximum likelihood.

Start from the GTM tree. Each round:
  1. RAxML-NG optimizes branch lengths of the current tree (fixed GTR+G model).
  2. Every SPR move within `radius` that keeps all subset trees induced (feasibility
     test in blend.py) is scored by a lazy-SPR likelihood (GTR+G Felsenstein pruning
     with directional partial likelihoods; only the 3 branches at the regraft point
     are adjusted on a small grid; all other branch lengths fixed).
  3. The top-K candidates are re-scored by RAxML-NG with full branch-length
     optimization; the best is accepted if it beats the current tree by > eps.
Stops when no candidate improves.
"""
import os
import re
import subprocess

import numpy as np
from scipy.stats import gamma as gamma_dist

from blend import BlendSearch
from phylo import parse_newick

BIN = os.environ.get("BIOBIN", "/opt/mm/root/envs/bio/bin")
TIP = np.zeros((16, 4))
for code in range(16):
    for k in range(4):
        if code >> k & 1:
            TIP[code, k] = 1.0
TIP[0] = 1.0


def split_lengths(newick, idx):
    """{normalized split bitset (bit 0 excluded): length} for a rooted/unrooted
    Newick with branch lengths; root edges are merged."""
    s = re.sub(r"\[[^\]]*\]", "", newick.strip())
    toks = re.findall(r"\(|\)|,|;|:[^(),;]+|[^(),;:]+", s)
    stack = [[]]
    items = []  # (bitset, length)
    last = None
    full = (1 << len(idx)) - 1
    for t in toks:
        if t == "(":
            stack.append([])
            last = None
        elif t == ",":
            last = None
        elif t == ")":
            ch = stack.pop()
            b = 0
            for c in ch:
                b |= c[0]
            node = [b, 0.0]
            stack[-1].append(node)
            items.append(node)
            last = node
        elif t.startswith(":"):
            last[1] = float(t[1:])
        elif t == ";":
            break
        elif last is None:
            node = [1 << idx[t.strip()], 0.0]
            stack[-1].append(node)
            items.append(node)
            last = node
    out = {}
    for b, l in items:
        if b == full or b == 0:
            continue
        k = b if not (b & 1) else full ^ b
        out[k] = out.get(k, 0.0) + l
    return out


def parse_model(path):
    txt = open(path).read().split(",")[0].strip()
    rates = [float(x) for x in re.search(r"GTR\{([^}]*)\}", txt).group(1).split("/")]
    m = re.search(r"\+F[CO]?\{([^}]*)\}", txt)
    freqs = [float(x) for x in m.group(1).split("/")] if m else [0.25] * 4
    alpha = float(re.search(r"\+G4m?\{([^}]*)\}", txt).group(1))
    return txt, rates, freqs, alpha


class GTRG:
    def __init__(self, rates, freqs, alpha, ncat=4):
        a, b, c, d, e, f = rates
        pi = np.array(freqs) / sum(freqs)
        R = np.array([[0, a, b, c], [a, 0, d, e], [b, d, 0, f], [c, e, f, 0]], float)
        Q = R * pi[None, :]
        np.fill_diagonal(Q, -Q.sum(1))
        Q /= -(np.diag(Q) * pi).sum()
        # symmetrize for a stable eigendecomposition
        sq = np.sqrt(pi)
        S = sq[:, None] * Q / sq[None, :]
        lam, U = np.linalg.eigh(S)
        self.lam = lam
        self.U = U / sq[:, None]          # left
        self.Ui = (U * sq[:, None]).T      # right
        self.pi = pi
        # mean-of-category discrete gamma (RAxML-NG "G4m")
        qs = gamma_dist.ppf(np.arange(1, ncat) / ncat, alpha, scale=1 / alpha)
        bounds = np.concatenate([[0], qs, [np.inf]])
        from scipy.special import gammainc
        cdf1 = gammainc(alpha + 1, bounds * alpha)
        self.r = (cdf1[1:] - cdf1[:-1]) * ncat
        self.ncat = ncat

    def P(self, t):
        """(G,4,4) transition matrices."""
        t = max(t, 1e-8)
        e = np.exp(self.lam[None, :] * self.r[:, None] * t)  # G x 4
        return np.einsum("ik,gk,kj->gij", self.U, e, self.Ui)


def prop(P, v):
    """v: (G,L,4) partial at child; returns (G,L,4) message at parent."""
    return np.einsum("gij,glj->gli", P, v)


class MLScorer:
    def __init__(self, bs, tipcodes, weights, model, lengths):
        self.bs = bs
        self.w = weights.astype(float)
        self.m = model
        self.tip = TIP[tipcodes]  # (n, L, 4)
        self.len = {}
        norm = lambda b: b if not (b & 1) else bs.ALL ^ b
        for v in range(bs.N):
            if bs.par[v] is not None:
                b = norm(bs.C[v])
                L = lengths.get(b)
                if L is None and bs.C[v].bit_count() == 1:
                    L = lengths.get(norm(bs.C[v]))
                self.len[(v, bs.par[v])] = self.len[(bs.par[v], v)] = L if L is not None else 0.01
        self._pcache = {}
        self.compute()

    def P(self, t):
        k = round(t, 7)
        if k not in self._pcache:
            if len(self._pcache) > 20000:
                self._pcache.clear()
            self._pcache[k] = self.m.P(t)
        return self._pcache[k]

    @staticmethod
    def _scale(v, s):
        mx = v.max(axis=(0, 2))
        mx[mx <= 0] = 1e-300
        return v / mx[None, :, None], s + np.log(mx)

    def leafvec(self, v):
        return np.broadcast_to(self.tip[v][None], (self.m.ncat,) + self.tip[v].shape)

    def compute(self):
        bs = self.bs
        G = self.m.ncat
        L = self.tip.shape[1]
        down, dsc = {}, {}
        for v in reversed(bs.order):
            if v < bs.n:
                down[v] = self.leafvec(v)
                dsc[v] = np.zeros(L)
                continue
            a, b = [u for u in bs.adj[v] if u != bs.par[v]]
            x = prop(self.P(self.len[(v, a)]), down[a]) * prop(self.P(self.len[(v, b)]), down[b])
            down[v], dsc[v] = self._scale(x, dsc[a] + dsc[b])
        up, usc = {}, {}
        r1 = bs.adj[0][0]
        up[r1], usc[r1] = self.leafvec(0), np.zeros(L)
        for v in bs.order[1:]:
            if v < bs.n:
                continue
            a, b = [u for u in bs.adj[v] if u != bs.par[v]]
            fromup = prop(self.P(self.len[(v, bs.par[v])]), up[v])
            up[a], usc[a] = self._scale(fromup * prop(self.P(self.len[(v, b)]), down[b]), usc[v] + dsc[b])
            up[b], usc[b] = self._scale(fromup * prop(self.P(self.len[(v, a)]), down[a]), usc[v] + dsc[a])
        self.down, self.dsc, self.up, self.usc = down, dsc, up, usc
        # total logL at edge (0, r1)
        self.logL = self.join3(None, None, None, pair=(prop(self.P(self.len[(0, r1)]), down[r1]), dsc[r1],
                                                       self.leafvec(0), np.zeros(L)))

    def D(self, u, v):
        """partial (at node v) of v's side of edge (u,v), with log scaler."""
        if self.bs.par[v] == u:
            return self.down[v], self.dsc[v]
        return self.up[u], self.usc[u]

    def join3(self, A, B, C, pair=None):
        """logL from messages already propagated to a common node."""
        if pair is not None:
            x = pair[0] * pair[2]
            s = pair[1] + pair[3]
        else:
            x = A[0] * B[0] * C[0]
            s = A[1] + B[1] + C[1]
        site = np.einsum("glk,k->l", x, self.m.pi) / self.m.ncat
        return float((self.w * (np.log(np.maximum(site, 1e-300)) + s)).sum())

    def score_moves(self, moves, fgrid=(0.25, 0.5, 0.75), xgrid=(0.5, 1.0, 2.0)):
        """moves: list of feasible (p, x, path). Lazy-SPR logL for each: messages are
        propagated along a DFS from the prune point; only feasible targets are scored
        and only branches leading to feasible targets are explored."""
        bs = self.bs
        want = {}
        prefixes = set()
        for p, x, path in moves:
            want[(p, x, tuple(path))] = None
            for k in range(2, len(path) + 1):
                prefixes.add((p, x) + tuple(path[:k]))
        out = []
        groups = {}
        for p, x, path in moves:
            groups.setdefault((p, x, path[1]), True)
        for (p, x, s) in groups:
            DX, sX = self.D(p, x)
            tx = self.len[(p, x)]
            o = [u for u in bs.adj[p] if u not in (x, s)][0]
            Do, so = self.D(p, o)
            M0 = (prop(self.P(self.len[(p, s)] + self.len[(p, o)]), Do), so)
            stack = [(p, s, M0, [p, s])]
            while stack:
                uprev, u, M, path = stack.pop()
                if u < bs.n:
                    continue
                others = [v for v in bs.adj[u] if v != uprev]
                for t in others:
                    np_ = path + [t]
                    if (p, x) + tuple(np_) not in prefixes:
                        continue
                    t2 = others[0] if others[1] == t else others[1]
                    D2, s2 = self.D(u, t2)
                    V = self._scale(M[0] * prop(self.P(self.len[(u, t2)]), D2), M[1] + s2)
                    F, sF = self.D(u, t)
                    l = self.len[(u, t)]
                    if (p, x, tuple(np_)) in want:
                        best = -np.inf
                        for f in fgrid:
                            A = (prop(self.P(f * l), V[0]), V[1])
                            B = (prop(self.P((1 - f) * l), F), sF)
                            for gx in xgrid:
                                C = (prop(self.P(gx * tx), DX), sX)
                                best = max(best, self.join3(A, B, C))
                        out.append((best, (p, x, np_)))
                    stack.append((u, t, (prop(self.P(l), V[0]), V[1]), np_))
        return out


def raxml_eval(trees, msa, model, prefix, threads=1):
    """Evaluate trees (list of newick) under a fixed model with branch-length
    optimization; returns (list of logL, best tree newick with lengths)."""
    tf = prefix + ".in.tre"
    with open(tf, "w") as f:
        f.write("\n".join(trees) + "\n")
    subprocess.run([f"{BIN}/raxml-ng", "--evaluate", "--msa", msa, "--tree", tf, "--model", model,
                    "--opt-model", "off", "--prefix", prefix, "--threads", str(threads), "--redo",
                    "--log", "RESULT"], check=True, capture_output=True)
    log = open(prefix + ".raxml.log").read()
    lls = [float(x) for x in re.findall(r"Tree #\d+, final logLikelihood: (-[\d.]+)", log)]
    if not lls:
        lls = [float(x) for x in re.findall(r"Final LogLikelihood: (-[\d.]+)", log)]
    best = open(prefix + ".raxml.bestTree").read().strip()
    return lls, best


def run(tree, names, seqs, subset_names, workdir, radius=4, topk=8, eps=0.05, threads=1,
        max_rounds=100, log=print, model_tree=None):
    os.makedirs(workdir, exist_ok=True)
    msa = f"{workdir}/full.fa"
    with open(msa, "w") as f:
        for n in names:
            f.write(">%s\n%s\n" % (n, seqs[n]))
    from blend import ENC
    A = np.vstack([ENC[np.frombuffer(seqs[n].encode(), dtype=np.uint8)] for n in names])
    A[A == 15] = 0
    u, w = np.unique(A, axis=1, return_counts=True)
    tipcodes = u.astype(np.int64)  # (n, L)
    states = np.ones((len(names), 1), np.uint8)
    bs = BlendSearch(tree, names, states, np.ones(1, np.int64), subset_names, radius=radius)
    # model on the starting tree
    mp = f"{workdir}/model"
    subprocess.run([f"{BIN}/raxml-ng", "--evaluate", "--msa", msa, "--tree", model_tree or write_tmp(bs, workdir),
                    "--model", "GTR+G", "--prefix", mp, "--threads", str(threads), "--redo", "--log", "ERROR"],
                   check=True, capture_output=True)
    mtxt, rates, freqs, alpha = parse_model(mp + ".raxml.bestModel")
    model = GTRG(rates, freqs, alpha)
    idx = {n: i for i, n in enumerate(names)}
    hist = []
    cur_ll, cur_nwk = raxml_eval([bs.to_newick()], msa, mtxt, f"{workdir}/cur", threads)
    cur_ll = cur_ll[0]
    hist.append(cur_ll)
    log(f"start logL {cur_ll:.2f}")
    for rnd in range(max_rounds):
        sc = MLScorer(bs, tipcodes, w, model, split_lengths(cur_nwk, idx))
        moves = [m for m in bs.oracle_candidates() if bs.feasible(bs.far(m[0], m[1]), bs.path_fars(m[2]))]
        scored = sc.score_moves(moves)
        scored.sort(key=lambda z: -z[0])
        # distinct topologies among top candidates
        top, seen = [], set()
        for ll, mv in scored:
            nw = apply_copy(bs, mv)
            key = frozenset(splits_of(nw, idx))
            if key in seen:
                continue
            seen.add(key)
            top.append((ll, mv, nw))
            if len(top) >= topk:
                break
        if not top or top[0][0] < sc.logL - 50:
            log(f"round {rnd + 1}: no promising candidate (lazy best {top[0][0] if top else None}, cur {sc.logL:.2f})")
            break
        lls, _ = raxml_eval([t[2] for t in top], msa, mtxt, f"{workdir}/cand", threads)
        j = int(np.argmax(lls))
        log(f"round {rnd + 1}: {len(moves)} feasible, lazy cur {sc.logL:.2f} best {top[0][0]:.2f}; "
            f"full best {lls[j]:.2f} vs cur {cur_ll:.2f}")
        if lls[j] <= cur_ll + eps:
            break
        p, x, path = top[j][1]
        bs.apply(p, x, path)
        bs.recompute()
        cur_ll, cur_nwk = raxml_eval([bs.to_newick()], msa, mtxt, f"{workdir}/cur", threads)
        cur_ll = cur_ll[0]
        hist.append(cur_ll)
    return bs, hist


def write_tmp(bs, workdir):
    p = f"{workdir}/start.tre"
    open(p, "w").write(bs.to_newick() + "\n")
    return p


def apply_copy(bs, mv):
    b2 = BlendSearch.__new__(BlendSearch)
    b2.__dict__ = dict(bs.__dict__)
    b2.adj = [list(a) for a in bs.adj]
    b2.apply(*mv)
    return b2.to_newick()


def splits_of(nwk, idx):
    from phylo import bipartitions
    return bipartitions(parse_newick(nwk), idx)
