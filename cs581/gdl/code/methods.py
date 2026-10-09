"""Distance-based species-tree methods for multi-copy gene-family trees.

* root_and_tag      : ASTRAL-Pro/DISCO-style rooting (min #duplications, ties ->
                      min #losses) and tagging (duplication iff children's species sets overlap)
* gene_distances    : per-gene species x species matrices for
                      'multi'     ASTRID-multi: all copy pairs, all internal nodes (as ASTRID's imap->average)
                      'pro'       ASTRID-Pro: orthologous pairs only (LCA = speciation), speciation nodes only
                      'ortho_all' ablation: orthologous pairs, all nodes
                      'spec_all'  ablation: all pairs, speciation nodes only
                      agg='mean' (per-gene average over pairs) or 'min' (closest copy, STAG-like)
* disco_decompose   : DISCO-style decomposition (re-implementation; original repo unreachable)
* fastme_tree       : species tree from an averaged matrix (FastME 2, BalME + NNI + SPR)
"""
import os
import subprocess
import tempfile

import numpy as np

from phylo import Tree, parse_newick

FASTME = os.environ.get("FASTME", "/opt/mm/root/envs/gdl/bin/fastme")
ASTRALPRO = os.environ.get("ASTRALPRO", "/opt/mm/root/envs/gdl/bin/astral-pro")


def simphy_species(label):
    return label.split("_")[0]


# ---------------------------------------------------------------- rooting and tagging

def _unrooted_adj(t):
    nb = [[] for _ in t.parent]
    for v, p in enumerate(t.parent):
        if p >= 0:
            nb[v].append(p)
            nb[p].append(v)
    r = t.root
    if len(nb[r]) == 2:  # suppress degree-2 root
        a, b = nb[r]
        nb[a].remove(r); nb[b].remove(r); nb[a].append(b); nb[b].append(a)
        nb[r] = []
    return nb


def _overlap(sets):
    acc = 0
    for s in sets:
        if acc & s:
            return True
        acc |= s
    return False


def root_and_tag(t, sp_of, sp_index, keep_root=False, truetags=False):
    """Returns (rooted Tree, tags list with True for duplication, leaf species idx list).
    truetags: keep the given root and use internal labels 'D' (simulator output)."""
    if truetags:
        tags, leafsp = tag(t, sp_of, sp_index)
        tags = [bool(t.children[v]) and t.label[v] == "D" for v in range(len(t.parent))]
        return t, tags, leafsp
    if keep_root:
        rt = t
    else:
        nb = _unrooted_adj(t)
        leaves = [v for v in range(len(nb)) if len(nb[v]) == 1]
        internal = [v for v in range(len(nb)) if len(nb[v]) >= 2]
        if not internal:
            rt = t
        else:
            r0 = internal[0]
            par = {r0: -1}
            order = [r0]
            for v in order:
                for q in nb[v]:
                    if q not in par:
                        par[q] = v
                        order.append(q)
            ch = {v: [q for q in nb[v] if q != par[v]] for v in order}
            D = {}
            for v in reversed(order):
                if not ch[v]:
                    D[v] = 1 << sp_index[sp_of(t.label[v])]
                else:
                    m = 0
                    for c in ch[v]:
                        m |= D[c]
                    D[v] = m
            U = {r0: 0}
            for v in order:
                cs = ch[v]
                for c in cs:
                    m = U[v]
                    for c2 in cs:
                        if c2 != c:
                            m |= D[c2]
                    U[c] = m

            def side(v, q):  # species set on q's side of edge v-q
                return D[q] if par.get(q) == v else U[v]

            def dup_given_parent(v, p):
                if not ch[v] and par[v] != -1 or len(nb[v]) == 1:
                    return 0
                return 1 if _overlap([side(v, q) for q in nb[v] if q != p]) else 0

            base = 0
            for w in order:
                if w != r0 and len(nb[w]) > 1:
                    base += dup_given_parent(w, par[w])
            g = {}
            for c in ch[r0]:
                g[c] = base + dup_given_parent(r0, c)
            for v in order[1:]:
                for c in ch[v]:
                    g[c] = g[v] - dup_given_parent(v, par[v]) + dup_given_parent(v, c)
            score = {v: g[v] + (1 if D[v] & U[v] else 0) for v in order[1:]}
            best = min(score.values())
            cands = [v for v in order[1:] if score[v] == best]
            if len(cands) > 1:
                cands = cands[:40]
                lossc = [(_loss_score(_build_rooted(t, nb, v, par[v]), sp_of, sp_index), i)
                         for i, v in enumerate(cands)]
                v = cands[min(lossc)[1]]
            else:
                v = cands[0]
            rt = _build_rooted(t, nb, v, par[v])
    tags, leafsp = tag(rt, sp_of, sp_index)
    return rt, tags, leafsp


def _build_rooted(t, nb, a, b):
    rt = Tree()
    root = rt.add(-1)
    rt.root = root
    stack = [(a, b, root), (b, a, root)]
    while stack:
        v, frm, p = stack.pop()
        i = rt.add(p)
        rt.label[i] = t.label[v] if len(nb[v]) == 1 else None
        for q in nb[v]:
            if q != frm:
                stack.append((q, v, i))
    return rt


def tag(rt, sp_of, sp_index):
    n = len(rt.parent)
    S = [0] * n
    tags = [False] * n
    leafsp = [-1] * n
    for v in rt.postorder():
        if not rt.children[v]:
            leafsp[v] = sp_index[sp_of(rt.label[v])]
            S[v] = 1 << leafsp[v]
        else:
            sets = [S[c] for c in rt.children[v]]
            tags[v] = _overlap(sets)
            m = 0
            for s in sets:
                m |= s
            S[v] = m
    rt._S = S
    return tags, leafsp




def _loss_score(rt, sp_of, sp_index):
    tags, _ = tag(rt, sp_of, sp_index)
    S = rt._S
    pc = lambda m: bin(m).count("1")
    tot = 0
    for v in range(len(rt.parent)):
        if tags[v]:
            for c in rt.children[v]:
                tot += pc(S[v]) - pc(S[c])
    return tot


# ---------------------------------------------------------------- distances

def gene_distances(rt, tags, leafsp, nsp, mode="pro", agg="mean"):
    """Per-gene species distance: returns (sum matrix, count matrix) for this gene,
    where for agg='mean' entries are (avg, 1) and for 'min' (min, 1)."""
    n = len(rt.parent)
    root = rt.root
    if mode in ("multi", "ortho_all"):
        counted = [True] * n
    else:
        counted = [not tg for tg in tags]
    for v in range(n):
        if not rt.children[v]:
            counted[v] = False
    if len(rt.children[root]) < 3:
        counted[root] = False  # degree-2 root is not a node of the unrooted tree
    need_ortho = mode in ("pro", "ortho_all")
    # cum[v] = number of counted nodes on root..v inclusive
    cum = [0] * n
    order = rt.postorder()[::-1]
    for v in order:
        p = rt.parent[v]
        cum[v] = (cum[p] if p >= 0 else 0) + (1 if counted[v] else 0)
    leaves_under = [None] * n
    tot = np.zeros((nsp, nsp))
    cnt = np.zeros((nsp, nsp))
    mn = np.full((nsp, nsp), np.inf)
    for v in rt.postorder():
        ch = rt.children[v]
        if not ch:
            leaves_under[v] = np.array([v])
            continue
        L = [leaves_under[c] for c in ch]
        leaves_under[v] = np.concatenate(L)
        if need_ortho and tags[v]:
            continue
        for i in range(len(ch)):
            for j in range(i + 1, len(ch)):
                a, b = L[i], L[j]
                ca = np.array([cum[rt.parent[x]] for x in a])
                cb = np.array([cum[rt.parent[x]] for x in b])
                d = ca[:, None] + cb[None, :] - 2 * cum[v] + (1 if counted[v] else 0)
                sa = np.array([leafsp[x] for x in a])
                sb = np.array([leafsp[x] for x in b])
                SA = np.broadcast_to(sa[:, None], d.shape).ravel()
                SB = np.broadcast_to(sb[None, :], d.shape).ravel()
                dd = d.ravel().astype(float)
                keep = SA != SB
                SA, SB, dd = SA[keep], SB[keep], dd[keep]
                if agg == "mean":
                    np.add.at(tot, (SA, SB), dd)
                    np.add.at(cnt, (SA, SB), 1)
                else:
                    np.minimum.at(mn, (SA, SB), dd)
    if agg == "mean":
        tot = tot + tot.T
        cnt = cnt + cnt.T
        has = cnt > 0
        out = np.zeros((nsp, nsp))
        out[has] = tot[has] / cnt[has]
        return out, has.astype(float)
    mn = np.minimum(mn, mn.T)
    has = np.isfinite(mn)
    out = np.where(has, mn, 0.0)
    return out, has.astype(float)


def average_matrix(per_gene):
    tot = cnt = None
    for m, h in per_gene:
        if tot is None:
            tot, cnt = m.copy(), h.copy()
        else:
            tot += m
            cnt += h
    D = np.zeros_like(tot)
    has = cnt > 0
    D[has] = tot[has] / cnt[has]
    np.fill_diagonal(D, 0)
    nmiss = int(((~has).sum() - D.shape[0]) // 2)
    if nmiss > 0:
        mx = D[has].max() if has.any() else 1.0
        D[~has] = mx  # simple fill-in for pairs never observed together
        np.fill_diagonal(D, 0)
    return D, nmiss


# ---------------------------------------------------------------- DISCO decomposition

def disco_decompose(rt, tags, sp_of, min_leaves=4):
    """Post-order: at each duplication node keep the largest child subtree and split
    the remaining children off as separate trees (leaf sets partition the input)."""
    out = []
    S = rt._S
    import sys
    sys.setrecursionlimit(100000)

    def rec(v):
        ch = rt.children[v]
        if not ch:
            return sp_of(rt.label[v]), 1, 1 << 0
        parts = [(rec(c), S[c]) for c in ch]
        if not tags[v]:
            nw = "(" + ",".join(p[0][0] for p in parts) + ")"
            return nw, sum(p[0][1] for p in parts), 0
        # duplication: keep children greedily by size while species-disjoint
        parts.sort(key=lambda p: -p[0][1])
        kept, acc, k = [], 0, 0
        for (nw, sz, _), s in parts:
            if not (acc & s):
                kept.append((nw, sz)); acc |= s
            else:
                out.append((nw, sz))
        if len(kept) == 1:
            return kept[0][0], kept[0][1], 0
        return "(" + ",".join(k[0] for k in kept) + ")", sum(k[1] for k in kept), 0

    nw, sz, _ = rec(rt.root)
    out.append((nw, sz))
    return [w + ";" for w, s in out if s >= min_leaves]


# ---------------------------------------------------------------- tree building

def fastme_tree(D, names, workdir=None):
    with tempfile.TemporaryDirectory(dir=workdir) as td:
        fin, fout = os.path.join(td, "d.phy"), os.path.join(td, "t.nwk")
        with open(fin, "w") as f:
            f.write("%d\n" % len(names))
            for i, nm in enumerate(names):
                f.write(nm + " " + " ".join("%.6f" % x for x in D[i]) + "\n")
        subprocess.run([FASTME, "-i", fin, "-o", fout, "-m", "B", "-n", "B", "-s", "-T", "1"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        with open(fout) as f:
            return parse_newick(f.read())


def astral_pro(gene_tree_file, mapping_file, out_file, threads=1):
    subprocess.run([ASTRALPRO, "-i", gene_tree_file, "-a", mapping_file, "-o", out_file,
                    "-t", str(threads)], check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL)
    with open(out_file) as f:
        return parse_newick(f.read())


# ---------------------------------------------------------------- pipeline helper

def species_tree_from_genes(gene_trees, sp_of, species, mode="pro", agg="mean",
                            keep_root=False, tagged=None, truetags=False):
    """gene_trees: list of Tree. tagged: optional precomputed list of (rt, tags, leafsp)."""
    sp_index = {s: i for i, s in enumerate(species)}
    if tagged is None:
        tagged = [root_and_tag(t, sp_of, sp_index, keep_root=keep_root, truetags=truetags) for t in gene_trees]
    per = [gene_distances(rt, tg, ls, len(species), mode=mode, agg=agg) for rt, tg, ls in tagged]
    D, nmiss = average_matrix(per)
    return fastme_tree(D, species), nmiss
