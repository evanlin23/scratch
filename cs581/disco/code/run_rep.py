"""Run baseline DISCO vs DISCO-R on one simulated replicate; write one JSON line.

python run_rep.py REPDIR GENEFILE OUTDIR [--s0root outgroup|mindl] [--tags] [--iter]

REPDIR holds s_tree.trees, l_trees.trees, g_true.trees, g_<bp>.trees (SimPhy layout of the DISCO and
DISCO+QR Data Bank sets; leaf labels species_locus_individual). GENEFILE is e.g. g_100.trees.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np
import treeswift

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import discor  # noqa: E402
from discor import species_of  # noqa: E402

TOOLS = "/opt/tools"


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def read_trees(path):
    return [line.strip() for line in open(path) if line.strip()]


# ---------------- RF ----------------
def bipartitions(newick, taxa=None):
    t = treeswift.read_tree_newick(newick)
    leaves = sorted(l.label for l in t.traverse_leaves())
    if taxa is not None:
        assert set(leaves) == set(taxa), (len(leaves), len(taxa))
    allset = frozenset(leaves)
    anchor = leaves[0]
    out = set()
    for n in t.traverse_postorder():
        n.ls = frozenset([n.label]) if n.is_leaf() else frozenset().union(*[c.ls for c in n.child_nodes()])
        if not n.is_leaf() and not n.is_root():
            s = n.ls if anchor not in n.ls else allset - n.ls
            if 1 < len(s) < len(allset) - 1:
                out.add(s)
    return out, leaves


def rf(est, true):
    bt, leaves = bipartitions(true)
    be, _ = bipartitions(est, leaves)
    n = len(leaves) - 3
    fn = len(bt - be) / n
    fp = len(be - bt) / n
    return {"rf": (fn + fp) / 2, "fn": fn, "fp": fp}


# ------------- species tree tools -------------
def astral_pro(genes_species_labels, out, exe="astral-pro2"):
    sh(f"{TOOLS}/{exe} -t 1 -u 0 -i {genes_species_labels} -o {out}")
    return open(out).readline().strip()


def astral(genes, out):
    sh(f"{TOOLS}/astral -t 1 -u 0 -i {genes} -o {out}")
    return open(out).readline().strip()


def astrid(genes, out):
    sh(f"{TOOLS}/wastrid --preset vanilla -i {genes} -o {out}")
    return open(out).readline().strip()


# ------------- rooting S0 -------------
def root_outgroup(newick, og):
    t = treeswift.read_tree_newick(newick)
    leaf = [l for l in t.traverse_leaves() if l.label == og][0]
    discor.disco.reroot_on_edge(t, leaf)
    return t.newick()


def root_mindl(newick, gene_newicks, nsample=200, seed=1):
    """Root unrooted S0 on the edge minimising total DL cost of a sample of gene trees."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(gene_newicks), min(nsample, len(gene_newicks)), replace=False)
    sample = [gene_newicks[i] for i in idx]
    t = treeswift.read_tree_newick(newick)
    cands = [n for n in t.traverse_preorder() if not n.is_root()]
    best, bestc = None, float("inf")
    for k in range(len(cands)):
        tt = treeswift.read_tree_newick(newick)
        node = [n for n in tt.traverse_preorder() if not n.is_root()][k]
        discor.disco.reroot_on_edge(tt, node)
        tt.suppress_unifurcations()
        S = discor.SpeciesTree(tt.newick())
        c = 0.0
        for g in sample:
            gt = treeswift.read_tree_newick(g)
            if gt.root.num_children() == 0:
                continue
            c += discor.best_dl_root(gt, S)[1]
        if c < bestc:
            best, bestc = tt.newick(), c
    return best


# ------------- tagging truth -------------
def heights(t):
    """height = max distance down to an extant (non 'Lost') leaf; None if no extant descendant."""
    for n in t.traverse_postorder():
        if n.is_leaf():
            n.h = None if n.label.startswith("Lost") else 0.0
        else:
            hs = [c.h + (c.edge_length or 0.0) for c in n.child_nodes() if c.h is not None]
            n.h = max(hs) if hs else None


def true_locus_tags(ltree_nw, stree_heights):
    """Return (locus leaf label -> index, matrix ortho[i,j] True iff locus LCA is a speciation)."""
    t = treeswift.read_tree_newick(ltree_nw)
    heights(t)
    labels = [l.label for l in t.traverse_leaves() if not l.label.startswith("Lost")]
    ix = {l: i for i, l in enumerate(labels)}
    M = np.zeros((len(labels), len(labels)), dtype=bool)
    for n in t.traverse_postorder():
        if n.is_leaf():
            n.L = [ix[n.label]] if n.label in ix else []
            continue
        kids = n.child_nodes()
        n.L = sum((c.L for c in kids), [])
        if n.h is None:
            continue
        spec = any(abs(n.h - h) <= 1e-6 * max(1.0, h) for h in stree_heights)
        for a in range(len(kids)):
            for b in range(a + 1, len(kids)):
                if kids[a].L and kids[b].L:
                    A, B = np.array(kids[a].L), np.array(kids[b].L)
                    M[np.ix_(A, B)] = spec
                    M[np.ix_(B, A)] = spec
    return ix, M


def tag_accuracy(tagged, locus_ix, truth):
    """Pairwise orthology accuracy over cross-species gene pairs of a rooted, tagged tree."""
    leaves = list(tagged.traverse_leaves())
    gi = {id(l): i for i, l in enumerate(leaves)}
    n = len(leaves)
    sp = np.array([species_of(l.label) for l in leaves])
    loc = np.array([locus_ix["_".join(l.label.split("_")[:2])] for l in leaves])
    P = np.zeros((n, n), dtype=bool)
    for v in tagged.traverse_postorder():
        if v.is_leaf():
            v.L = [gi[id(v)]]
            continue
        kids = v.child_nodes()
        v.L = sum((c.L for c in kids), [])
        o = getattr(v, "tag", "S") == "S"
        for a in range(len(kids)):
            for b in range(a + 1, len(kids)):
                A, B = np.array(kids[a].L), np.array(kids[b].L)
                P[np.ix_(A, B)] = o
                P[np.ix_(B, A)] = o
    T = truth[np.ix_(loc, loc)]
    mask = sp[:, None] != sp[None, :]
    mask = np.triu(mask, 1)
    tp = int((P & T & mask).sum())
    fp = int((P & ~T & mask).sum())
    fn = int((~P & T & mask).sum())
    tot = int(mask.sum())
    return tot, tot - fp - fn, tp, fp, fn


def root_split(t):
    a, b = t.root.child_nodes()[:2]
    la = frozenset(l.label for l in a.traverse_leaves())
    lb = frozenset(l.label for l in b.traverse_leaves())
    return frozenset([la, lb])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("repdir")
    ap.add_argument("genefile")
    ap.add_argument("outdir")
    ap.add_argument("--s0root", default="mindl")
    ap.add_argument("--outgroup", default="0")
    ap.add_argument("--tags", action="store_true", help="also compute tagging accuracy")
    ap.add_argument("--iter", action="store_true", help="one extra DISCO-R round from the ASTRAL-DISCO-R tree")
    ap.add_argument("--maxgenes", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    W = a.outdir
    res = {"rep": a.repdir, "genes": a.genefile, "s0root": a.s0root, "time": {}, "rf": {}}
    strue = open(os.path.join(a.repdir, "s_tree.trees")).readline().strip()
    # true species tree with species labels only (strip nothing: SimPhy s_tree uses plain ids)
    genes = read_trees(os.path.join(a.repdir, a.genefile))
    if a.maxgenes:
        genes = genes[: a.maxgenes]
    # remove any information carried by the input rooting (SimPhy true trees are rooted at the truth,
    # and DISCO breaks score ties in favour of the current root): reroot each tree on a random edge
    rng = np.random.default_rng(12345)
    rr = []
    for g in genes:
        t = treeswift.read_tree_newick(g)
        nodes = [n for n in t.traverse_preorder() if not n.is_root()]
        if len(nodes) > 2:
            discor.disco.reroot_on_edge(t, nodes[int(rng.integers(len(nodes)))])
            t.suppress_unifurcations()
        rr.append(t.newick())
    genes = rr
    # relabelled multi-copy trees for ASTRAL-Pro
    relab = os.path.join(W, "genes_sp.tre")
    with open(relab, "w") as f:
        for g in genes:
            t = treeswift.read_tree_newick(g)
            for l in t.traverse_leaves():
                l.label = species_of(l.label)
            f.write(t.newick() + "\n")

    def timed(key, fn, *args):
        t0 = time.time()
        r = fn(*args)
        res["time"][key] = res["time"].get(key, 0.0) + time.time() - t0
        return r

    # 1. ASTRAL-Pro 2 (also S0)
    s0 = timed("astral-pro2", astral_pro, relab, os.path.join(W, "apro2.tre"))
    res["rf"]["ASTRAL-Pro2"] = rf(s0, strue)

    # 2. Baseline DISCO
    def disco_run(name, mode, S=None, tagmode="overlap"):
        t0 = time.time()
        trees = []
        for g in genes:
            t = treeswift.read_tree_newick(g)
            discor.root_and_tag(t, mode, S, tagmode)
            trees.append(t)
        dec = discor.decompose_trees(trees)
        path = os.path.join(W, f"{name}.tre")
        with open(path, "w") as f:
            f.write("\n".join(dec) + "\n")
        res["time"][f"{name}-decomp"] = time.time() - t0
        res.setdefault("ntrees", {})[name] = len(dec)
        out = {}
        out["ASTRID"] = timed(f"{name}-astrid", astrid, path, path + ".astrid")
        out["ASTRAL"] = timed(f"{name}-astral", astral, path, path + ".astral")
        res["rf"][f"ASTRID-{name}"] = rf(out["ASTRID"], strue)
        res["rf"][f"ASTRAL-{name}"] = rf(out["ASTRAL"], strue)
        return out

    base = disco_run("DISCO", "disco")

    # 3. root S0
    t0 = time.time()
    if a.s0root == "outgroup":
        s0r = root_outgroup(s0, a.outgroup)
    else:
        s0r = root_mindl(s0, genes)
    res["time"]["s0-root"] = time.time() - t0
    t = treeswift.read_tree_newick(s0r)
    res["s0_root_correct"] = root_split(t) == root_split(treeswift.read_tree_newick(strue))
    S0 = discor.SpeciesTree(s0r)
    Strue = discor.SpeciesTree(strue)

    r1 = disco_run("DISCOR", "dl", S0, "overlap")
    disco_run("DISCOR-lca", "dl", S0, "lca")
    disco_run("DISCOR-oracle", "dl", Strue, "overlap")
    if a.iter:
        t1 = r1["ASTRAL"]
        s1r = root_outgroup(t1, a.outgroup) if a.s0root == "outgroup" else root_mindl(t1, genes)
        disco_run("DISCOR-it2", "dl", discor.SpeciesTree(s1r), "overlap")

    # 4. tagging accuracy vs truth
    if a.tags:
        st = treeswift.read_tree_newick(strue)
        heights(st)
        sh_ = [n.h for n in st.traverse_postorder(leaves=False)]
        ltrees = read_trees(os.path.join(a.repdir, "l_trees.trees"))
        gtrue = read_trees(os.path.join(a.repdir, "g_true.trees"))
        acc = {}
        for name, mode, S, tm in [("DISCO", "disco", None, "overlap"), ("DISCOR", "dl", S0, "overlap"),
                                  ("DISCOR-lca", "dl", S0, "lca"), ("DISCOR-oracle", "dl", Strue, "overlap")]:
            agg = np.zeros(5, dtype=np.int64)
            rootok = rootn = 0
            for i, g in enumerate(genes):
                if i >= len(ltrees):
                    break
                lix, M = true_locus_tags(ltrees[i], sh_)
                t = treeswift.read_tree_newick(g)
                if t.num_nodes(internal=False) < 2:
                    continue
                discor.root_and_tag(t, mode, S, tm)
                agg += np.array(tag_accuracy(t, lix, M))
                sps = [species_of(l.label) for l in t.traverse_leaves()]
                if a.genefile == "g_true.trees" and len(set(sps)) < len(sps):  # multi-copy trees only
                    rootn += 1
                    rootok += root_split(t) == root_split(treeswift.read_tree_newick(gtrue[i]))
            tot, ok, tp, fp, fn = agg.tolist()
            acc[name] = {"pairs": tot, "acc": ok / tot, "orth_prec": tp / max(1, tp + fp),
                         "orth_rec": tp / max(1, tp + fn)}
            if rootn:
                acc[name]["root_acc"] = rootok / rootn
        res["tags"] = acc

    with open(os.path.join(W, "result.json"), "w") as f:
        json.dump(res, f)
    print(json.dumps(res))


if __name__ == "__main__":
    main()
