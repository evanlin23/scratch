"""Simulation + methods for the Forest+DTM pilot.

Methods (all share the same alignment; distance methods share the same capped distance matrix):
  NJ        FastME -m N
  BIONJ     FastME -m I
  FastME    FastME -m B -n B -s   (balanced ME + BalME NNI + SPR)
  FastTree  FastTree -nt (JC+CAT; -gtr for K2P data)
  Forest    DMR forest, hyperparameters picked without the true tree (max #splits over a grid
            built from the data, following the guidance in Kim et al. 2026); scored as a partial tree
  F+GTM[guide,refine]
            forest components -> polytomies refined (refine=sub: by the guide method re-run on the
            component's own sub-matrix, DCM-style; refine=ind: by the full guide tree's induced
            splits) -> GTM (default 'convex' mode) with the guide tree (NJ or FastME)
"""
import os, sys, subprocess, tempfile, time, random
import numpy as np
import dendropy
from scipy.sparse.csgraph import minimum_spanning_tree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import forest_fast as F  # noqa: E402
from common import (read_fasta, jc_distance, k2p_distance, cap_distances, tree_splits,  # noqa: E402
                    splits_from_newick, score_tree, score_forest)

FASTME = os.environ.get("FASTME", os.path.expanduser("~/mamba/envs/phy/bin/fastme"))
FASTTREE = os.environ.get("FASTTREE", "fasttree")
IQTREE = os.environ.get("IQTREE", os.path.expanduser("~/mamba/envs/phy/bin/iqtree3"))
GTM = os.environ.get("GTM", os.path.expanduser("~/ext/GTM/gtm.py"))


# ---------------------------------------------------------------- simulation
def sim_tree(n, regime, rng):
    """Random birth-death topology (as in Kim et al. 2026), branch lengths by regime."""
    pyrng = random.Random(int(rng.integers(1 << 30)))
    t = dendropy.model.birthdeath.birth_death_tree(birth_rate=1.0, death_rate=0.5, num_extant_tips=n,
                                                   rng=pyrng, gsa_ntax=None)
    t.prune_extinct_lineages() if hasattr(t, "prune_extinct_lineages") else None
    for i, lf in enumerate(t.leaf_node_iter()):
        lf.taxon.label = f"t{i}"
    t.is_rooted = False
    t.collapse_basal_bifurcation()
    edges = [e for e in t.postorder_edge_iter() if e.tail_node is not None]
    kind, *par = regime.split(":")
    if kind == "U":  # uniform branch lengths U[a,b]
        a, b = float(par[0]), float(par[1])
        for e in edges:
            e.length = rng.uniform(a, b)
    elif kind == "UH":  # ultrametric BD tree scaled to root height h, lognormal(0, s) rate multipliers
        h, s = float(par[0]), float(par[1])
        # rescale to height h
        t.calc_node_ages(is_force_max_age=True)
        root_age = max(nd.age for nd in t.preorder_node_iter())
        for e in edges:
            e.length = e.length / root_age * h * rng.lognormal(0.0, s)
    else:
        raise ValueError(regime)
    return t


def simulate(n, k, regime, model, seed, workdir):
    rng = np.random.default_rng(seed)
    t = sim_tree(n, regime, rng)
    tf = os.path.join(workdir, "true.tre")
    with open(tf, "w") as fh:
        fh.write(t.as_string(schema="newick", suppress_rooting=True).strip() + "\n")
    prefix = os.path.join(workdir, "aln")
    m = {"JC": "JC", "K2P": "K2P{4}"}[model]
    subprocess.run([IQTREE, "--alisim", prefix, "-t", tf, "-m", m, "--length", str(k), "--seed", str(seed),
                    "--out-format", "fasta", "-redo", "-quiet"], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return tf, prefix + ".fa"


# ---------------------------------------------------------------- distance methods
def write_phylip(D, names, path):
    with open(path, "w") as fh:
        fh.write(f"{len(names)}\n")
        for i, nm in enumerate(names):
            fh.write(nm + " " + " ".join(f"{max(0.0, x) + 0.0:.8f}" for x in D[i]) + "\n")


def fastme(D, names, method, workdir, tag):
    """method: 'N' (NJ), 'I' (BIONJ), 'B' (BME+NNI+SPR)."""
    pin = os.path.join(workdir, f"{tag}.phy")
    pout = os.path.join(workdir, f"{tag}.tre")
    write_phylip(D, names, pin)
    args = [FASTME, "-i", pin, "-o", pout, "-m", method, "-T", "1"]
    if method == "B":
        args += ["-n", "B", "-s"]
    subprocess.run(args, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   cwd=workdir)
    with open(pout) as fh:
        s = fh.read().strip()
    for f in (pin, pout, pin + "_fastme_stat.txt"):
        if os.path.exists(f):
            os.remove(f)
    return s


def fasttree(aln, workdir, model):
    out = os.path.join(workdir, "ft.tre")
    args = [FASTTREE, "-nt", "-quiet", "-nopr"]
    if model == "K2P":
        args.append("-gtr")
    with open(aln) as fi, open(out, "w") as fo:
        subprocess.run(args, stdin=fi, stdout=fo, stderr=subprocess.DEVNULL, check=True)
    return open(out).read().strip()


# ---------------------------------------------------------------- forest + hyperparameters
def nj_internal_lengths(nwk, names):
    tns = dendropy.TaxonNamespace(names)
    t = dendropy.Tree.get(data=nwk, schema="newick", taxon_namespace=tns, rooting="force-unrooted")
    return np.array([e.length for e in t.postorder_edge_iter()
                     if e.length is not None and not e.head_node.is_leaf() and e.tail_node is not None])


def mst_bottleneck(D):
    Df = np.where(np.isfinite(D), D, 0.0)  # 0 = no edge for scipy
    Df = Df + 1e-12 * (np.isfinite(D) & (D == 0))
    T = minimum_spanning_tree(Df)
    return T.max() if T.nnz else 0.0


def forest_grid(D, nj_nwk, names, small=False):
    il = nj_internal_lengths(nj_nwk, names)
    il = il[il > 1e-9]
    if len(il) == 0:
        il = np.array([0.01])
    qs = [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7] if not small else [0.05, 0.2, 0.5]
    taus = sorted(set(round(float(x), 5) for x in 0.5 * np.quantile(il, qs) if x > 1e-5))
    b = mst_bottleneck(D)
    ms = [b * f for f in ((1.0001, 0.8, 0.6, 0.45, 0.3) if not small else (1.0001, 0.6, 0.3))]
    grid = []
    for m in ms:
        tt = [t for t in taus if m > 3 * t]
        if not tt or m <= 0:
            continue
        M1 = 2 * m + 3 * max(tt) + 1e-3
        for M in ([M1, 1.5 * M1] if not small else [M1]):
            grid.append((m, M, tt))
    return grid


def best_forest(D, grid):
    """Max total #splits among valid forests (ties: larger tau, then fewer components)."""
    best = None
    tried = valid = 0
    for (m, M, taus) in grid:
        res = F.run_forest(D, m, M, taus)
        for tau in taus:
            r = res[tau]
            tried += 1
            if not r["valid"]:
                continue
            valid += 1
            ns = sum(len(s) for s in r["splits"])
            key = (ns, tau, -len(r["comps"]))
            if best is None or key > best[0]:
                best = (key, dict(m=m, M=M, tau=tau, comps=r["comps"], splits=r["splits"]))
    return (best[1] if best else None), tried, valid


# ---------------------------------------------------------------- merge
def induced_compatible(base, cand, full):
    """add every split of cand compatible with all of base (and with already-added ones)."""
    out = set(base)
    for s in sorted(cand, key=lambda x: -min(F.popcount(x), F.popcount(full & ~x))):
        if all(F.compatible(s, t, full) for t in out):
            out.add(s)
    return out


def restrict_global(splits, comp):
    sub = F.mask_of(comp)
    low = 1 << int(min(comp))
    c = len(comp)
    out = set()
    for s in splits:
        a = s & sub
        k = F.popcount(a)
        if 2 <= k <= c - 2:
            if a & low:
                a = sub & ~a
            out.add(a)
    return out


def refine_component(comp, fsplits, guide_splits, sub_splits_fn):
    """forest component splits + (subset-tree or induced guide) splits that are compatible."""
    full = F.mask_of(comp)
    cand = sub_splits_fn(comp) if sub_splits_fn else restrict_global(guide_splits, comp)
    s = induced_compatible(fsplits, cand, full)
    if len(s) < len(comp) - 3 and sub_splits_fn:  # top up with induced guide splits
        s = induced_compatible(s, restrict_global(guide_splits, comp), full)
    return s


def forest_gtm(forest, guide_nwk, names, D, guide_method, refine, workdir):
    n = len(names)
    G = splits_from_newick(guide_nwk, names)
    sub_fn = None
    if refine == "sub":
        def sub_fn(comp):
            comp = [int(i) for i in comp]
            if len(comp) < 4:
                return set()
            sn = [names[i] for i in comp]
            nw = fastme(D[np.ix_(comp, comp)], sn, guide_method, workdir, "sub")
            loc = splits_from_newick(nw, sn)  # bitmasks over local indices
            out = set()
            for s in loc:
                g = 0
                for li, gi in enumerate(comp):
                    if (s >> li) & 1:
                        g |= 1 << gi
                out.add(g)
            return restrict_global(out, comp)
    cdir = tempfile.mkdtemp(dir=workdir)
    files = []
    unresolved = 0
    for ci, comp in enumerate(forest["comps"]):
        S = refine_component(comp, forest["splits"][ci], G, sub_fn) if len(comp) >= 4 else set()
        unresolved += max(0, len(comp) - 3 - len(S))
        nwk = F.splits_to_newick(S, comp, names)
        p = os.path.join(cdir, f"c{ci}.tre")
        with open(p, "w") as fh:
            fh.write(nwk + "\n")
        files.append(p)
    gp = os.path.join(cdir, "guide.tre")
    with open(gp, "w") as fh:
        fh.write(guide_nwk.strip() + "\n")
    op = os.path.join(cdir, "out.tre")
    subprocess.run([sys.executable, GTM, "-s", gp, "-t", *files, "-o", op], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = open(op).read().strip()
    for f in files + [gp, op]:
        os.remove(f)
    os.rmdir(cdir)
    return out, unresolved


# ---------------------------------------------------------------- one replicate
def run_replicate(n, k, regime, model, seed, workdir, small_grid=False, do_fasttree=True):
    rec = dict(n=n, k=k, regime=regime, model=model, seed=seed)
    tf, aln = simulate(n, k, regime, model, seed, workdir)
    names, X = read_fasta(aln)
    tns = dendropy.TaxonNamespace(names)
    tt = dendropy.Tree.get(path=tf, schema="newick", taxon_namespace=tns, rooting="force-unrooted")
    T = tree_splits(tt, names)
    t0 = time.time()
    D0 = jc_distance(X) if model == "JC" else k2p_distance(X)
    rec["frac_saturated"] = float((~np.isfinite(D0)).sum() / (n * (n - 1)))
    D = cap_distances(D0)
    rec["t_dist"] = time.time() - t0
    trees = {}
    for meth, code in (("NJ", "N"), ("BIONJ", "I"), ("FastME", "B")):
        t0 = time.time()
        trees[meth] = fastme(D, names, code, workdir, meth)
        rec[f"t_{meth}"] = time.time() - t0
    if do_fasttree:
        t0 = time.time()
        trees["FastTree"] = fasttree(aln, workdir, model)
        rec["t_FastTree"] = time.time() - t0
    for meth, nw in trees.items():
        fn, fp, _ = score_tree(splits_from_newick(nw, names), T, n)
        rec[f"FN_{meth}"], rec[f"FP_{meth}"] = fn, fp
    # forest (D0: unreliable distances stay infinite -> never below m or M)
    t0 = time.time()
    grid = forest_grid(D0, trees["NJ"], names, small=small_grid)
    fo, tried, valid = best_forest(D0, grid)
    rec["t_Forest"] = time.time() - t0
    rec["forest_grid_tried"], rec["forest_grid_valid"] = tried, valid
    if fo is None:  # no valid forest: every leaf a singleton (Forest+GTM == guide tree)
        comps = [np.array([i]) for i in range(n)]
        fo = dict(m=None, M=None, tau=None, comps=comps, splits=[set() for _ in comps])
    fn, fp, nfalse, ntot = score_forest(fo["comps"], fo["splits"], T, n)
    rec.update(FN_Forest=fn, FP_Forest=fp, forest_false=nfalse, forest_splits=ntot,
               forest_ncomp=len(fo["comps"]), forest_maxcomp=int(max(len(c) for c in fo["comps"])),
               forest_m=fo["m"], forest_M=fo["M"], forest_tau=fo["tau"])
    for gm, gcode in (("NJ", "N"), ("FastME", "B")):
        for refine in ("sub", "ind"):
            t0 = time.time()
            nw, unres = forest_gtm(fo, trees[gm], names, D, gcode, refine, workdir)
            name = f"FGTM_{gm}_{refine}"
            rec[f"t_{name}"] = time.time() - t0 + rec["t_Forest"] + rec[f"t_{gm}"]
            fn, fp, _ = score_tree(splits_from_newick(nw, names), T, n)
            rec[f"FN_{name}"], rec[f"FP_{name}"] = fn, fp
            rec[f"unres_{name}"] = unres
    return rec
