"""Supertree methods: baselines (ASTRAL-III, ASTER astral4, TREE-QMC) and the divide-and-conquer
pipeline DC-GTM (guide tree -> centroid-edge decomposition -> subset supertrees -> GTM merge)."""
import os, sys, time
import treeswift
from common import ENV, ASTRAL3, GTM, read_trees, leaves, restrict, strip_newick, timed


def _clean_input(src, dst):
    """Topology-only copy of source trees (ASTRAL-III chokes on some branch-length formats)."""
    with open(dst, "w") as f:
        for t in read_trees(src):
            f.write(strip_newick(t) + "\n")


def astral3(src, out, wd, mem="6g"):
    os.makedirs(wd, exist_ok=True)
    inp = os.path.join(wd, "in.tre")
    _clean_input(src, inp)
    return timed(f"java -Xmx{mem} -jar {ASTRAL3} -i {inp} -o {out}", log=os.path.join(wd, "astral3.log"))


def astral4(src, out, wd, threads=1, extra=""):
    os.makedirs(wd, exist_ok=True)
    return timed(f"{ENV}/astral4 -t {threads} {extra} -u 0 -o {out} {src}", log=os.path.join(wd, "astral4.log"))


def tqmc(src, out, wd, extra=""):
    os.makedirs(wd, exist_ok=True)
    return timed(f"{ENV}/tree-qmc {extra} -i {src} -o {out}", log=os.path.join(wd, "tqmc.log"))


# ---------------------------------------------------------------- decomposition

def centroid_decompose(tree, max_size):
    """Recursively split an unrooted tree at its centroid edge until every leaf set has
    <= max_size leaves. Returns a list of leaf-label sets (disjoint, covering)."""
    all_leaves = set(leaves(tree))
    out = []
    stack = [(tree, all_leaves)]
    while stack:
        t, ls = stack.pop()
        if len(ls) <= max_size:
            out.append(ls)
            continue
        n = len(ls)
        size = {}
        for nd in t.traverse_postorder():
            size[nd] = 1 if nd.is_leaf() else sum(size[c] for c in nd.children)
        best, best_val = None, None
        for nd in t.traverse_preorder():
            if nd is t.root:
                continue
            val = max(size[nd], n - size[nd])
            if best_val is None or val < best_val:
                best, best_val = nd, val
        side = set(l.label for l in best.traverse_leaves())
        other = ls - side
        for part in (side, other):
            sub = restrict(t, part) if len(part) >= 3 else None
            if sub is None or len(part) <= max_size:
                out.append(part)
            else:
                stack.append((sub, part))
    return out


# ---------------------------------------------------------------- DC-GTM

def dc_gtm(src, guide_path, out, wd, max_size=100, subset_method="astral4", threads=1, workers=4,
           guide_for_merge=None):
    """Divide-and-conquer supertree. Returns dict with timings."""
    from concurrent.futures import ThreadPoolExecutor
    os.makedirs(wd, exist_ok=True)
    t0 = time.time()
    guide = treeswift.read_tree_newick(open(guide_path).read().strip())
    parts = centroid_decompose(guide, max_size)
    srcs = read_trees(src)
    t_decomp = time.time() - t0

    guide_added = []

    def solve(i):
        part = parts[i]
        sd = os.path.join(wd, f"sub{i}")
        os.makedirs(sd, exist_ok=True)
        sp = os.path.join(sd, "src.tre")
        nkept = 0
        covered = set()
        with open(sp, "w") as f:
            for t in srcs:
                r = restrict(t, part)
                if r is not None and len(leaves(r)) >= 4:
                    f.write(strip_newick(r) + "\n")
                    covered.update(leaves(r))
                    nkept += 1
            # taxa that no restricted source tree (>=4 leaves) covers would be dropped by
            # ASTRAL; add the guide tree restricted to the subset as one extra input tree
            if nkept and covered != part and len(part) >= 4:
                f.write(strip_newick(restrict(guide, part)) + "\n")
                guide_added.append(i)
        so = os.path.join(sd, "out.tre")
        if len(part) <= 3 or nkept == 0:
            # trivial subset: star / restricted guide tree
            g = restrict(guide, part)
            with open(so, "w") as f:
                f.write((strip_newick(g) if g is not None else "(" + ",".join(sorted(part)) + ");") + "\n")
            return 0.0, 0.0
        if subset_method == "astral4":
            w, m = astral4(sp, so, sd, threads=threads)
        elif subset_method == "astral4R":
            w, m = astral4(sp, so, sd, threads=threads, extra="-R")
        elif subset_method == "astral3":
            w, m = astral3(sp, so, sd, mem="2g")
        elif subset_method == "tqmc":
            w, m = tqmc(sp, so, sd)
        else:
            raise ValueError(subset_method)
        return w, m

    t1 = time.time()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        res = list(ex.map(solve, range(len(parts))))
    t_sub = time.time() - t1
    sub_cpu = sum(r[0] for r in res)
    sub_peak = max((r[1] or 0) for r in res)

    # GTM needs subset trees with >=? leaves; single/two-leaf subsets are written as-is
    subfiles = []
    for i in range(len(parts)):
        p = os.path.join(wd, f"sub{i}", "out.tre")
        s = open(p).read().strip()
        if len(parts[i]) == 1:
            s = "(" + list(parts[i])[0] + ");"
            open(p, "w").write(s + "\n")
        subfiles.append(p)
    # topology-only guide for GTM
    gclean = os.path.join(wd, "guide_clean.tre")
    g2 = treeswift.read_tree_newick(open(guide_for_merge or guide_path).read().strip())
    open(gclean, "w").write(strip_newick(g2) + "\n")
    w_gtm, m_gtm = timed(f"python3 {GTM} -s {gclean} -t {' '.join(subfiles)} -o {out}",
                         log=os.path.join(wd, "gtm.log"))
    return {"n_subsets": len(parts), "max_subset": max(len(p) for p in parts),
            "t_decomp": t_decomp, "t_sub_wall": t_sub, "t_sub_cpu": sub_cpu, "t_gtm": w_gtm,
            "peak_sub_mb": sub_peak, "peak_gtm_mb": m_gtm,
            "n_guide_added": len(guide_added), "wall": time.time() - t0}


# ---------------------------------------------------------------- MRL via FastTree (fast guide tree)

def mrp_matrix(src, out_fasta):
    """MRP matrix as a pseudo-DNA alignment (A=0, C=1, -=missing), one column per non-trivial
    source-tree bipartition. Returns (ntaxa, ncols)."""
    trees = read_trees(src)
    taxa = sorted(set(l for t in trees for l in leaves(t)))
    ti = {l: i for i, l in enumerate(taxa)}
    cols = []
    for t in trees:
        tl = [ti[l] for l in leaves(t)]
        below = {}
        for nd in t.traverse_postorder():
            if nd.is_leaf():
                below[nd] = [ti[nd.label]]
            else:
                b = []
                for c in nd.children:
                    b.extend(below[c])
                below[nd] = b
                if nd is not t.root and 2 <= len(b) <= len(tl) - 2:
                    cols.append((tl, b))
                for c in nd.children:  # free memory of children lists
                    if not c.is_leaf():
                        below[c] = None
    n = len(taxa)
    rows = [bytearray(b"-" * len(cols)) for _ in range(n)]
    for j, (tl, b) in enumerate(cols):
        for i in tl:
            rows[i][j] = 65  # A
        for i in b:
            rows[i][j] = 67  # C
    with open(out_fasta, "w") as f:
        for i, l in enumerate(taxa):
            f.write(f">{l}\n{rows[i].decode()}\n")
    return n, len(cols)


def mrl_fasttree(src, out, wd, threads=1):
    os.makedirs(wd, exist_ok=True)
    fa = os.path.join(wd, "mrp.fa")
    t0 = time.time()
    mrp_matrix(src, fa)
    tm = time.time() - t0
    exe = f"OMP_NUM_THREADS={threads} {ENV}/FastTreeMP" if threads > 1 else f"{ENV}/FastTree"
    w, m = timed(f"{exe} -nt -nosupport -quiet {fa} > {out}", log=os.path.join(wd, "ft.log"))
    return w + tm, m
