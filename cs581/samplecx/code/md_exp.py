"""Missing-data experiment: corrected ASTRID variants vs ASTRID, Asteroid, ASTRAL.

Data: Nute et al. 2018 (doi:10.13012/B2IDB-7735354_V1), 25 ingroup + outgroup '0'.
  * 'rand-30', 'rand-60': the published i.i.d. missing-data replicates (true + RAxML gene trees,
    first 50/200/1000 genes, as in the paper).
  * 'sim-50', 'sim-70': our own i.i.d. deletion (p = 0.5, 0.7) applied to the complete
    ('full') replicates, 1000 genes; seed = replicate.
Methods: astrid (unrooted NJst distance + FastME), astrid-ht (exact unbiased weights, global p_hat),
  astrid-plugin (bounded plug-in weights), astrid-cmp (gene trees completed against the first-pass
  ASTRID tree, then ASTRID), asteroid, astral (ASTER), plus the published astrid/astral trees.
Output: /opt/runs/md_exp.jsonl (restartable).
"""
import sys, os, json, time
from multiprocessing import Pool
import numpy as np
import treeswift as ts
sys.path.insert(0, os.path.dirname(__file__))
import stree

ROOT = "/opt/data/nute"
OUT = "/opt/runs/md_exp.jsonl"
HEIGHTS = ["10M", "2M", "500K"]
RATES = ["1E-7", "1E-6"]


def load(h, r, pat, ng, rep, gt):
    if pat.startswith("rand"):
        d = f"{ROOT}/25tax-1000gen-0bps-{h}-{r}-{pat}/{rep:02d}"
        f = "true-genes-w-missing-data-na.tre" if gt == "true" else "raxml-genes.tre"
        genes = stree.read_trees(f"{d}/{f}", limit=ng)
    else:
        p = int(pat.split("-")[1]) / 100
        d = f"{ROOT}/25tax-1000gen-0bps-{h}-{r}-full/{rep:02d}"
        f = "true-genes.tre" if gt == "true" else "raxml-genes.tre"
        genes0 = stree.read_trees(f"{d}/{f}", limit=ng)
        rng = np.random.default_rng(1000 * rep + int(p * 100))
        genes = []
        for g in genes0:
            t = stree.parse(g)
            labs = list(t.labels(internal=False))
            keep = {x for x in labs if rng.random() > p}
            if len(keep) < 4:
                continue
            t = t.extract_tree_with(keep)
            t.suppress_unifurcations()
            genes.append(t.newick())
    genes = [stree.strip(g) for g in genes if g.count(",") >= 3]  # >= 4 leaves
    true = open(f"{ROOT}/25tax-1000gen-0bps-{h}-{r}-full/{rep:02d}/true-species.tre").read().strip()
    return genes, true


def run(job):
    h, r, pat, ng, rep, gt = job
    genes, true = load(h, r, pat, ng, rep, gt)
    taxa = sorted(stree.bipartitions(true)[1])
    n = len(taxa)
    ph = stree.estimate_p(genes, n)
    ests, secs = {}, {}

    def timed(name, fn):
        t0 = time.time()
        try:
            ests[name] = fn()
        except Exception as e:
            ests[name] = None
        secs[name] = time.time() - t0

    def astrid_mode(mode):
        M, C = stree.distance_matrix(genes, taxa, mode=mode, p=ph)
        return stree.fastme(M, taxa)
    timed("astrid", lambda: astrid_mode("plain"))
    timed("astrid-ht", lambda: astrid_mode("ht"))
    timed("astrid-plugin", lambda: astrid_mode("plugin"))

    def cmp():
        t0 = time.time()
        ref = ts.read_tree_newick(ests["astrid"])
        root_leaf = "0"
        ref.reroot(next(l for l in ref.traverse_leaves() if l.label == root_leaf)) if False else None
        # root the reference at leaf '0' by rebuilding the newick around it
        ref = ts.read_tree_newick(reroot_at_leaf(ests["astrid"], root_leaf))
        comp = [stree.complete_gene(g, ref, taxa, root_leaf) for g in genes]
        M, C = stree.distance_matrix(comp, taxa)
        return stree.fastme(M, taxa)
    timed("astrid-cmp", cmp)
    secs["astrid-cmp"] += secs["astrid"]
    timed("asteroid", lambda: stree.asteroid(genes))
    timed("astral", lambda: stree.astral(genes))
    rows = []
    for m, e in ests.items():
        if e is None:
            continue
        fn, fp, _ = stree.fn_fp(true, e)
        rows.append(dict(height=h, rate=r, pattern=pat, ngenes=ng, rep=rep, genetrees=gt,
                         method=m, FN=fn, FP=fp, sec=secs[m], p_hat=ph, ngenes_used=len(genes)))
    if pat.startswith("rand"):
        d = f"{ROOT}/25tax-{ng}gen-0bps-{h}-{r}-{pat}/{rep:02d}"
        for m in ["astrid", "astral"]:
            f = f"{d}/{m}-{gt}-genes.tre"
            if os.path.exists(f) and os.path.getsize(f) > 0:
                e = open(f).read().strip().split("\n")[-1]
                fn, fp, _ = stree.fn_fp(true, e)
                rows.append(dict(height=h, rate=r, pattern=pat, ngenes=ng, rep=rep, genetrees=gt,
                                 method=m + "-pub", FN=fn, FP=fp, sec=None, p_hat=ph,
                                 ngenes_used=len(genes)))
    return rows


def reroot_at_leaf(nwk, leaf):
    labels, adj, leaves = stree.to_adj(stree.parse(nwk))
    lab = {lf: l for lf, l in zip(leaves, labels)}
    a = next(lf for lf, l in lab.items() if l == leaf)
    nb = adj[a][0]
    # suppress degree-2 nodes first by walking; simple recursive writer from nb excluding a
    def rec(u, p):
        ch = [w for w in adj[u] if w != p]
        if not ch:
            return lab[u]
        if len(ch) == 1:
            return rec(ch[0], u)
        return "(" + ",".join(rec(w, u) for w in ch) + ")"
    return f"({rec(nb, a)},{leaf});"


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    done = set()
    if os.path.exists(OUT):
        for line in open(OUT):
            d = json.loads(line)
            done.add((d["height"], d["rate"], d["pattern"], d["ngenes"], d["rep"], d["genetrees"]))
    jobs = []
    for rep in range(1, 21):
        for h in HEIGHTS:
            for r in RATES:
                for gt in ["true", "raxml"]:
                    if which in ("all", "rand"):
                        for pat in ["rand-30", "rand-60"]:
                            for ng in [50, 200, 1000]:
                                jobs.append((h, r, pat, ng, rep, gt))
                    if which in ("all", "sim"):
                        for pat in ["sim-50", "sim-70"]:
                            jobs.append((h, r, pat, 1000, rep, gt))
    jobs = [j for j in jobs if j not in done]
    print(len(jobs), "jobs", flush=True)
    with Pool(int(os.environ.get("NPROC", 4))) as pool, open(OUT, "a") as fh:
        for rows in pool.imap_unordered(run, jobs):
            for x in rows:
                fh.write(json.dumps(x) + "\n")
            fh.flush()
