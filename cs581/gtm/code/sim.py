"""Simulation: where must an unblended merger (GTM) fail and can blending help?

Per replicate: model tree (caterpillar or Yule) -> AliSim GTR+G alignment ->
FastTree guide tree -> centroid decomposition of the guide (max subset size) ->
subset trees (true restricted = 'exact', or FastTree/IQ-TREE on the subset) ->
mergers: GTM (guide), constrained-SPR parsimony from GTM, parsimony constrained
insertion (+SPR polish). Error = FN vs model tree.
Usage: python3 sim.py SHAPE REP OUTDIR [n] [maxsub] [subtrees: exact|fasttree|iqtree]
"""
import json
import os
import random
import subprocess
import sys
import time
from phylo import Tree, read_tree, fn_fp, is_induced, parse_newick
from datasets import read_fasta
from blend import BlendSearch, compress
import insert

GTM = os.environ.get("GTM_SRC", "/opt/gtmdata/GTM_src")
IQ = os.environ.get("IQTREE", "/opt/mm/root/envs/bio/bin/iqtree3")


def model_tree(shape, n, rng, internal_mean, pendant_mean):
    """Newick with branch lengths. shape: cat | yule."""
    names = ["t%d" % i for i in range(n)]
    ex = lambda m: rng.expovariate(1.0 / m)
    if shape == "cat":
        s = "(%s:%.5f,%s:%.5f)" % (names[0], ex(pendant_mean), names[1], ex(pendant_mean))
        for i in range(2, n):
            s = "(%s:%.5f,%s:%.5f)" % (s, ex(internal_mean), names[i], ex(pendant_mean))
        return s + ";"
    # Yule: random joins of lineages
    nodes = [(nm, True) for nm in names]
    while len(nodes) > 1:
        i, j = rng.sample(range(len(nodes)), 2)
        a, b = nodes[i], nodes[j]
        la = ex(pendant_mean) if a[1] else ex(internal_mean)
        lb = ex(pendant_mean) if b[1] else ex(internal_mean)
        new = ("(%s:%.5f,%s:%.5f)" % (a[0], la, b[0], lb), False)
        nodes = [x for k, x in enumerate(nodes) if k not in (i, j)] + [new]
    return nodes[0][0] + ";"


def centroid_decomp(t, maxsize):
    """Recursive centroid-edge decomposition of phylo.Tree t -> list of leaf-name lists."""
    names = [t.label[v] for v in t.leaves()]
    if len(names) <= maxsize:
        return [names]
    # leaf counts below each node, rooted at first leaf
    root = t.leaves()[0]
    par = {root: None}
    order = []
    st = [root]
    while st:
        v = st.pop()
        order.append(v)
        for u in t.adj[v]:
            if u not in par:
                par[u] = v
                st.append(u)
    cnt = {}
    for v in reversed(order):
        cnt[v] = (1 if v in t.label else 0) + sum(cnt[u] for u in t.adj[v] if u != par[v])
    N = len(names)
    best = min((v for v in order if par[v] is not None), key=lambda v: abs(N - 2 * cnt[v]))
    side = set()
    st = [best]
    while st:
        v = st.pop()
        if v in t.label:
            side.add(t.label[v])
        for u in t.adj[v]:
            if u != par[v]:
                st.append(u)
    a = t.copy().restrict(side)
    b = t.copy().restrict(set(names) - side)
    return centroid_decomp(a, maxsize) + centroid_decomp(b, maxsize)


def write_fasta(seqs, names, path):
    with open(path, "w") as f:
        for n in names:
            f.write(">%s\n%s\n" % (n, seqs[n]))


def main():
    shape, rep, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    n = int(sys.argv[4]) if len(sys.argv) > 4 else 500
    imean = float(os.environ.get("INTERNAL_MEAN", "0"))
    maxsub = int(sys.argv[5]) if len(sys.argv) > 5 else 100
    stmode = sys.argv[6] if len(sys.argv) > 6 else "exact"
    imean = imean or (0.01 if shape == "cat" else 0.02)
    d = f"{out}/{shape}_n{n}_m{maxsub}_i{imean}/{rep}"
    os.makedirs(d, exist_ok=True)
    rng = random.Random(1000 * rep + (0 if shape == "cat" else 1))
    if not os.path.exists(f"{d}/true.tre"):
        nwk = model_tree(shape, n, rng, internal_mean=imean, pendant_mean=0.1)
        open(f"{d}/true.tre", "w").write(nwk + "\n")
        subprocess.run([IQ, "--alisim", f"{d}/aln", "-t", f"{d}/true.tre", "-m", "GTR{1.0,3.0,1.0,1.0,3.0}+F{0.25,0.25,0.25,0.25}+G4{0.5}",
                        "--length", "1000", "-af", "fasta", "-seed", str(rep + 1), "--no-unaligned", "-redo"],
                       check=True, capture_output=True)
        with open(f"{d}/guide.tre", "w") as f:
            subprocess.run(["FastTree", "-nt", "-gtr", "-gamma", "-quiet", "-nopr", f"{d}/aln.fa"],
                           check=True, stdout=f, stderr=subprocess.DEVNULL)
    T = read_tree(f"{d}/true.tre")
    GU = read_tree(f"{d}/guide.tre")
    seqs = read_fasta([f"{d}/aln.fa"])
    names = sorted(seqs)
    subsets = centroid_decomp(GU, maxsub)
    sd = f"{d}/sub_{stmode}"
    os.makedirs(sd, exist_ok=True)
    subfiles = []
    for i, s in enumerate(subsets):
        p = f"{sd}/s{i}.tre"
        if not os.path.exists(p):
            if stmode == "exact":
                open(p, "w").write(T.copy().restrict(s).to_newick() + "\n")
            else:
                write_fasta(seqs, s, f"{sd}/s{i}.fa")
                if stmode == "fasttree":
                    with open(p, "w") as f:
                        subprocess.run(["FastTree", "-nt", "-gtr", "-gamma", "-quiet", "-nopr", f"{sd}/s{i}.fa"],
                                       check=True, stdout=f, stderr=subprocess.DEVNULL)
                else:
                    subprocess.run([IQ, "-s", f"{sd}/s{i}.fa", "-m", "GTR+G", "-T", "1", "--prefix", f"{sd}/s{i}",
                                    "-seed", "1", "-quiet", "-redo"], check=True, capture_output=True)
                    os.replace(f"{sd}/s{i}.treefile", p)
        subfiles.append(p)
    subs = [read_tree(p) for p in subfiles]
    res = dict(shape=shape, rep=rep, n=n, maxsub=maxsub, subtrees=stmode, k=len(subsets),
               fn_guide=fn_fp(GU, T)[0],
               fn_subsets=sum(fn_fp(s, T)[2] for s in subs) / max(1, sum(fn_fp(s, T)[4] for s in subs)))
    # GTM
    subprocess.run([sys.executable, f"{GTM}/gtm.py", "-s", f"{d}/guide.tre", "-t", *subfiles, "-o", f"{sd}/gtm.tre"],
                   check=True, capture_output=True)
    G = read_tree(f"{sd}/gtm.tre")
    res["fn_gtm"] = fn_fp(G, T)[0]
    # GTM with true tree as guide = best unblended merger
    subprocess.run([sys.executable, f"{GTM}/gtm.py", "-s", f"{d}/true.tre", "-t", *subfiles, "-o", f"{sd}/gtm_trueguide.tre"],
                   check=True, capture_output=True)
    res["fn_best_unblended"] = fn_fp(read_tree(f"{sd}/gtm_trueguide.tre"), T)[0]
    states, w = compress(seqs, names)
    if os.environ.get("QUICK"):
        o = insert.run(names, states, w, subs, "oracle", ref_tree=T)
        res["fn_oracle_ins"] = fn_fp(parse_newick(o.to_newick()), T)[0]
        json.dump(res, open(f"{sd}/quick.json", "w"), indent=1)
        print(json.dumps(res), flush=True)
        return
    subnames = [[t.label[v] for v in t.leaves()] for t in subs]
    # constrained parsimony SPR from GTM
    t0 = time.time()
    bs = BlendSearch(G, names, states, w, subnames, radius=6)
    res["pars_gtm"] = bs.length
    res["spr_moves"] = bs.run_parsimony()
    B = parse_newick(bs.to_newick())
    res.update(fn_gtm_spr=fn_fp(B, T)[0], pars_gtm_spr=bs.length, sec_gtm_spr=time.time() - t0)
    open(f"{sd}/gtm_spr.tre", "w").write(bs.to_newick() + "\n")
    # parsimony constrained insertion (+ polish)
    t0 = time.time()
    g = insert.run(names, states, w, subs, "pars", ref_tree=GU)
    I = parse_newick(g.to_newick())
    bi = BlendSearch(I, names, states, w, subnames, radius=6)
    res.update(fn_ins=fn_fp(I, T)[0], pars_ins=bi.length)
    bi.run_parsimony()
    P = parse_newick(bi.to_newick())
    res.update(fn_ins_spr=fn_fp(P, T)[0], pars_ins_spr=bi.length, sec_ins=time.time() - t0)
    open(f"{sd}/ins_spr.tre", "w").write(bi.to_newick() + "\n")
    # oracle insertion (headroom)
    o = insert.run(names, states, w, subs, "oracle", ref_tree=T)
    res["fn_oracle_ins"] = fn_fp(parse_newick(o.to_newick()), T)[0]
    # model-tree parsimony (sanity: is MP a good criterion here?)
    bt = BlendSearch(T, names, states, w, [], radius=1)
    res["pars_true"] = bt.length
    res["constraints_ok"] = all(is_induced(X, s) for X in (B, P) for s in subs)
    json.dump(res, open(f"{sd}/result.json", "w"), indent=1)
    print(json.dumps(res), flush=True)


if __name__ == "__main__":
    main()
