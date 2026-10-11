"""Large-k exploration: plain vs corrected internode distances under iid deletion (true gene trees)."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import msc, stree
from msc import Node

def delete(genes, p, rng, taxa):
    out = []
    for g in genes:
        keep = [t for t in taxa if rng.random() > p]
        if len(keep) < 4: continue
        t = stree.parse(g).extract_tree_with(set(keep)); t.suppress_unifurcations()
        out.append(t.newick())
    return out

def mixed_cat(n, short, long_):
    # caterpillar: alternating short / long internal branches
    node = Node([Node(label="t0"), Node(label="t1")], length=short)
    for i in range(2, n):
        node = Node([node, Node(label=f"t{i}")], length=short if i % 2 == 0 else long_)
    node.length = np.inf
    return node

def unbal(nA, nB, f, L):
    # a short branch f separating a big caterpillar clade from a small clade, long elsewhere
    A = msc.caterpillar(nA, L); A.length = f
    B = msc.balanced(nB, L) if nB > 1 else Node(label="b0")
    for lf in [x for x in []]: pass
    return A, B

if __name__ == "__main__":
    rng = np.random.default_rng(1)
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    cases = []
    for n in [8, 12, 16]:
        for f in [0.05, 0.1, 0.3]:
            cases.append((f"cat{n}-f{f}", msc.caterpillar(n, f)))
        cases.append((f"mixcat{n}-0.05/2", mixed_cat(n, 0.05, 2.0)))
    cases.append(("bal16-0.1", msc.balanced(16, 0.1)))
    for name, sp in cases:
        true = msc.species_newick(sp)
        n = true.count(",") + 1
        taxa = [f"t{i}" for i in range(n)]
        g = msc.sim_genes(sp, K, int(rng.integers(1e9)))
        res = []
        for p in [0.0, 0.5, 0.7]:
            gd = g if p == 0 else delete(g, p, rng, taxa)
            ph = stree.estimate_p(gd, n) if p > 0 else 0.0
            row = [f"p={p}"]
            for mode in ["plain", "ht", "plugin"]:
                if p == 0 and mode != "plain": continue
                M, C = stree.distance_matrix(gd, taxa, mode=mode, p=p)
                row.append(f"{mode}:{stree.fn_fp(true, stree.fastme(M, taxa))[0]:.2f}/{stree.fn_fp(true, stree.nj(M, taxa))[0]:.2f}")
            res.append(" ".join(row))
        print(name, " | ".join(res), flush=True)
