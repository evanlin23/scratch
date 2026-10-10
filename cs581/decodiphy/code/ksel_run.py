"""Paired k-selection pilot on the DecoDiPhy E1 setup.

For each (biological tree, true k, seed, noise level): simulate queries + (noisy) mixture distances
with the authors' own generator (MISC/simulate_exp.py), then run the authors' greedy+refinement
search (decodiphy.search.hill_climbing with the OSQP small-problem) for EVERY k = 1..k_true+3,
recording loss/anchors/p/x/ybar per k. Stopping rules are then compared offline on identical
trajectories (paired). Restartable: one JSON per run.

usage: python ksel_run.py DECODIPHY_REPO TREES_DIR OUTDIR [nproc]
"""
import os, sys, json, time, random, io, contextlib, itertools
os.environ["OMP_NUM_THREADS"] = "1"; os.environ["OPENBLAS_NUM_THREADS"] = "1"
import numpy as np
from multiprocessing import Pool

REPO, TREES, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
NPROC = int(sys.argv[4]) if len(sys.argv) > 4 else 4
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, "MISC"))

TREE_FILES = {
    "bees": "bees-bosseret/bees-bosseret.tre",
    "mammals-song": "mammals-song/mammals-song-castlespro.tre",
    "birds-jarvis": "birds-jarvis/birds-jarvis.tre.rerooted",
    "tilapia": "tilapia-Ciezarek/tilapia-Ciezarek.tre.rerooted",
    "1kp": "1kp/1kp-concat-fig2.tre",
    "pancrustacean": "pancrustacean-Bernot/pancrustacean-Bernot.tre.rerooted",
    "beetles": "beetles-Johnson/beetles-Johnson.tre",
    "fish-troyer": "fish-Troyer/fish-Troyer.tre.rerooted",
    "hemipteroid": "hemipteroid-johnson/hemipteroid-johnson.tre.rerooted",
}
KS = [2, 3, 5, 7, 10]
SEEDS = range(1, 11)
NOISE = {"noise0": (0, 1.0), "noise1": (1, 1.0), "noise2": (1, 2.0)}  # (add_noise, exp_scale)
READS = 10 ** 5


def one(args):
    try:
        return _one(args)
    except Exception as e:  # e.g. the authors' OSQP failure path (UnboundLocalError in optimize.py)
        tname, k, seed, nname = args
        out = os.path.join(OUT, f"{tname}_k{k}_s{seed}_{nname}.failed")
        open(out, "w").write(repr(e))
        return out


def _one(args):
    tname, k, seed, nname = args
    out = os.path.join(OUT, f"{tname}_k{k}_s{seed}_{nname}.json")
    if os.path.exists(out) or os.path.exists(out[:-5] + ".failed"):
        return out
    import simulate_exp as se
    from treeswift import read_tree_newick
    from decodiphy.tree_utils import label_tree, preprocess_input
    from decodiphy.matrix import get_input_matrices
    from decodiphy.neighbors import get_neighbors
    from decodiphy.search import hill_climbing
    from decodiphy.optimize import AnchorOSQP

    work = out[:-5] + "_sim"
    os.makedirs(work, exist_ok=True)
    random.seed(seed); np.random.seed(seed)
    with open(os.path.join(TREES, TREE_FILES[tname])) as f:
        nwk = f.read().strip().split("\n")[0]
    add_noise, es = NOISE[nname]
    with contextlib.redirect_stdout(io.StringIO()):
        se.simulate_queries(read_tree_newick(nwk), None, k, work, add_noise=add_noise, read_count=READS, exp_scale=es)
    true = [l.split() for l in open(os.path.join(work, "true_queries.txt"))]
    distances = {a: float(b) for a, b in (l.split() for l in open(os.path.join(work, "distances.txt")))}
    tree_obj = read_tree_newick(open(os.path.join(work, "pruned_tree.trees")).read().strip())
    label_tree(tree_obj)
    preprocess_input(tree_obj, distances)
    d, l, C, D, index_to_node, node_to_index, index_to_leaf, bl, tl = get_input_matrices(tree_obj, distances)
    D_T, C_T, l_T = D.T, C.T, l.reshape(-1)
    nbrs = get_neighbors(tree_obj, node_to_index, 2)
    # node adjacency (for structural features of each solution)
    lab2node = tree_obj.label_to_node(selection="all")
    parent = {lab: (nd.parent.label if nd.parent is not None else None) for lab, nd in lab2node.items()}
    np.random.seed(1142)
    rounds, anchors = [], None
    K = k + 3
    for kk in range(1, min(K, len(index_to_node)) + 1):
        if anchors is None:
            anchors = np.random.choice(len(index_to_node), kk, replace=False)
        else:
            new = np.random.choice([i for i in range(len(index_to_node)) if i not in anchors], 1, replace=False)
            anchors = np.array(list(anchors) + list(new))
        t0 = time.time()
        solver = AnchorOSQP(d, l_T, C_T, D_T, anchors)
        with contextlib.redirect_stdout(io.StringIO()):
            anchors, obj, p, x, y, info = hill_climbing(d, l_T, C_T, D_T, index_to_node, kk, solver,
                                                         anchors=anchors, quick="1", neighbors=nbrs)
        rounds.append(dict(k=kk, loss=float(obj), anchors=[index_to_node[i] for i in anchors],
                           p=[float(v) for v in p], x=[float(v) for v in x], ybar=float(y),
                           time=time.time() - t0))
    res = dict(tree=tname, k=k, seed=seed, noise=nname, n=len(index_to_leaf),
               true_anchors=[t[0] for t in true], true_p=[float(t[1]) for t in true],
               true_x=[float(t[2]) for t in true], parent=parent, rounds=rounds)
    with open(out + ".tmp", "w") as f:
        json.dump(res, f)
    os.replace(out + ".tmp", out)
    return out


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    jobs = [(t, k, s, nn) for t in TREE_FILES for k in KS for s in SEEDS for nn in NOISE]
    jobs.sort(key=lambda j: (j[2], list(TREE_FILES).index(j[0])))  # finish whole seeds first
    with Pool(NPROC) as pool:
        for i, o in enumerate(pool.imap_unordered(one, jobs)):
            print(i, len(jobs), o, flush=True)
