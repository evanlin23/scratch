import json, sys, glob, numpy as np
from scipy.sparse.csgraph import shortest_path
from idclass import URTree
D0 = "/tmp/claude-0/-home-user-scratch/78aaac8f-d2b2-5da5-a76a-2c029e089703/scratchpad/DecoDiPhy-Data/biotrees"
bad = 0; tot = 0
for f in sorted(glob.glob(D0 + "/birds-jarvis/*noise0"))[:40]:
    R = URTree(open(f + "/pruned_tree.trees").read().strip().split("\n")[0])
    rounds = json.load(open(f + "/all_rounds.json"))
    tr = rounds[0]; fin = [r for r in rounds if r["rounds"] == "final"][-1]
    Dm = shortest_path(R.L.copy().tocsr().multiply(0) + __import__("scipy").sparse.csr_matrix(
        ([R.adj[u][w] for u in range(R.N) for w in R.adj[u]], ([u for u in range(R.N) for w in R.adj[u]], [w for u in range(R.N) for w in R.adj[u]])), shape=(R.N, R.N)))
    def dvec(r):
        pls = [R.placement(a, x) for a, x in zip(r["anchors"], r["x"])]
        m = R.node_measure(pls, r["p"])
        return Dm[R.leaves] @ m + r["y"]
    e = np.abs(dvec(tr) - dvec(fin)).max(); tot += 1; bad += e > 1e-3
    print(f.split("/")[-1], tr["k"], fin["k"], f"{e:.2e}", fin["loss"])
print(bad, tot)
