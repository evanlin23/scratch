"""How often are random placements (k distinct uniformly random edges) adjacent / closed-claw on real trees?"""
import sys, collections
import numpy as np
from treeswift import read_tree_newick
T = sys.argv[1]
files = {"birds-jarvis (48)": "birds-jarvis/birds-jarvis.tre.rerooted", "1kp (103)": "1kp/1kp-concat-fig2.tre",
         "hemipteroid (193)": "hemipteroid-johnson/hemipteroid-johnson.tre.rerooted"}
rng = np.random.default_rng(0)
print("| tree | k | P(some adjacent pair) | P(closed claw) | P(any claw, open or closed) |")
print("|---|---|---|---|---|")
for name, f in files.items():
    t = read_tree_newick(open(f"{T}/{f}").read().strip().split("\n")[0])
    nodes = list(t.traverse_preorder()); idx = {nd: i for i, nd in enumerate(nodes)}
    edges = [(idx[nd], idx[nd.parent]) for nd in nodes if nd.parent is not None]
    nbr = collections.defaultdict(set)
    for a, b in edges:
        nbr[a].add(b); nbr[b].add(a)
    internal = [v for v in nbr if len(nbr[v]) >= 3]
    for k in (3, 5, 10, 20, 50):
        if 2 * k > len(edges):
            continue
        adj = cl = anycl = 0; R = 2000
        for _ in range(R):
            S = rng.choice(len(edges), k, replace=False)
            VS = [x for e in S for x in edges[e]]
            sV = set(VS)
            adj += len(sV) < len(VS)
            cl += any(nbr[v] <= sV and v in sV for v in internal)
            anycl += any(nbr[v] <= sV for v in internal)
        print(f"| {name} | {k} | {adj/R:.3f} | {cl/R:.3f} | {anycl/R:.3f} |")
