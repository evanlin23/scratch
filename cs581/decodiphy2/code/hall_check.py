"""Is the continuum condition (B) [L[V\\W, I] rank-deficient, W = V(S)] purely combinatorial?
Compare numerical rank (unit + 3 random length draws) with the structural rank (max matching in the
bipartite pattern graph rows V\\W x columns I, edge iff u == v or u ~ v), i.e. the Hall condition
'exists U subset I with |N[U] \\ W| < |U|'.  usage: python hall_check.py NMAX OUT"""
import sys, os, itertools, collections
import numpy as np, networkx as nx
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "decodiphy", "code"))
from config_enum import shapes, laplacian
nmax, out = int(sys.argv[1]), sys.argv[2]
rng = np.random.default_rng(3)
c = collections.Counter()
for n in range(4, nmax + 1):
    for E in shapes(n):
        N = 2 * n - 2
        Ls = [laplacian(N, E, np.ones(len(E)))] + [laplacian(N, E, rng.exponential(1, len(E)) + 0.01) for _ in range(3)]
        nbr = collections.defaultdict(set)
        for a, b in E: nbr[a].add(b); nbr[b].add(a)
        for k in range(1, 6):
            for S in itertools.combinations(range(len(E)), k):
                W = {x for e in S for x in E[e]}
                rows = [u for u in range(N) if u not in W]
                num = [np.linalg.matrix_rank(L[rows][:, n:], tol=1e-9) < N - n for L in Ls]
                G = nx.Graph()
                G.add_nodes_from([("c", v) for v in range(n, N)])
                for v in range(n, N):
                    for u in nbr[v] | {v}:
                        if u not in W: G.add_edge(("c", v), ("r", u))
                mm = len(nx.bipartite.maximum_matching(G, top_nodes=[("c", v) for v in range(n, N)])) // 2
                struct = mm < N - n
                c["S"] += 1; c["deficient"] += struct
                c["agree_all"] += all(x == struct for x in num)
                if not all(x == struct for x in num): c[f"disagree n={n}"] += 1
open(out, "w").write("\n".join(f"{k}: {v}" for k, v in c.items()) + "\n")
print(c)
