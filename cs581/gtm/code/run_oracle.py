"""Reachability: constrained SPR hill-climb from GTM scored by RF to the TRUE tree.
Usage: python3 run_oracle.py COND REP GUIDE radius"""
import sys, time
from datasets import inputs, read_fasta
from phylo import read_tree, fn_fp, is_induced, parse_newick
from blend import BlendSearch, compress, run_oracle
import numpy as np
cond, rep, guide, radius = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
start = sys.argv[5] if len(sys.argv) > 5 else None
x = inputs(cond, rep, guide)
T = read_tree(x["true"]); G = read_tree(start or x["published_gtm"])
subs = [read_tree(p) for p in x["subsets"]]
names = sorted(G.label.values())
states = np.zeros((len(names), 1), dtype=np.uint8) + 1; w = np.ones(1, dtype=np.int64)
bs = BlendSearch(G, names, states, w, [[t.label[v] for v in t.leaves()] for t in subs], radius=radius)
from phylo import bipartitions
bs.oracle = bipartitions(T, bs.idx)
t0 = time.time()
m = run_oracle(bs, log=lambda m, rf, d: print(m, rf, f"{time.time()-t0:.0f}s", flush=True) if m % 20 == 0 else None)
B = parse_newick(bs.to_newick())
print(cond, rep, guide, "moves", m, "FN gtm %.4f -> oracle-search %.4f" % (fn_fp(G, T)[0], fn_fp(B, T)[0]),
      "constraints_ok", all(is_induced(B, t) for t in subs))
