"""Achievable blending headroom on published data: oracle-guided constrained insertion
followed by oracle-scored constrained SPR (radius 6), best of that and oracle SPR from GTM.
Usage: python3 oracle_pub.py COND OUT_TSV"""
import sys
import numpy as np
from datasets import REPS, inputs, read_fasta
from phylo import read_tree, fn_fp, parse_newick, bipartitions, is_induced
from blend import BlendSearch, run_oracle
import insert
cond, out = sys.argv[1], sys.argv[2]
rows = []
for rep in REPS[cond]:
    for g in ("FT", "IQ"):
        x = inputs(cond, rep, g)
        T = read_tree(x["true"]); G = read_tree(x["published_gtm"])
        subs = [read_tree(p) for p in x["subsets"]]
        names = sorted(G.label.values())
        st, w = np.ones((len(names), 1), np.uint8), np.ones(1, np.int64)
        sn = [[t.label[v] for v in t.leaves()] for t in subs]
        res = []
        for start in ("gtm", "ins"):
            S = G if start == "gtm" else parse_newick(insert.run(names, st, w, subs, "oracle", ref_tree=T).to_newick())
            bs = BlendSearch(S, names, st, w, sn, radius=6)
            bs.oracle = bipartitions(T, bs.idx)
            run_oracle(bs)
            B = parse_newick(bs.to_newick())
            assert all(is_induced(B, t) for t in subs)
            res.append(fn_fp(B, T)[0])
        rows.append((cond, rep, g, fn_fp(G, T)[0], res[0], res[1], min(res)))
        print(*rows[-1], sep="\t", flush=True)
with open(out, "w") as f:
    f.write("condition\treplicate\tguide\tFN_GTM\tFN_oracleSPR_from_GTM\tFN_oracle_insert+SPR\tFN_oracle_best\n")
    for r in rows:
        f.write("%s\t%s\t%s\t%.5f\t%.5f\t%.5f\t%.5f\n" % r)
