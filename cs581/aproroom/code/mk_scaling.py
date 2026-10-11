"""Inputs for the runtime-scaling experiment.
Taxa: DISCO species_1000 rep 01 (1000 species, estimated gene trees from 100 bp); keep the first NG gene
trees and induce them on random species subsets of size k (nested), dropping trees with < 4 leaves.
Writes /opt/data/scaling/taxa_{k}/{genes.trees,s_tree.trees}.  Usage: python mk_scaling.py NG"""
import os, random, sys
import treeswift
NG = int(sys.argv[1])
src = "/opt/data/disco/trees/species_1000/01"
out = "/opt/data/scaling"
st = treeswift.read_tree_newick(open(f"{src}/s_tree.trees").read().strip())
sp = sorted(l.label for l in st.traverse_leaves())
rng = random.Random(1)
rng.shuffle(sp)
genes = [l.strip() for l in open(f"{src}/g_100.trees") if ";" in l][:NG]
for k in (50, 100, 200, 500, 1000):
    keep = set(sp[:k])
    d = f"{out}/taxa_{k}"
    os.makedirs(d, exist_ok=True)
    s2 = st.extract_tree_with(keep) if k < len(sp) else st
    s2.suppress_unifurcations()
    open(f"{d}/s_tree.trees", "w").write(s2.newick() + "\n")
    with open(f"{d}/genes.trees", "w") as f:
        for g in genes:
            t = treeswift.read_tree_newick(g)
            labs = [l.label for l in t.traverse_leaves() if l.label.split("_")[0] in keep]
            if len(labs) < 4:
                continue
            if len(labs) < t.num_nodes(internal=False):
                t = t.extract_tree_with(set(labs))
                t.suppress_unifurcations()
            f.write(t.newick() + "\n")
    print(k, sum(1 for _ in open(f"{d}/genes.trees")))
