# AI-assisted (Claude), exploration code for CS581 project
"""Copied from h10_trees.py on claude/cs581-gcmvote-h10.
FastTree (-gtr -gamma -nt) on named alignments; FN/FP/RF (nRF) against the true tree.
Adapted from cs581/protbench/code/trees.py (which uses -lg for proteins).

    python h10_trees.py DATASET TRUE_TREE OUT.jsonl METHOD=ALN.fasta [...]

Run with a Python that has DendroPy (e.g. /opt/mm/root/envs/pasta183/bin/python).
"""
import json, os, subprocess, sys, time
import dendropy
from dendropy.calculate import treecompare

name, true_tree, out = sys.argv[1:4]
done = set()
if os.path.exists(out):
    done = {(r["dataset"], r["method"]) for r in map(json.loads, open(out))}
for arg in sys.argv[4:]:
    m, aln = arg.split("=", 1)
    if (name, m) in done or not os.path.exists(aln):
        continue
    tre = aln + ".fasttree.nwk"
    start = time.time()
    with open(tre, "w") as f, open(tre + ".log", "w") as e:
        subprocess.run(["FastTree", "-gtr", "-gamma", "-nt", "-quiet", aln], stdout=f, stderr=e, check=True)
    wall = round(time.time() - start, 1)
    tns = dendropy.TaxonNamespace()
    t = dendropy.Tree.get(path=true_tree, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    e = dendropy.Tree.get(path=tre, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    t.is_rooted = e.is_rooted = False
    t.encode_bipartitions()
    e.encode_bipartitions()
    fp, fn = treecompare.false_positives_and_negatives(t, e)
    nint = len(t.leaf_nodes()) - 3
    ti = len([b for b in t.bipartition_encoding if not b.is_trivial()])
    ei = len([b for b in e.bipartition_encoding if not b.is_trivial()])
    row = {"dataset": name, "method": m, "FN": round(fn / ti, 4), "FP": round(fp / max(ei, 1), 4),
           "RF": round((fn + fp) / (2 * nint), 4), "true_internal": ti, "est_internal": ei, "fasttree_wall": wall}
    with open(out, "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)
