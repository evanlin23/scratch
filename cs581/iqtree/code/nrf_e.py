# AI-assisted (Claude), exploration code for CS581 project
"""nRF of an IQ-TREE tree vs the true tree (same computation as cs581/protbench/code/trees.py); append a row.

    python nrf_e.py TRUE_TREE EST_TREE OUT.jsonl REP METHOD WALL VERSION
Run with a Python that has DendroPy (e.g. /opt/mm/root/envs/pasta183/bin/python).
"""
import json, sys
import dendropy
from dendropy.calculate import treecompare

true_tree, est, out, rep, method, wall, ver = sys.argv[1:8]
tns = dendropy.TaxonNamespace()
t = dendropy.Tree.get(path=true_tree, schema="newick", taxon_namespace=tns, preserve_underscores=True)
e = dendropy.Tree.get(path=est, schema="newick", taxon_namespace=tns, preserve_underscores=True)
assert len(tns) == len(t.leaf_nodes()) == len(e.leaf_nodes()), "taxon mismatch"
t.is_rooted = e.is_rooted = False
t.encode_bipartitions()
e.encode_bipartitions()
fp, fn = treecompare.false_positives_and_negatives(t, e)
nint = len(t.leaf_nodes()) - 3
ti = len([b for b in t.bipartition_encoding if not b.is_trivial()])
ei = len([b for b in e.bipartition_encoding if not b.is_trivial()])
row = {"rep": rep, "method": method, "nRF": round((fn + fp) / (2 * nint), 4), "FN": round(fn / ti, 4),
       "FP": round(fp / max(ei, 1), 4), "wall": float(wall), "iqtree": ver,
       "settings": "-m LG+G4 --fast -T 2 -seed 1"}
with open(out, "a") as f:
    f.write(json.dumps(row) + "\n")
print(json.dumps(row), flush=True)
