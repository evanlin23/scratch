"""FastTree (-lg -gamma) on each alignment of a simulated dataset; FN/FP/RF rates against the true tree.

    python trees.py WORKDIR TRUE_TREE OUT.jsonl [method ...]

WORKDIR is a basemeth.py dataset directory (true.fasta, merge_<method>.fasta). Adapted from the protbench branch.
Run with a Python that has DendroPy (e.g. /opt/mm/root/envs/pasta183/bin/python).
"""
import json, os, subprocess, sys, time
import dendropy
from dendropy.calculate import treecompare

work, true_tree, out = sys.argv[1], sys.argv[2], sys.argv[3]
methods = sys.argv[4:] or ["true", "merge_linsi"]
done = set()
if os.path.exists(out):
    done = {(r["dataset"], r["method"]) for r in map(json.loads, open(out))}
name = os.path.basename(work.rstrip("/"))
for m in methods:
    aln = os.path.join(work, m + ".fasta")
    if (name, m) in done or not os.path.exists(aln):
        continue
    tre = os.path.join(work, m + ".fasttree.nwk")
    start = time.time()
    with open(tre, "w") as f, open(tre + ".log", "w") as e:
        subprocess.run(["FastTree", "-lg", "-gamma", "-quiet", aln], stdout=f, stderr=e, check=True)
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
