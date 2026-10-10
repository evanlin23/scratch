"""IQ-TREE 3 (-m LG+G4 --fast, 1 thread, seed 1) on alignments of one replicate; nRF vs the true tree.

    python iqtrees.py TREEDIR TRUE_TREE OUT.jsonl [method ...]      (DendroPy python, e.g. pasta183)

TREEDIR is a trees.sh directory (true.fasta, magus.fasta, recipe.fasta, ...). Restartable.
"""
import json, os, subprocess, sys, time
import dendropy
from dendropy.calculate import treecompare

work, true_tree, out = sys.argv[1], sys.argv[2], sys.argv[3]
methods = sys.argv[4:]
done = {(r["dataset"], r["method"]) for r in map(json.loads, open(out))} if os.path.exists(out) else set()
name = os.path.basename(work.rstrip("/")).rsplit("_d", 1)[0]
for m in methods:
    aln = os.path.join(work, m + ".fasta")
    if (name, m) in done or not os.path.exists(aln):
        continue
    pre = os.path.join(work, "iq_" + m)
    start = time.time()
    subprocess.run(["/opt/mm/root/envs/bio/bin/iqtree3", "-s", aln, "-m", "LG+G4", "--fast", "-T", "1", "-seed", "1",
                    "--prefix", pre, "-redo", "-quiet"], check=True)
    wall = round(time.time() - start, 1)
    tns = dendropy.TaxonNamespace()
    t = dendropy.Tree.get(path=true_tree, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    e = dendropy.Tree.get(path=pre + ".treefile", schema="newick", taxon_namespace=tns, preserve_underscores=True)
    t.is_rooted = e.is_rooted = False
    t.encode_bipartitions(); e.encode_bipartitions()
    fp, fn = treecompare.false_positives_and_negatives(t, e)
    nint = len(t.leaf_nodes()) - 3
    ti = len([b for b in t.bipartition_encoding if not b.is_trivial()])
    ei = len([b for b in e.bipartition_encoding if not b.is_trivial()])
    lnl = None
    for l in open(pre + ".iqtree"):
        if l.startswith("Log-likelihood of the tree"):
            lnl = float(l.split(":")[1].split()[0])
    row = {"dataset": name, "method": m, "FN": round(fn / ti, 4), "FP": round(fp / max(ei, 1), 4),
           "RF": round((fn + fp) / (2 * nint), 4), "lnL": lnl, "iqtree_wall": wall}
    with open(out, "a") as f:
        f.write(json.dumps(row) + "\n")
    print(json.dumps(row), flush=True)
    for ext in (".ckp.gz", ".mldist", ".bionj", ".log", ".model.gz", ".uniqueseq.phy"):
        if os.path.exists(pre + ext):
            os.remove(pre + ext)
