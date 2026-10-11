# AI-assisted (Claude), exploration code for CS581 project
"""IQ-TREE 3 (-m LG+G4 --fast -T 2 -seed 1) on one alignment; nRF vs the true tree (as protbench trees.py).

    python iq_b.py REP METHOD ALN TRUE_TREE WORKDIR OUT.jsonl     (DendroPy python, e.g. pasta183)

Appends {rep, method, nRF, FN, FP, lnL, wall, iqtree} to OUT.jsonl. Restartable: (rep, method) rows already
in OUT.jsonl are skipped; an unfinished IQ-TREE run is restarted from scratch (-redo).
"""
import fcntl, json, os, subprocess, sys, time
import dendropy
from dendropy.calculate import treecompare

IQ = "/opt/mm/root/envs/bio/bin/iqtree3"
rep, m, aln, true_tree, work, out = sys.argv[1:7]
if os.path.exists(out) and any((r["rep"], r["method"]) == (rep, m) for r in map(json.loads, open(out))):
    sys.exit(0)
ver = subprocess.run([IQ, "--version"], capture_output=True, text=True).stdout.split("\n")[0].strip()
os.makedirs(work, exist_ok=True)
pre = os.path.join(work, "iq_" + m)
start = time.time()
subprocess.run([IQ, "-s", aln, "-m", "LG+G4", "--fast", "-T", "2", "-seed", "1", "--prefix", pre, "-redo", "-quiet"],
               check=True)
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
row = {"rep": rep, "method": m, "nRF": round((fn + fp) / (2 * nint), 4), "FN": round(fn / ti, 4),
       "FP": round(fp / max(ei, 1), 4), "lnL": lnl, "wall": wall, "over_90min": wall > 5400, "iqtree": ver}
with open(out, "a") as f:
    fcntl.flock(f, fcntl.LOCK_EX)
    f.write(json.dumps(row) + "\n")
print(json.dumps(row), flush=True)
