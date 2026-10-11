# AI-assisted (Claude), exploration code for CS581 project
"""IQ-TREE (-m LG+G4 --fast -T 2 -seed 1) on true / magus / es4 / hard-bb of SIMHIGH replicates, NPAR runs at once (env NPAR, default 1:
two concurrent runs of ~5.5 GB each were OOM-killed by the sandbox memory cgroup);
nRF vs the true tree (as cs581/protbench/code/trees.py: (FN + FP) / (2 (n - 3))) -> OUT.jsonl.

    python iq.py WORK OUT.jsonl REP [REP ...]

WORK/reps/REP holds the unpacked bank replicate with vote/{magus,es4,hard-bb}/out.fasta (gcmvote run.py).
Restartable: rows already in OUT are skipped; an interrupted IQ-TREE run resumes from its checkpoint.
Run with a Python that has DendroPy (e.g. /opt/mm/root/envs/pasta183/bin/python).
"""
import json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
import dendropy
from dendropy.calculate import treecompare

IQ = "/opt/mm/root/envs/bio/bin/iqtree3"
METHODS = ["true", "magus", "es4", "hard-bb"]
work, out, reps = sys.argv[1], sys.argv[2], sys.argv[3:]
version = subprocess.run([IQ, "--version"], capture_output=True, text=True).stdout.splitlines()[0].strip()
done = {(r["rep"], r["method"]) for r in map(json.loads, open(out))} if os.path.exists(out) else set()


def nrf(true_tree, est):
    tns = dendropy.TaxonNamespace()
    t = dendropy.Tree.get(path=true_tree, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    e = dendropy.Tree.get(path=est, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    t.is_rooted = e.is_rooted = False
    t.encode_bipartitions()
    e.encode_bipartitions()
    fp, fn = treecompare.false_positives_and_negatives(t, e)
    nint = len(t.leaf_nodes()) - 3
    return round((fn + fp) / (2 * nint), 4), fn, fp


def job(rep, m):
    rd = os.path.join(work, "reps", rep)
    aln = os.path.join(rd, "true.fasta") if m == "true" else os.path.join(rd, "vote", m, "out.fasta")
    td = os.path.join(work, "iq", rep)
    os.makedirs(td, exist_ok=True)
    pre = os.path.join(td, m)
    resumed = os.path.exists(pre + ".ckp.gz")
    start = time.time()
    with open(pre + ".stdout", "a") as f:
        subprocess.run([IQ, "-s", aln, "-m", "LG+G4", "--fast", "-T", "2", "-seed", "1", "--prefix", pre],
                       stdout=f, stderr=subprocess.STDOUT, check=True)
    wall = round(time.time() - start, 1)
    r, fn, fp = nrf(os.path.join(rd, "true_tree.nwk"), pre + ".treefile")
    row = {"rep": rep, "method": m, "nRF": r, "FN": fn, "FP": fp, "wall": wall, "resumed": resumed,
           "over_90min": wall > 5400, "iqtree": version, "settings": "-m LG+G4 --fast -T 2 -seed 1"}
    return row


jobs = [(rep, m) for rep in reps for m in METHODS if (rep, m) not in done]
left = {rep: sum(1 for j in jobs if j[0] == rep) for rep in reps}
with ThreadPoolExecutor(int(os.environ.get("NPAR", "1"))) as ex:
    futs = [(j, ex.submit(job, *j)) for j in jobs]
    for (rep, m), fu in futs:
        try:
            row = fu.result()
        except Exception as exc:  # e.g. OOM kill: report and continue; a rerun resumes from the checkpoint
            print("FAILED", rep, m, exc, flush=True)
            left[rep] -= 1
            continue
        with open(out, "a") as f:
            f.write(json.dumps(row) + "\n")
        print(json.dumps(row), flush=True)
        left[rep] -= 1
        if left[rep] == 0:
            print("REPDONE", rep, flush=True)
print("ALLDONE", flush=True)
