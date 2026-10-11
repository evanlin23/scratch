# AI-assisted (Claude), exploration code for CS581 project
"""FastTree on alignment files; normalized RF to the true tree (same settings as protbench/code/trees.py).

    python3 trees.py JOBS.tsv OUT.jsonl [--lanes 4]

JOBS.tsv lines: key<TAB>label<TAB>alignment path. Proteins: FastTree -lg -gamma; 1000M*: -nt -gtr -gamma.
Trees go to /opt/work/treecrit/trees/<key>/<label>.nwk. Restartable: (key, label) rows in OUT are skipped.
"""
import argparse
import json
import os
import subprocess
import time
from multiprocessing import Pool

W = "/opt/work/treecrit"


def tree_file(rep):
    for f in ("true_tree.nwk", "tree.nwk", "true.tree"):
        if os.path.exists(os.path.join(rep, f)):
            return os.path.join(rep, f)


def rf(true_tree, est_tree):
    import dendropy
    from dendropy.calculate import treecompare
    tns = dendropy.TaxonNamespace()
    t = dendropy.Tree.get(path=true_tree, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    e = dendropy.Tree.get(path=est_tree, schema="newick", taxon_namespace=tns, preserve_underscores=True)
    t.is_rooted = e.is_rooted = False
    t.encode_bipartitions(); e.encode_bipartitions()
    fp, fn = treecompare.false_positives_and_negatives(t, e)
    nint = len(t.leaf_nodes()) - 3
    ti = len([b for b in t.bipartition_encoding if not b.is_trivial()])
    ei = len([b for b in e.bipartition_encoding if not b.is_trivial()])
    return dict(FN=round(fn / ti, 4), FP=round(fp / max(ei, 1), 4), RF=round((fn + fp) / (2 * nint), 4))


def job(a):
    key, label, aln = a
    d = os.path.join(W, "trees", key)
    os.makedirs(d, exist_ok=True)
    tre = os.path.join(d, label + ".nwk")
    flags = ["-nt", "-gtr", "-gamma"] if key.startswith("1000") else ["-lg", "-gamma"]
    start = time.time()
    with open(tre + ".tmp", "w") as f, open(tre + ".log", "w") as e:
        r = subprocess.run(["nice", "-n", "5", "fasttree", "-quiet", "-nosupport"] + flags + [aln], stdout=f, stderr=e)
    if r.returncode:
        return {"key": key, "label": label, "error": r.returncode}
    os.rename(tre + ".tmp", tre)
    row = {"key": key, "label": label, "aln": aln, "fasttree_wall": round(time.time() - start, 1),
           **rf(tree_file(os.path.join(W, "reps", key)), tre)}
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("jobs")
    ap.add_argument("out")
    ap.add_argument("--lanes", type=int, default=4)
    a = ap.parse_args()
    jobs = [tuple(l.rstrip("\n").split("\t")) for l in open(a.jobs) if l.strip() and not l.startswith("#")]
    done = set()
    if os.path.exists(a.out):
        done = {(r["key"], r["label"]) for r in map(json.loads, open(a.out)) if "error" not in r}
    jobs = [j for j in jobs if (j[0], j[1]) not in done]
    print(len(jobs), "trees", flush=True)
    with Pool(a.lanes) as p, open(a.out, "a") as f:
        for row in p.imap_unordered(job, jobs):
            f.write(json.dumps(row) + "\n"); f.flush()
            print(json.dumps(row), flush=True)


if __name__ == "__main__":
    main()
