"""Pilot: disjoint-tree-merger divide-and-conquer ML (GTM pipeline, Smirnov & Warnow 2020).

    python dnc.py DATASET REP ALN MAXSIZE [--polish]

1. guide tree = FastTree (-gtr -gamma) tree on ALN (from runtrees.py; its time is added);
2. centroid-edge decomposition of the guide tree into disjoint subsets of <= MAXSIZE leaves;
3. RAxML-ng (GTR+G, one parsimony start) on each subset alignment, 4 subsets in parallel;
4. GTM merges the subset trees using the guide tree (git clone of vlasmirnov/GTM at /opt/tools/GTM);
5. --polish: RAxML-ng search started from the merged tree (lets the merge "blend").
Appends a row (method "dnc<MAXSIZE>" / "dnc<MAXSIZE>_polish") to results/dnc.jsonl.
"""
import json
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

import dendropy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runtrees  # noqa: E402
import treeerr  # noqa: E402

GTM = "/opt/tools/GTM/gtm.py"


def decompose(tree, maxsize):
    """Recursive centroid-edge decomposition; returns a list of leaf-label sets."""
    labels = [l.taxon.label for l in tree.leaf_node_iter()]
    if len(labels) <= maxsize:
        return [labels]
    tree.is_rooted = True
    counts = {}
    for nd in tree.postorder_node_iter():
        counts[nd] = 1 if nd.is_leaf() else sum(counts[c] for c in nd.child_node_iter())
    n = len(labels)
    best = min((nd for nd in tree.preorder_node_iter() if nd is not tree.seed_node),
               key=lambda nd: abs(n - 2 * counts[nd]))
    side = {l.taxon.label for l in best.leaf_iter()}
    out = []
    for part in (side, set(labels) - side):
        sub = tree.extract_tree_with_taxa_labels(part)
        out += decompose(sub, maxsize)
    return out


def main(ds, rep, aln, maxsize, polish):
    d = os.path.join(runtrees.MLDATA, ds, "R%s" % rep)
    clean = os.path.join(d, aln + ".clean.fasta")
    guide = os.path.join(d, "trees", "%s.fasttree.tre" % aln)
    rows = [json.loads(l) for l in open(os.path.join(os.path.dirname(__file__), "..", "results", "baseline.jsonl"))]
    ft_sec = [r["seconds"] for r in rows if (r["dataset"], str(r["rep"]), r["aln"], r["method"]) == (ds, str(rep), aln, "fasttree")][0]
    work = tempfile.mkdtemp(prefix="dnc_%s_%s_%s_" % (ds, rep, aln))
    t0 = time.time()
    tree = dendropy.Tree.get(path=guide, schema="newick", preserve_underscores=True)
    subsets = decompose(tree, maxsize)
    names, seqs = runtrees.read_fasta(clean)
    seqd = dict(zip(names, seqs))
    t_dec = time.time() - t0

    def sub(i):
        sd = os.path.join(work, "s%d" % i)
        os.makedirs(sd)
        a = os.path.join(sd, "aln.fasta")
        runtrees.write_fasta(a + ".raw", subsets[i], [seqd[x] for x in subsets[i]])
        runtrees.clean_alignment(a + ".raw", a)
        out = os.path.join(sd, "tree.tre")
        runtrees.estimate("raxmlng", a, out, sd)
        return out
    t1 = time.time()
    with ThreadPoolExecutor(4) as ex:
        subtrees = list(ex.map(sub, range(len(subsets))))
    t_sub = time.time() - t1  # wall time with 4 concurrent single-thread jobs
    t2 = time.time()
    merged = os.path.join(d, "trees", "%s.dnc%d.tre" % (aln, maxsize))
    subprocess.run([sys.executable, GTM, "-s", guide, "-t"] + subtrees + ["-o", merged], check=True,
                   stdout=subprocess.DEVNULL)
    t_merge = time.time() - t2
    sizes = sorted(len(s) for s in subsets)
    res = []
    err = treeerr.error(os.path.join(d, "true_tree.tre"), merged)
    res.append({"dataset": ds, "rep": rep, "aln": aln, "method": "dnc%d" % maxsize, "subsets": sizes,
                "seconds_wall": round(ft_sec + t_dec + t_sub + t_merge, 1), "seconds_fasttree": ft_sec,
                "seconds_subsets_wall": round(t_sub, 1), "seconds_merge": round(t_merge, 1),
                **{k: err[k] for k in ("fn_rate", "fp_rate", "rf_rate")}, "tree": merged})
    if polish:
        pt = os.path.join(d, "trees", "%s.dnc%d_polish.tre" % (aln, maxsize))
        pw = os.path.join(work, "polish")
        os.makedirs(pw)
        sec, lnl = runtrees.estimate("raxmlng_ft", clean, pt, pw, start_tree=merged)
        err = treeerr.error(os.path.join(d, "true_tree.tre"), pt)
        res.append({**res[0], "method": "dnc%d_polish" % maxsize, "seconds_polish": round(sec, 1),
                    "seconds_wall": round(res[0]["seconds_wall"] + sec, 1), "lnl_tool": lnl,
                    **{k: err[k] for k in ("fn_rate", "fp_rate", "rf_rate")}, "tree": pt})
    with open(os.path.join(os.path.dirname(__file__), "..", "results", "dnc.jsonl"), "a") as f:
        for r in res:
            f.write(json.dumps(r) + "\n")
            print(json.dumps(r), flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), "--polish" in sys.argv)
