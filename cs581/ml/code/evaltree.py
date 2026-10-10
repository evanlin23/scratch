"""Log-likelihood of a fixed topology under the same model RAxML-NG searches with (GTR+G),
optimising branch lengths and model parameters (raxml-ng --evaluate). Used for the "search gap":
lnL(true tree) - lnL(best tree found). Positive = the search missed a better-scoring, more accurate
tree (search headroom); negative = the ML criterion itself prefers wrong trees (no search headroom).
Model-agnostic: only the --model string is data-type specific.

    python evaltree.py ALN.clean.fasta TREE [--model GTR+G]
prints the log-likelihood. Polytomies are resolved randomly first (a lower bound on the best resolution).
"""
import argparse
import os
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runtrees as rt  # noqa: E402


def binary(tree, out):
    """Resolve polytomies randomly (zero-length edges: likelihood can only stay or improve)."""
    import random
    import dendropy
    t = dendropy.Tree.get(path=tree, schema="newick", preserve_underscores=True)
    t.resolve_polytomies(limit=2, update_bipartitions=False, rng=random.Random(1))
    t.write(path=out, schema="newick", suppress_rooting=True, suppress_internal_node_labels=True)


def evaluate(aln, tree, model="GTR+G"):
    work = tempfile.mkdtemp(prefix="eval_")
    binary(tree, os.path.join(work, "bin.tre"))
    tree = os.path.join(work, "bin.tre")
    subprocess.run([rt.RAXMLNG, "--evaluate", "--msa", aln, "--tree", tree, "--model", model, "--threads", "1",
                    "--seed", "1", "--prefix", os.path.join(work, "ev"), "--nofiles", "interim"],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for line in open(os.path.join(work, "ev.raxml.log")):
        if line.startswith("Final LogLikelihood:"):
            return float(line.split(":")[1])
    raise RuntimeError("evaluate failed: " + work)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("aln")
    ap.add_argument("tree")
    ap.add_argument("--model", default="GTR+G")
    a = ap.parse_args()
    print(evaluate(a.aln, a.tree, a.model))
