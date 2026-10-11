"""GTM-Blend-ML (prior pilot, cs581/gtm/code/mlspr.py) on a pipeline directory.
Usage: python3 run_blendml.py ALN GTM_TREE OUT_TREE SUBTREE...   (env MLSPR_MAX_ROUNDS, default 10)"""
import os
import sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "../../gtm/code"))
from phylo import read_tree  # noqa: E402
from datasets import read_fasta  # noqa: E402
import mlspr  # noqa: E402
aln, gtm, out, *subf = sys.argv[1:]
seqs = read_fasta([aln]); names = sorted(seqs)
subs = [read_tree(p) for p in subf]
logf = open(out + ".log", "w")
bs, h = mlspr.run(read_tree(gtm), names, seqs, [[t.label[v] for v in t.leaves()] for t in subs], out + "_work",
                  radius=4, model_tree=gtm, max_rounds=int(os.environ.get("MLSPR_MAX_ROUNDS", "10")),
                  log=lambda s: (logf.write(s + "\n"), logf.flush()))
open(out, "w").write(bs.to_newick() + "\n")
logf.write("logL path %s\n" % h)
