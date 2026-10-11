# AI-assisted (Claude), exploration code for CS581 project
"""Oracle replicate: MAGUS's subsets and backbones, but each subset alignment replaced by the TRUE alignment
restricted to that subset (all-gap columns dropped). Merging it with GCM isolates the tree cost of the
merge step (cross-subset homology) from the cost of subset-alignment errors (mostly between close relatives).

    python3 truesub.py KEY [...]   -> /opt/work/treecrit/reps/<KEY>_ts/ (then vote.py <rep> magus|es4)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from measures import read_fasta, W  # noqa: E402

for key in sys.argv[1:]:
    src = os.path.join(W, key)
    dst = os.path.join(W, key + "_ts")
    os.makedirs(os.path.join(dst, "inputs", "subalignments"), exist_ok=True)
    if not os.path.exists(os.path.join(dst, "inputs", "backbones")):
        os.symlink(os.path.join(src, "inputs", "backbones"), os.path.join(dst, "inputs", "backbones"))
    for f in ("true.fasta", "true_tree.nwk", "tree.nwk", "true.tree"):
        if os.path.exists(os.path.join(src, f)) and not os.path.exists(os.path.join(dst, f)):
            os.symlink(os.path.join(src, f), os.path.join(dst, f))
    true = read_fasta(os.path.join(src, "true.fasta"))
    sd = os.path.join(src, "inputs", "subalignments")
    for f in sorted(os.listdir(sd)):
        taxa = list(read_fasta(os.path.join(sd, f)))
        rows = [true[t] for t in taxa]
        keep = [j for j in range(len(rows[0])) if any(r[j] not in "-." for r in rows)]
        with open(os.path.join(dst, "inputs", "subalignments", f), "w") as out:
            for t, r in zip(taxa, rows):
                out.write(">{}\n{}\n".format(t, "".join(r[j] for j in keep)))
    print(dst)
