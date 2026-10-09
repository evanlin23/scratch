"""Consensus MSA from a collection of alignments of the same sequences (GCM with soft constraints).

    python -m gcmx.consensus OUT.fasta WORKDIR ALN1 ALN2 [ALN3 ...] [--group-size 13] [--primary auto|INDEX]

1. Primary input: the input with the highest mean SP agreement with the other
   inputs ("most central"; needs no reference alignment).
2. Groups: average-linkage clustering of the sequences on p-distances from the
   primary alignment, into about n / group-size groups.
3. Constraints: each group's alignment induced from the primary (all-gap
   columns removed). Groups are small, so these are "soft" constraints: GCM
   can re-decide homology between groups.
4. Evidence: every input alignment is a GCM backbone; GCM (MCL + MAGUS's
   minclusters trace) merges all groups at once.
"""

import argparse
import itertools
import json
import os
import shutil
import subprocess
import sys

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from . import fasta, score

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def p_distance_matrix(aln):
    names = list(aln)
    rows = np.array([np.frombuffer(aln[n].upper().encode(), dtype=np.uint8) for n in names])
    gap = (rows == ord("-")) | (rows == ord("."))
    present = (~gap).astype(np.float32)
    shared = present @ present.T
    matches = np.zeros_like(shared)
    for letter in np.unique(rows[~gap]):
        x = (rows == letter).astype(np.float32)
        matches += x @ x.T
    with np.errstate(invalid="ignore", divide="ignore"):
        dist = np.where(shared > 0, 1.0 - matches / np.maximum(shared, 1), 1.0)
    np.fill_diagonal(dist, 0.0)
    return names, dist


def agreement(a, b):
    s = score.fastsp(a, b)  # symmetric: shared / mean(homologies) = 1 - (SPFN + SPFP) / 2 approximately
    return 1.0 - s["avgErr"]


def build(out, workdir, inputs, group_size=13, primary="auto", threads=1, weight_power=0.0):
    os.makedirs(workdir, exist_ok=True)
    evidence = os.path.join(workdir, "evidence")
    os.makedirs(evidence, exist_ok=True)
    alns = []
    for i, path in enumerate(inputs):
        aln = fasta.upper(fasta.read(path))
        dest = os.path.join(evidence, "input_{}.fasta".format(i))
        fasta.write(aln, dest)
        alns.append((dest, aln))

    centrality = None
    if primary == "auto" or weight_power:
        k = len(alns)
        agree = np.zeros((k, k))
        for i, j in itertools.combinations(range(k), 2):
            agree[i, j] = agree[j, i] = agreement(alns[i][0], alns[j][0])
        centrality = agree.sum(axis=1) / max(k - 1, 1)
    p = int(np.argmax(centrality)) if primary == "auto" else int(primary)
    env = dict(os.environ)
    if weight_power:
        # M-Coffee-style: trust inputs that agree more with the others (no reference needed)
        w = (centrality / centrality.max()) ** weight_power
        weights_path = os.path.join(workdir, "weights.json")
        json.dump({os.path.basename(path): float(x) for (path, _), x in zip(alns, w)}, open(weights_path, "w"))
        env["GCMX_BACKBONE_WEIGHTS"] = weights_path
    names, dist = p_distance_matrix(alns[p][1])
    n_groups = max(1, round(len(names) / group_size))
    labels = fcluster(linkage(squareform(dist, checks=False), method="average"), t=n_groups, criterion="maxclust")

    groups_dir = os.path.join(workdir, "groups")
    shutil.rmtree(groups_dir, ignore_errors=True)
    os.makedirs(groups_dir)
    members = {}
    for name, label in zip(names, labels):
        members.setdefault(label, []).append(name)
    for label, taxa in members.items():
        fasta.write(fasta.restrict(alns[p][1], taxa), os.path.join(groups_dir, "group_{}.fasta".format(label)))

    magus_dir = os.path.join(workdir, "gcm")
    shutil.rmtree(magus_dir, ignore_errors=True)
    with open(os.path.join(workdir, "gcm.log"), "w") as log:
        subprocess.run([sys.executable, "-m", "gcmx.run_magus", "-np", str(threads), "-d", magus_dir,
                        "-s", groups_dir, "-b", evidence, "-o", out], cwd=CODE, stdout=log, stderr=subprocess.STDOUT,
                       check=True, env=env)
    return {"primary": p, "centrality": None if centrality is None else centrality.round(4).tolist(),
            "groups": len(members)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("out")
    parser.add_argument("workdir")
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--group-size", type=int, default=13)
    parser.add_argument("--primary", default="auto")
    parser.add_argument("--weight-power", type=float, default=0.0,
                        help="weight inputs by (centrality / max centrality) ** power; 0 = equal weights")
    args = parser.parse_args()
    print(json.dumps(build(args.out, args.workdir, args.inputs, args.group_size, args.primary,
                           weight_power=args.weight_power)))


if __name__ == "__main__":
    main()
