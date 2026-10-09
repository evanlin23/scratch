"""Soft constraints: split each constraint alignment into smaller groups.

    python -m gcmx.split SUBALIGNMENT_DIR OUT_DIR GROUPS_PER_SUBSET [--method linkage|random] [--seed S]

MAGUS treats each subset alignment as a hard constraint, so an error inside
a subset alignment propagates to every cross-subset pair. Splitting each
subset alignment into m groups (induced sub-alignments, all-gap columns
removed) and merging all k*m groups at once lets GCM re-decide homology
between groups of the same subset. Pass the original subset alignments as
extra backbones so their evidence is kept, but can be outvoted:

    magus -s OUT_DIR -b BACKBONE_DIR SUBALIGNMENT_DIR ...

m = 1 reproduces MAGUS; m = subset size is unconstrained MAGUS.
Groups are formed by average-linkage clustering on p-distances computed
from the subset alignment (so each group holds similar sequences).
"""

import argparse
import os
import random

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from . import fasta


def p_distances(aln):
    names = list(aln)
    codes = np.array([[ord(c) for c in aln[n].upper()] for n in names], dtype=np.int16)
    gap = (codes == ord("-")) | (codes == ord("."))
    n = len(names)
    dist = np.zeros((n, n))
    for i in range(n):
        both = ~gap[i] & ~gap
        shared = both.sum(axis=1)
        diff = ((codes[i] != codes) & both).sum(axis=1)
        dist[i] = np.where(shared > 0, diff / np.maximum(shared, 1), 1.0)
    return names, (dist + dist.T) / 2


def split_groups(aln, m, method, rng):
    names = list(aln)
    if m <= 1 or len(names) <= 1:
        return [names]
    m = min(m, len(names))
    if method == "random":
        shuffled = names[:]
        rng.shuffle(shuffled)
        return [shuffled[i::m] for i in range(m)]
    names, dist = p_distances(aln)
    np.fill_diagonal(dist, 0)
    labels = fcluster(linkage(squareform(dist, checks=False), method="average"), t=m, criterion="maxclust")
    groups = {}
    for name, label in zip(names, labels):
        groups.setdefault(label, []).append(name)
    return list(groups.values())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("subalignments")
    parser.add_argument("outdir")
    parser.add_argument("groups", type=int)
    parser.add_argument("--method", default="linkage", choices=("linkage", "random"))
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)
    total = 0
    for filename in sorted(os.listdir(args.subalignments)):
        aln = fasta.read(os.path.join(args.subalignments, filename))
        for g, taxa in enumerate(split_groups(aln, args.groups, args.method, rng)):
            fasta.write(fasta.restrict(aln, taxa), os.path.join(args.outdir, "{}_g{}.txt".format(filename[:-4], g)))
            total += 1
    print("wrote {} group alignments".format(total))


if __name__ == "__main__":
    main()
