"""PASTA-style local backbones for GCM.

    python -m gcmx.local_backbones SUBALIGNMENT_DIR GUIDE_ALIGNMENT OUT_DIR
        [--neighbors 2] [--size 0] [--jobs 4] [--seed 0]

MAGUS's backbones each sample a few sequences from *every* subset, so every
backbone is an alignment of a very divergent set, which MAFFT aligns poorly.
PASTA instead aligned pairs of *adjacent* subsets (similar sequences, accurate
alignments) and merged them by transitivity. Here we compute such local
alignments and give them to GCM as additional backbones, so GCM still merges
all subset alignments at once, but with accurate local evidence.

Subset adjacency: the subsets' pairwise mean p-distance is estimated from
GUIDE_ALIGNMENT (any alignment of all sequences, e.g. MAGUS's own output);
pairs = maximum-similarity spanning tree plus each subset's `--neighbors`
nearest subsets. Each local backbone is MAFFT-L-INS-i (MAGUS's exact
command) on the union of the two subsets (`--size N`: at most N sequences
per subset, 0 = all).
"""

import argparse
import concurrent.futures
import itertools
import os
import random
import subprocess

import numpy as np
from scipy.sparse.csgraph import minimum_spanning_tree

from . import fasta
from .subset_bakeoff import MAGUS_MAFFT


def subset_distances(subsets, guide, rng, samples=8):
    names = {n: np.frombuffer(guide[n].upper().encode(), dtype=np.uint8) for s in subsets for n in s}
    gap = ord("-")
    k = len(subsets)
    dist = np.zeros((k, k))
    for i, j in itertools.combinations(range(k), 2):
        ds = []
        for a in rng.sample(subsets[i], min(samples, len(subsets[i]))):
            for b in rng.sample(subsets[j], min(samples, len(subsets[j]))):
                x, y = names[a], names[b]
                both = (x != gap) & (y != gap)
                if both.any():
                    ds.append(float((x[both] != y[both]).mean()))
        dist[i, j] = dist[j, i] = np.mean(ds) if ds else 1.0
    return dist


def choose_pairs(dist, neighbors):
    k = len(dist)
    mst = minimum_spanning_tree(dist + 1e-9).tocoo()  # +eps keeps zero-distance edges
    pairs = {tuple(sorted((int(i), int(j)))) for i, j in zip(mst.row, mst.col)}
    for i in range(k):
        order = [j for j in np.argsort(dist[i]) if j != i]
        pairs.update(tuple(sorted((i, int(j)))) for j in order[:neighbors])
    return sorted(pairs)


def align(seqs, out_path):
    tmp = out_path + ".in.fa"
    fasta.write(seqs, tmp)
    with open(out_path + ".tmp", "w") as out:
        subprocess.run([MAGUS_MAFFT, "--localpair", "--maxiterate", "1000", "--ep", "0.123", "--quiet",
                        "--thread", "1", "--anysymbol", tmp], stdout=out, stderr=subprocess.DEVNULL, check=True)
    aln = fasta.upper(fasta.read(out_path + ".tmp"))
    fasta.write(aln, out_path)
    os.remove(tmp)
    os.remove(out_path + ".tmp")
    return out_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("subalignments")
    parser.add_argument("guide")
    parser.add_argument("outdir")
    parser.add_argument("--neighbors", type=int, default=2)
    parser.add_argument("--size", type=int, default=0)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    files = sorted(os.listdir(args.subalignments))
    subs = [fasta.read(os.path.join(args.subalignments, f)) for f in files]
    subsets = [list(s) for s in subs]
    guide = fasta.read(args.guide)
    pairs = choose_pairs(subset_distances(subsets, guide, rng), args.neighbors)
    os.makedirs(args.outdir, exist_ok=True)
    print("{} subsets -> {} local backbones".format(len(subsets), len(pairs)))

    jobs = []
    with concurrent.futures.ThreadPoolExecutor(args.jobs) as pool:
        for i, j in pairs:
            seqs = {}
            for idx in (i, j):
                taxa = subsets[idx][:]
                if args.size:
                    rng.shuffle(taxa)
                    taxa = taxa[:args.size]
                seqs.update({t: subs[idx][t].replace("-", "").replace(".", "").upper() for t in taxa})
            jobs.append(pool.submit(align, seqs, os.path.join(args.outdir, "local_{}_{}.txt".format(i, j))))
        for job in jobs:
            job.result()
    print("wrote", len(jobs), "local backbones to", args.outdir)


if __name__ == "__main__":
    main()
