"""Column support of alignment A from one or more other alignments of the same sequences.

For column c of A with residue set R_c, its support from alignment B is the fraction of the
homology pairs asserted by c (pairs of residues in c) that B also asserts:
    sum_g C(|g|, 2) / C(|R_c|, 2),   g ranging over the B-columns that receive R_c's residues.
Singleton columns (|R_c| < 2) get support 1. Computed in O(total residues) per pair (A, B).

    python colsupport.py A.fasta B.fasta [C.fasta ...] > A.support   # one float per A column
    python colsupport.py --mask T A.fasta B.fasta ... > A.masked.fasta  # keep columns with support >= T
"""
import sys
from collections import Counter

import numpy as np


def read_fasta(path):
    names, seqs, cur = [], [], []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if names:
                seqs.append("".join(cur))
            names.append(line[1:].split()[0])
            cur = []
        elif line:
            cur.append(line)
    seqs.append("".join(cur))
    return dict(zip(names, (s.upper().replace(".", "-") for s in seqs)))


def residue_columns(aln, names):
    """For each sequence, the alignment column of its k-th residue (numpy int arrays)."""
    out = {}
    for n in names:
        a = np.frombuffer(aln[n].encode(), dtype=np.uint8)
        out[n] = np.nonzero(a != ord("-"))[0]
    return out


def support(a, b):
    names = sorted(a)
    assert set(names) == set(b), "alignments must have the same sequences"
    ca, cb = residue_columns(a, names), residue_columns(b, names)
    ncol = len(next(iter(a.values())))
    # stack (colA, colB) for every residue
    colA = np.concatenate([ca[n] for n in names])
    colB = np.concatenate([cb[n] for n in names])
    assert len(colA) == len(colB), "sequences differ after removing gaps"
    sizeA = np.bincount(colA, minlength=ncol).astype(float)
    pairs = Counter(zip(colA.tolist(), colB.tolist()))
    agree = np.zeros(ncol)
    for (x, _), g in pairs.items():
        agree[x] += g * (g - 1) / 2
    total = sizeA * (sizeA - 1) / 2
    with np.errstate(invalid="ignore", divide="ignore"):
        s = np.where(total > 0, agree / np.maximum(total, 1), 1.0)
    return s, sizeA


def mean_support(a, others):
    return np.mean([support(a, b)[0] for b in others], axis=0)


if __name__ == "__main__":
    args = sys.argv[1:]
    mask = None
    if args[0] == "--mask":
        mask, args = float(args[1]), args[2:]
    a = read_fasta(args[0])
    s = mean_support(a, [read_fasta(p) for p in args[1:]])
    if mask is None:
        for v in s:
            print("%.4f" % v)
    else:
        keep = np.nonzero(s >= mask)[0]
        for n, seq in a.items():
            arr = np.frombuffer(seq.encode(), dtype=np.uint8)[keep]
            print(">%s\n%s" % (n, arr.tobytes().decode()))
        print("kept %d of %d columns" % (len(keep), len(s)), file=sys.stderr)
