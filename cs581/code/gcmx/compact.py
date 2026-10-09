"""Majority-vote column compaction (post-processing for over-split alignments).

    python -m gcmx.compact IN.fasta OUT.fasta EVIDENCE [EVIDENCE ...] [--threshold 0.5] [--window 12]

EVIDENCE: alignment files or directories of them (full or partial alignments of
the same sequences; HMMER lower-case insertion letters carry no evidence).

Two columns c1 < c2 of the input alignment are merged when
  * no sequence has letters in both (compatible),
  * the merge keeps every sequence's letters in order: either no sequence of c2
    has a letter strictly between c1 and c2 (move c2 into c1) or no sequence of
    c1 does (move c1 into c2), and
  * support >= threshold, where support = fraction of (letter of c1, letter of
    c2, evidence alignment) triples -- over triples where that evidence contains
    both letters -- in which the evidence puts the two letters in one column.
Each column keeps a per-evidence "signature" (counts of the evidence columns
its letters fall in), so support is a sparse dot product. Left-to-right
passes, merging each column with its best valid partner within `window`
columns, repeat until nothing changes.
"""

import argparse
import collections
import os

import numpy as np

from . import fasta


def letter_columns(aln_path, names):
    """sequence -> array: residue index -> evidence column (-1 if absent or an insertion)."""
    aln = fasta.read(aln_path)
    out = {}
    for name in names:
        if name not in aln:
            continue
        cols, col = [], 0
        for ch in aln[name]:
            if ch == "-":
                col += 1
            elif ch == ".":
                continue
            elif ch.islower():
                cols.append(-1)
            else:
                cols.append(col)
                col += 1
        out[name] = np.array(cols, dtype=np.int64)
    return out


def compact(inp, out, evidence, threshold=0.5, window=12):
    aln = fasta.upper(fasta.read(inp))
    names = list(aln)
    raw = fasta.ungap(aln)
    grid = np.array([np.frombuffer(aln[n].encode(), dtype=np.uint8) for n in names])
    present = grid != ord("-")
    res_index = np.cumsum(present, axis=1) - 1
    L = grid.shape[1]

    files = []
    for e in evidence:
        files += [os.path.join(e, f) for f in sorted(os.listdir(e))] if os.path.isdir(e) else [e]
    evid = [letter_columns(f, names) for f in files]

    letters = [list(zip(np.nonzero(present[:, c])[0].tolist(), res_index[present[:, c], c].tolist())) for c in range(L)]
    seqs = [set(s for s, _ in col) for col in letters]
    sigs = []
    for col in letters:
        sig = []
        for ev in evid:
            counter = collections.Counter()
            for s, r in col:
                arr = ev.get(names[s])
                if arr is not None and arr[r] >= 0:
                    counter[int(arr[r])] += 1
            sig.append(counter)
        sigs.append(sig)

    def support(a, b):
        hit = tot = 0
        for ca, cb in zip(sigs[a], sigs[b]):
            na, nb = sum(ca.values()), sum(cb.values())
            if na and nb:
                tot += na * nb
                small, large = (ca, cb) if len(ca) < len(cb) else (cb, ca)
                hit += sum(v * large.get(k, 0) for k, v in small.items())
        return hit / tot if tot else 0.0

    alive = list(range(L))
    merged, changed = 0, True
    while changed:
        changed = False
        i = 0
        while i < len(alive):
            c1 = alive[i]
            between, best = set(), None
            for j in range(i + 1, min(len(alive), i + 1 + window)):
                c2 = alive[j]
                if not (seqs[c1] & seqs[c2]):
                    into_c1 = not (seqs[c2] & between)
                    if into_c1 or not (seqs[c1] & between):
                        sup = support(c1, c2)
                        if sup >= threshold and (best is None or sup > best[0]):
                            best = (sup, j, into_c1)
                between |= seqs[c2]
            if best:
                _, j, into_c1 = best
                c2 = alive[j]
                keep, drop = (c1, c2) if into_c1 else (c2, c1)
                letters[keep] += letters[drop]
                seqs[keep] |= seqs[drop]
                for k in range(len(evid)):
                    sigs[keep][k] = sigs[keep][k] + sigs[drop][k]
                alive.remove(drop)
                merged += 1
                changed = True
                if keep == c1:
                    continue  # try to grow the same column further
            i += 1

    rows = {n: bytearray(b"-" * len(alive)) for n in names}
    for j, c in enumerate(alive):
        for s, r in letters[c]:
            rows[names[s]][j] = ord(raw[names[s]][r])
    fasta.write({n: rows[n].decode() for n in names}, out)
    return merged, L, len(alive)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inp")
    parser.add_argument("out")
    parser.add_argument("evidence", nargs="+")
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--window", type=int, default=12)
    args = parser.parse_args()
    merged, before, after = compact(args.inp, args.out, args.evidence, args.threshold, args.window)
    print("merged {} column pairs: {} -> {} columns".format(merged, before, after))


if __name__ == "__main__":
    main()
