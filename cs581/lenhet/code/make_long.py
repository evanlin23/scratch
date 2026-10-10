"""Build controlled length-heterogeneity datasets from a ROSE true alignment.

A fraction `frac` of sequences ("long" sequences) get extra sequence attached,
split at random between the 5' and 3' ends, of total length `mult` x their own
length. Two flank types:

  rand : i.i.d. letters with the dataset's base composition (non-homologous to
         everything; in the reference every flank letter is its own column).
  dom  : a second "domain": the flank is a sequence from a different ROSE
         replicate (a different family, `--dom-src`), so flanks of different
         long sequences ARE homologous to each other. The reference has the
         family-B columns (true alignment of B) after/before the A block.
         The whole B sequence is appended at the 3' end (same side for all,
         so all B pieces share one column block); mult is ignored.

Outputs (in OUT): unaligned.fasta, true.fasta (reference), long.txt (names).

    python make_long.py ALN OUT --frac 0.1 --mult 1.0 --type rand --seed 1
"""
import argparse
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "code"))
from gcmx import fasta  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("aln")
    ap.add_argument("out")
    ap.add_argument("--frac", type=float, default=0.1)
    ap.add_argument("--mult", type=float, default=1.0)
    ap.add_argument("--type", choices=["rand", "dom"], default="rand")
    ap.add_argument("--dom-src", help="true alignment of family B (dom type)")
    ap.add_argument("--n", type=int, default=0, help="subsample to n sequences (0 = all)")
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    aln = fasta.read(a.aln)
    names = list(aln)
    if a.n:
        names = sorted(rng.sample(names, a.n), key=names.index)
        aln = fasta.restrict(aln, names)
    nlong = round(a.frac * len(names))
    long_ = set(rng.sample(names, nlong))
    ncol = len(next(iter(aln.values())))
    letters = "".join(s.replace("-", "") for s in aln.values())
    comp = {c: letters.count(c) for c in set(letters)}
    alpha, weights = zip(*sorted(comp.items()))

    if a.type == "rand":
        # columns: [5' flank block][A block][3' flank block]; each flank letter
        # gets its own column (non-homologous to everything)
        left, right, raw = {}, {}, {}
        for n in names:
            core = aln[n].replace("-", "")
            if n in long_:
                tot = int(round(a.mult * len(core)))
                k = rng.randint(0, tot)
                left[n] = "".join(rng.choices(alpha, weights, k=k))
                right[n] = "".join(rng.choices(alpha, weights, k=tot - k))
            else:
                left[n] = right[n] = ""
            raw[n] = left[n] + core + right[n]
        L = sum(len(v) for v in left.values())
        R = sum(len(v) for v in right.values())
        ref, lo, ro = {}, 0, 0
        for n in names:
            ref[n] = ("-" * lo + left[n] + "-" * (L - lo - len(left[n])) + aln[n]
                      + "-" * ro + right[n] + "-" * (R - ro - len(right[n])))
            lo += len(left[n])
            ro += len(right[n])
    else:
        b = fasta.read(a.dom_src)
        bnames = rng.sample(list(b), nlong)
        bmap = dict(zip(sorted(long_, key=names.index), bnames))
        bb = fasta.restrict(b, bnames)
        bcol = len(next(iter(bb.values())))
        raw, ref = {}, {}
        for n in names:
            core = aln[n].replace("-", "")
            if n in long_:
                bs = bb[bmap[n]]
                braw = bs.replace("-", "")
                raw[n] = core + braw
                ref[n] = aln[n] + bs
            else:
                raw[n] = core
                ref[n] = aln[n] + "-" * bcol
        ref = fasta.restrict(ref, names)
    os.makedirs(a.out, exist_ok=True)
    fasta.write(raw, os.path.join(a.out, "unaligned.fasta"))
    fasta.write(ref, os.path.join(a.out, "true.fasta"))
    with open(os.path.join(a.out, "long.txt"), "w") as f:
        f.write("\n".join(n for n in names if n in long_) + "\n")


if __name__ == "__main__":
    main()
