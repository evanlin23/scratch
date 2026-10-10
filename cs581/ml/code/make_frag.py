"""Create a "half fragmentary" (HF) version of a simulated dataset, following the protocol of
Smirnov & Warnow 2021 / Park et al. 2021 as measured on the published 1000M1-HF data: 50% of the
sequences (chosen at random) are replaced by one contiguous fragment each, of length
~ Normal(0.25 * median ungapped length, 60) (clipped to >= 20), at a uniformly random start.
Other sequences and the true alignment columns are untouched (fragment residues keep their
true-alignment columns, the rest become gaps).

    python make_frag.py CONDITION REP     # CONDITION: ROSE name (1000M1..1000M4, 1000L1..) or RNASim
writes $MLDATA/<CONDITION>HF/R<rep>/{true_align.fasta,true_tree.tre}; seed = hash of (condition, rep).
"""
import os
import random
import statistics
import sys
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runtrees as rt  # noqa: E402

DATA = "/opt/data/Datasets"


def sources(cond, rep):
    if cond == "RNASim":
        d = os.path.join(DATA, "RNASim", "1000", "R%s" % rep)
        return os.path.join(d, "true_align.txt"), os.path.join(d, "true_tree.tre")
    d = os.path.join(DATA, "ROSE", cond, "R%s" % rep)
    return os.path.join(d, "rose.aln.true.fasta"), os.path.join(d, "rose.tt")


def main(cond, rep):
    aln, tree = sources(cond, rep)
    out = os.path.join(rt.MLDATA, cond + "HF", "R%s" % rep)
    os.makedirs(out, exist_ok=True)
    rng = random.Random(zlib.crc32(("%s/%s" % (cond, rep)).encode()))
    names, seqs = rt.read_fasta(aln)
    seqs = [s.upper() for s in seqs]
    med = statistics.median(len(s.replace("-", "")) for s in seqs)
    frag = set(rng.sample(range(len(names)), len(names) // 2))
    new = []
    for i, s in enumerate(seqs):
        if i not in frag:
            new.append(s)
            continue
        pos = [j for j, c in enumerate(s) if c != "-"]
        L = max(20, min(len(pos), int(round(rng.gauss(0.25 * med, 60)))))
        a = rng.randint(0, len(pos) - L)
        keep = set(pos[a:a + L])
        new.append("".join(c if j in keep else "-" for j, c in enumerate(s)))
    rt.write_fasta(os.path.join(out, "true_align.fasta"), names, new)
    if not os.path.exists(os.path.join(out, "true_tree.tre")):
        os.symlink(tree, os.path.join(out, "true_tree.tre"))
    print(out, "median", med, "fragments", len(frag))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
