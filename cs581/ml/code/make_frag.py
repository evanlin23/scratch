"""Create a "half fragmentary" (HF) version of a simulated dataset, following the protocol of
Smirnov & Warnow 2021 / Park et al. 2021 as measured on the published 1000M1-HF data: 50% of the
sequences (chosen at random) are replaced by one contiguous fragment each, of length
~ Normal(0.25 * median ungapped length, 60) (clipped to >= 20), at a uniformly random start.
Other sequences and the true alignment columns are untouched (fragment residues keep their
true-alignment columns, the rest become gaps).

    python make_frag.py CONDITION REP [--amplicon]
CONDITION: ROSE name (1000M1..1000M4, 1000L1..), RNASim, or 16S.M (biological: CRW 16S.M reference
alignment, 901 sequences, reference tree = CRW tree with only well-supported edges, 47% resolved).
--amplicon: instead of a uniformly random start, every fragment starts at the same relative position
(45% of the sequence length, +- N(0, 15) residues), mimicking amplicon reads from one primer pair.
Writes $MLDATA/<CONDITION>HF[amp]/R<rep>/{true_align.fasta,true_tree.tre,unaligned.fasta};
seed = hash of (condition, rep, protocol). For 16S.M, REP only changes which sequences are fragmented.
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
    if cond == "16S.M":
        d = os.path.join(DATA, "Gutell", "16S.M", "R0")
        return os.path.join(d, "cleaned.alignment.fasta"), os.path.join(d, "16S.M.reference.nwk")
    if cond == "RNASim":
        d = os.path.join(DATA, "RNASim", "1000", "R%s" % rep)
        return os.path.join(d, "true_align.txt"), os.path.join(d, "true_tree.tre")
    d = os.path.join(DATA, "ROSE", cond, "R%s" % rep)
    return os.path.join(d, "rose.aln.true.fasta"), os.path.join(d, "rose.tt")


def main(cond, rep, amplicon=False):
    aln, tree = sources(cond, rep)
    out = os.path.join(rt.MLDATA, cond.replace(".", "") + ("HFamp" if amplicon else "HF"), "R%s" % rep)
    os.makedirs(out, exist_ok=True)
    rng = random.Random(zlib.crc32(("%s/%s%s" % (cond, rep, "/amp" if amplicon else "")).encode()))
    names, seqs = rt.read_fasta(aln)
    seqs = [s.upper().replace(".", "-") for s in seqs]
    med = statistics.median(len(s.replace("-", "")) for s in seqs)
    frag = set(rng.sample(range(len(names)), len(names) // 2))
    new = []
    for i, s in enumerate(seqs):
        if i not in frag:
            new.append(s)
            continue
        pos = [j for j, c in enumerate(s) if c != "-"]
        L = max(20, min(len(pos), int(round(rng.gauss(0.25 * med, 60)))))
        if amplicon:
            a = max(0, min(len(pos) - L, int(round(0.45 * len(pos) + rng.gauss(0, 15)))))
        else:
            a = rng.randint(0, len(pos) - L)
        keep = set(pos[a:a + L])
        new.append("".join(c if j in keep else "-" for j, c in enumerate(s)))
    rt.write_fasta(os.path.join(out, "true_align.fasta"), names, new)
    rt.write_fasta(os.path.join(out, "unaligned.fasta"), names, [x.replace("-", "") for x in new])
    if not os.path.exists(os.path.join(out, "true_tree.tre")):
        os.symlink(tree, os.path.join(out, "true_tree.tre"))
    print(out, "median", med, "fragments", len(frag))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], "--amplicon" in sys.argv)
