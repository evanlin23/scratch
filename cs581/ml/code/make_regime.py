"""Data regimes for the headroom survey (model-agnostic data transformations of ROSE replicates).

    python make_regime.py COND REP            -> $MLDATA/<COND>/R<rep>  (full-length true alignment)
    python make_regime.py COND REP --sites K  -> $MLDATA/<COND>s<K>/R<rep> (K random non-empty columns
                                                  of the true alignment: "many taxa, few sites")
True tree = ROSE rose.tt. Column sample seeded by (COND, REP, K).
"""
import os
import random
import sys
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import runtrees as rt  # noqa: E402


def main(cond, rep, sites=None):
    src = os.path.join("/opt/data/Datasets/ROSE", cond, "R%s" % rep)
    out = os.path.join(rt.MLDATA, cond + ("s%d" % sites if sites else ""), "R%s" % rep)
    os.makedirs(out, exist_ok=True)
    if not os.path.exists(os.path.join(out, "true_tree.tre")):
        os.symlink(os.path.join(src, "rose.tt"), os.path.join(out, "true_tree.tre"))
    if os.path.exists(os.path.join(out, "true_align.fasta")):
        return
    names, seqs = rt.read_fasta(os.path.join(src, "rose.aln.true.fasta"))
    seqs = [s.upper() for s in seqs]
    if sites:
        cols = [i for i in range(len(seqs[0])) if any(s[i] != "-" for s in seqs)]
        keep = sorted(random.Random(zlib.crc32(("%s/%s/%d" % (cond, rep, sites)).encode())).sample(cols, sites))
        seqs = ["".join(s[i] for i in keep) for s in seqs]
    rt.write_fasta(os.path.join(out, "true_align.fasta"), names, seqs)
    print(out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[4]) if "--sites" in sys.argv else None)
