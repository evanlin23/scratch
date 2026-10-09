"""Derived alignments for the pilots, written next to the cached inputs.

    python make_variants.py DATASET REP

  gcm_maskT.fasta     MAGUS alignment, columns whose homology support from the PASTA
                      alignment is < T/100 removed (T in 50, 70)
  gcm_oracleT.fasta   same, but support measured against the TRUE alignment (upper bound)
  gcm_pasta.fasta     MAGUS and PASTA alignments concatenated column-wise: every homology
                      is counted once per alignment that asserts it (a cheap soft weighting)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import colsupport as cs  # noqa: E402

MLDATA = os.environ.get("MLDATA", "/opt/data/mlcache")


def write_cols(path, aln, keep):
    with open(path, "w") as f:
        for n, seq in aln.items():
            f.write(">%s\n%s\n" % (n, np.frombuffer(seq.encode(), dtype=np.uint8)[keep].tobytes().decode()))


def main(ds, rep):
    d = os.path.join(MLDATA, ds, "R%s" % rep)
    gcm = cs.read_fasta(os.path.join(d, "gcm.fasta"))
    pasta = cs.read_fasta(os.path.join(d, "pasta_align.fasta"))
    true = cs.read_fasta(os.path.join(d, "true_align.fasta"))
    sp, _ = cs.support(gcm, pasta)
    st, _ = cs.support(gcm, true)
    ncol = len(sp)
    for t in (50, 70):
        for name, s in (("mask", sp), ("oracle", st)):
            keep = np.nonzero(s >= t / 100.0)[0]
            write_cols(os.path.join(d, "gcm_%s%d.fasta" % (name, t)), gcm, keep)
            print(ds, rep, "gcm_%s%d" % (name, t), "kept %d/%d" % (len(keep), ncol))
    with open(os.path.join(d, "gcm_pasta.fasta"), "w") as f:
        for n in gcm:
            f.write(">%s\n%s%s\n" % (n, gcm[n], pasta[n]))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
