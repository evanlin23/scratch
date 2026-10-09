"""Build oracle inputs that isolate where MAGUS's error comes from.

    python -m gcmx.oracle TRUE.fa SUBALIGNMENT_DIR BACKBONE_DIR OUTDIR

Writes, using the true alignment restricted to the *same* taxa MAGUS used:
  OUTDIR/true_subalignments/   true alignment induced on each subset
  OUTDIR/true_backbones/       true alignment induced on each backbone's taxa
  OUTDIR/true_full/true.fasta  the full true alignment (a single perfect backbone)

Mixing estimated/true subalignments with estimated/true backbones separates
(1) constraint (subset alignment) error, (2) backbone alignment error and
(3) error introduced by GCM's clustering + ordering itself.
"""

import os
import sys

from . import fasta


def main(true_path, subalignment_dir, backbone_dir, outdir):
    true = fasta.upper(fasta.read(true_path))
    for name, src in (("true_subalignments", subalignment_dir), ("true_backbones", backbone_dir)):
        os.makedirs(os.path.join(outdir, name), exist_ok=True)
        for filename in sorted(os.listdir(src)):
            taxa = list(fasta.read(os.path.join(src, filename)))
            fasta.write(fasta.restrict(true, taxa), os.path.join(outdir, name, filename))
    os.makedirs(os.path.join(outdir, "true_full"), exist_ok=True)
    fasta.write(true, os.path.join(outdir, "true_full", "true.fasta"))


if __name__ == "__main__":
    main(*sys.argv[1:5])
