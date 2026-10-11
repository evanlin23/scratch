"""Build the inputs under $MLDATA (default /opt/data/mlcache).

    python prep.py

M1HF/R0-R4: the exact 1000M1-HF inputs of Park, Zaharias & Warnow 2021 (IDB-7008049,
  1000M1_HF_Analysis.tar.gz unpacked in /opt/data/park2021/1000M1_HF), rebuilt from their
  full-width per-subset alignments (as in cs581/ml/code/validate_park.py); true tree = ROSE rose.tt.
M1HF/R5-R9: ROSE 1000M1 R5-R9 (MAGUS paper Datasets.zip, IDB-2643961) fragmented with
  make_frag.py (same protocol, measured on the published 1000M1-HF).
RNASimHF/R0, RNASim10KHF/R0-R1: RNASim 1K / 10K (Datasets.zip) fragmented with make_frag.py.
"""
import glob
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_frag  # noqa: E402
import runtrees as rt  # noqa: E402

H = "/opt/data/park2021/1000M1_HF"
ROSE = "/opt/data/Datasets/ROSE/1000M1"


def link(src, dst):
    if not os.path.exists(dst):
        os.symlink(src, dst)


for r in range(5):
    d = os.path.join(rt.MLDATA, "M1HF", "R%d" % r)
    os.makedirs(d, exist_ok=True)
    if not os.path.exists(os.path.join(d, "true_align.fasta")):
        with open(os.path.join(d, "true_align.fasta"), "w") as out:
            for f in sorted(glob.glob(os.path.join(H, "CreateConstraintTrees/FastTree/500/R%d/output/sequence_partition_*.out" % r))):
                out.write(open(f).read().rstrip("\n") + "\n")
    link(os.path.join(ROSE, "R%d" % r, "rose.tt"), os.path.join(d, "true_tree.tre"))
    n, _ = rt.read_fasta(os.path.join(d, "true_align.fasta"))
    assert len(n) == len(set(n)) == 1000, (r, len(n))
for r in range(5, 10):
    make_frag.main("1000M1", r)
    src = os.path.join(rt.MLDATA, "1000M1HF", "R%d" % r)
    dst = os.path.join(rt.MLDATA, "M1HF", "R%d" % r)
    if not os.path.exists(dst):
        shutil.move(src, dst)
for r in range(1):
    make_frag.main("RNASim", r)
for r in range(2):
    make_frag.main("RNASim10K", r)

# RNASim10KHF: analysis alignment = columns with < 95% gaps among the full-length sequences
# (ungapped length >= 0.5 x median), as in the epang pilot (memory for 5,000-tip EPA-ng); all arms use it.
for r in range(2):
    d = os.path.join(rt.MLDATA, "RNASim10KHF", "R%d" % r)
    if os.path.exists(os.path.join(d, "true_mask.fasta")):
        continue
    import statistics
    n, s = rt.read_fasta(os.path.join(d, "true_align.fasta"))
    s = [x.upper().replace(".", "-") for x in s]
    L = [len(x.replace("-", "")) for x in s]
    med = statistics.median(L)
    full = [x for x, l in zip(s, L) if l >= 0.5 * med]
    keep = [j for j in range(len(s[0])) if sum(x[j] == "-" for x in full) < 0.95 * len(full)]
    rt.write_fasta(os.path.join(d, "true_mask.fasta"), n, ["".join(x[j] for j in keep) for x in s])
    print(d, "masked", len(s[0]), "->", len(keep))
