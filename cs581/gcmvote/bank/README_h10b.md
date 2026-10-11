AI-assisted (Claude), exploration code for CS581 project

1000M1_R2.tar.gz (helper h10b): gcmgen-layout rep dir for ROSE 1000M1 R2, one fresh MAGUS draw (paper flags).
inputs/subalignments (25 subsets), inputs/backbones (MAGUS's 10 L-INS-i backbones, 200 seqs each),
sets/s0 (unaligned backbone sets) and aligned/linsi/s0 (gg's copy of the aligned backbones), true.fasta (reference),
true.tree (ROSE rose.tt), unaligned.fasta, magus.json (MAGUS scores). No merged outputs. No bb5/bb20 dirs were made.
Merge-only redo: python3 cs581/gcmvote/code/run.py <extracted>/1000M1_R2 hard
