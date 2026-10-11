AI-assisted (Claude), exploration code for CS581 project

SIMHIGH_R19.tar.gz (helper t19): gcmgen-layout rep dir from the gcmtrees MAGUS draw (paper flags), without merged
outputs. Contents under SIMHIGH_R19/: inputs/subalignments (25 MAGUS subsets), inputs/backbones (MAGUS's 10
L-INS-i backbones, 200 seqs each), aligned/{linsi,fftns2,fftns2-op3}/s0 (backbone alignments used by gg.py
variants), sets/s0 (unaligned backbone sets), true.fasta (AliSim reference, uppercased), true_tree.nwk (AliSim
tree), unaligned.fasta, magus.json (MAGUS run scores). Use: untar, then `python3 cs581/gcmvote/code/vote.py
SIMHIGH_R19 VARIANT` or `gg.py run SIMHIGH_R19 VARIANT`. No bb5/bb20 dirs were made for this rep.
