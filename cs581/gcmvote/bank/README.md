AI-assisted (Claude), exploration code for CS581 project

Rep bank (helper h5d). Each `<REP>.tar.gz` unpacks to `<REP>/`, a gcmgen-layout replicate without merged outputs:
`inputs/subalignments` (MAGUS's 25 subset alignments), `inputs/backbones` (MAGUS's 10 L-INS-i backbones),
`sets/s0` (the backbones' unaligned sequence sets), `true.fasta` (reference used for scoring), `unaligned.fasta`,
`magus.json` (fresh MAGUS draw stats), `true_tree.tre` (simulated true tree, when present).
Use: `tar xzf REP.tar.gz -C /opt/work/gcmvote/reps && python3 cs581/gcmvote/code/vote.py /opt/work/gcmvote/reps/REP VARIANT`.
Not included: aligned/linsi (gg.py's cache, a copy of inputs/backbones), variants/, vote/.
