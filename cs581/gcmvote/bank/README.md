AI-assisted (Claude), exploration code for CS581 project

Rep bank (helper h10c). Each tarball is a gcmgen-layout rep dir with only the inputs needed for a merge-only
rerun (vote.py / run.py / gg.py): inputs/subalignments (MAGUS's 25 subset alignments), inputs/backbones
(MAGUS's 10 L-INS-i backbones), aligned/linsi + sets (gg's backbone alignments and backbone sequence sets),
true.fasta (reference used for scoring), true.tree (ROSE true tree, rose.tt), unaligned.fasta, magus.json
(MAGUS draw stats). Merged outputs (variants/, vote/, trees) are excluded. vote.py writes subsets.json itself.
Paths inside the tar are listed in MANIFEST.tsv.
