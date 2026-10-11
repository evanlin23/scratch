AI-assisted (Claude), exploration code for CS581 project

# Rep bank (helper t10)

One tarball per replicate with what a merge-only run (gcmvote vote.py / gcmgen gg.py) needs, so that
new voting models can be tested without re-running MAGUS. Untar into a reps dir and point vote.py/run.py at it.

SIMHIGH_R10.tar.gz -> SIMHIGH_R10/
- inputs/subalignments/   MAGUS's 25 subset alignments (subalignment_subset_N.txt)
- inputs/backbones/       MAGUS's 10 L-INS-i backbones x 200 seqs (backbone_N_mafft.txt)
- aligned/{linsi,fftns2,fftns2-op3}/s0/  gg.py's backbone alignments of the same 10 sets (used by the recipe /
                          hard variants)
- sets/s0/                the 10 unaligned backbone sequence sets
- true.fasta              true (AliSim) alignment, used for FastSP scoring
- true_tree.nwk           true tree (/opt/data/sim/SIMHIGH/R10/tree.nwk)
- unaligned.fasta         the 1000 unaligned sequences
- magus.json              MAGUS draw-0 score row
Excluded: variants/, vote/, trees. No bb5/bb20 directories were made for this replicate.
Restore check: untarred, vote.py magus reproduces the original vote/magus alignment row-for-row.
