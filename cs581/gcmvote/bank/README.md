AI-assisted (Claude), exploration code for CS581 project

# Rep bank (SIMHIGH_R3, helper t03)

Inputs to redo merge-only runs without re-running MAGUS. Each tarball unpacks to <REPNAME>/ in the gcmgen layout:
inputs/{subalignments,backbones} (MAGUS's 25 subsets, 10 L-INS-i backbones), aligned/<tool>/s0 (the extra backbone
alignments gg.py built: fftns2, fftns2-op3, linsi), sets/s0 (backbone sequence sets), true.fasta (reference),
unaligned.fasta, true_tree.nwk (AliSim tree), magus.json (MAGUS draw scores). See MANIFEST.tsv.
