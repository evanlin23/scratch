AI-assisted (Claude), exploration code for CS581 project

# Rep bank (helper p47)

One tarball per replicate (SIMHIGH_R47, SIMHIGH_R48). Each holds the inputs for a merge-only rerun (vote.py / gg.py)
without rerunning MAGUS: `inputs/subalignments` (MAGUS's 25 subset alignments), `inputs/backbones` (MAGUS's 10
L-INS-i backbones), `aligned/<tool>/s0/` (gcmgen backbone alignments: linsi, fftns2, fftns2-op3, same 10 sets),
`sets/s0/` (unaligned backbone sets), `true.fasta` (FastSP reference), `true_tree.nwk` (AliSim true tree),
`unaligned.fasta`, and `magus.json` (MAGUS end-to-end score). No subsets.json is written by this pipeline.
Merged outputs are left out. Columns are listed in MANIFEST.tsv.
