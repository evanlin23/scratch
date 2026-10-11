AI-assisted (Claude), exploration code for CS581 project

# Rep bank (helpers t20, p25)

One tarball per replicate. Each holds the inputs for a merge-only rerun (vote.py / gg.py) without rerunning MAGUS:
`inputs/subalignments` (MAGUS's 25 subset alignments), `inputs/backbones` (MAGUS's 10 L-INS-i backbones),
`aligned/<tool>/s0/` (gcmgen backbone alignments by tool: linsi, fftns2, fftns2-op3, on the same 10 sets),
`sets/s0/` (the unaligned backbone sequence sets), `true.fasta` (the reference used for FastSP scoring),
`true_tree.nwk` (the AliSim true tree), `unaligned.fasta`, and `magus.json` (MAGUS end-to-end score), `subsets.json` (p25 tarballs only). Merged
outputs (variants/, vote/, trees) are left out. Extract it, then point vote.py/gg.py at the extracted dir.
Columns are listed in MANIFEST.tsv.
