AI-assisted (Claude), exploration code for CS581 project

# Rep bank (helper p23)

One tarball per replicate, same layout as helper t20's bank: `inputs/subalignments` (MAGUS's 25 subset alignments),
`inputs/backbones` (MAGUS's 10 L-INS-i backbones), `aligned/<tool>/s0/` (gcmgen backbone alignments: linsi,
fftns2, fftns2-op3, on the same 10 sets), `sets/s0/` (unaligned backbone sequence sets), `true.fasta` (FastSP
reference), `true_tree.nwk` (AliSim true tree), `unaligned.fasta`, `magus.json` (MAGUS end-to-end score). Merged
outputs are left out. Extract, then point vote.py/gg.py at the extracted dir. Columns in MANIFEST.tsv.
