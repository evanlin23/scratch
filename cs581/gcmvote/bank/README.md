AI-assisted (Claude), exploration code for CS581 project

# Rep bank (helper p27)

One tarball per replicate, same layout as the other helpers' banks. Each holds the inputs for a merge-only rerun
(vote.py / gg.py) without rerunning MAGUS: `inputs/subalignments` (MAGUS's 25 subset alignments),
`inputs/backbones` (MAGUS's 10 L-INS-i backbones), `aligned/<tool>/s0/` (gcmgen backbone alignments: linsi,
fftns2, fftns2-op3 on the same 10 sets), `sets/s0/` (unaligned backbone sets), `true.fasta` (FastSP reference),
`true_tree.nwk` (AliSim true tree), `unaligned.fasta`, `magus.json`. Columns are listed in MANIFEST.tsv.
