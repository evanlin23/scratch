AI-assisted (Claude), exploration code for CS581 project

# Rep bank (helper p37)

One tarball per replicate (same layout as helper t20). Each holds the inputs for a merge-only rerun (vote.py /
gg.py) without rerunning MAGUS: `inputs/subalignments` (MAGUS's 25 subset alignments), `inputs/backbones`
(MAGUS's 10 L-INS-i backbones), `aligned/<tool>/s0/` (gcmgen backbone alignments by tool: linsi, fftns2,
fftns2-op3, on the same 10 sets), `sets/s0/` (the unaligned backbone sequence sets), `true.fasta` (FastSP
reference), `true_tree.nwk` (AliSim true tree), `unaligned.fasta`, `magus.json` (MAGUS end-to-end score) and
`subsets.json` (vote.py's subset membership, from vote/magus). Merged outputs are left out. Columns in MANIFEST.tsv.
