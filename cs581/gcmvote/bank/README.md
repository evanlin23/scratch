AI-assisted (Claude), exploration code for CS581 project

# Rep bank (helper p21)

One tarball per replicate (SIMHIGH_R21, SIMHIGH_R22), same layout as helper t20's bank. Each holds the inputs
for a merge-only rerun (vote.py / gg.py) without rerunning MAGUS: `inputs/subalignments` (MAGUS's 25 subset
alignments), `inputs/backbones` (MAGUS's 10 L-INS-i backbones), `aligned/<tool>/s0/` (gcmgen backbone alignments
by tool: linsi, fftns2, fftns2-op3, on the same 10 sets), `sets/s0/` (the unaligned backbone sequence sets),
`true.fasta` (the FastSP reference), `true_tree.nwk` (the AliSim true tree), `unaligned.fasta`, `magus.json`
(MAGUS end-to-end score) and `subsets.json` (vote.py's subset map, from vote/magus). Merged outputs are left out.
Columns are listed in MANIFEST.tsv.
