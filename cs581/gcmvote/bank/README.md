AI-assisted (Claude), exploration code for CS581 project

# Rep bank (helper p31)

One tarball per replicate (SIMHIGH_R31, SIMHIGH_R32), same layout as helper t20's bank: inputs for a merge-only
rerun (vote.py / gg.py) without rerunning MAGUS: `inputs/subalignments` (MAGUS's 25 subset alignments),
`inputs/backbones` (MAGUS's 10 L-INS-i backbones), `aligned/<tool>/s0/` (gcmgen backbone alignments: linsi,
fftns2, fftns2-op3), `sets/s0/`, `true.fasta`, `true_tree.nwk`, `unaligned.fasta`, `magus.json`. (No
subsets.json is written by this pipeline.) Merged outputs are left out. Columns are listed in MANIFEST.tsv.
