AI-assisted (Claude), exploration code for CS581 project

# Rep bank

Each `<REP>.tar.gz` unpacks to `<REP>/`, a gcmgen-layout rep dir holding only the inputs for a merge-only rerun (vote.py / gg.py):
`inputs/subalignments` (MAGUS's 25 subset alignments), `inputs/backbones` (MAGUS's L-INS-i backbones, `backbone_N_mafft.txt`),
`sets/s0` (unaligned backbone sequence sets), `true.fasta` (reference alignment used for FastSP), `unaligned.fasta`,
`magus.json` (MAGUS run score), and `true_tree.nwk` for simulated sets. No merged outputs (variants/, vote/, trees).
`aligned/` is left out because it repeats the backbones as `aligned/linsi/s0/backbone_N.fa`.
The columns are described in MANIFEST.tsv.
