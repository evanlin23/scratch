# IQ-TREE check of MAGUS merge-step evidence filtering (helper a)

AI-assisted (Claude), exploration code for CS581 project.

- Data: SIMHIGH_R1, SIMHIGH_R2 from the gcmtrees rep bank (branch claude/cs581-gcmtrees, sha256 as in its MANIFEST.tsv).
- Alignments rebuilt with `cs581/gcmvote/code/run.py` (vote.py) for `magus`, `es4`, `hard-bb`.
  `magus` reproduces the bank's recorded MAGUS SP error (avgErr R1 0.25321, R2 0.22914).
- Trees: IQ-TREE 3.1.4 (bioconda, via cs581/code/setup.sh), `-m LG+G4 --fast -T 2 -seed 1`, 2 runs at a time
  on 4 cores (`code/iq.py`). nRF vs the true tree as in cs581/protbench/code/trees.py.
- Rows: `results_a/iqtree.jsonl` (rep, method, nRF, FN, FP, wall seconds, IQ-TREE version).
