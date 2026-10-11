AI-assisted (Claude), exploration code for CS581 project

Helper h1 rows. Dataset: SIMHIGH_R5 only (R6-R8 reassigned to helpers t06/t07/t08 by the orchestrator).

- aln.jsonl: gcmgen gg.py variants on the gcmtrees rep (linsi = MAGUS merge, recipe, es3, hard, linsi#es4).
- vote_aln.jsonl: gcmvote vote.py/run.py variants (magus, es4, hard, soft, soft2, soft4, hard-bb, soft-bb,
  hard+mask), built in a separate rep dir that links the gcmtrees rep's inputs/ and true.fasta.
- trees.jsonl: FastTree -lg -gamma (protbench trees.py), nRF vs the true tree, run 4 at a time.
  Methods: true, magus, es4, recipe, es3, hard, vote_hard, vote_soft, vote_soft2, vote_soft4, vote_hard-bb,
  vote_soft-bb, vote_hard_mask_masked (= vote hard+mask, out.masked.fasta; 3 of 10028 columns masked).
  All 13 alignments differ byte-wise (cmp), so no tree was reused. vote_magus and vote_es4 were not treed;
  they are not byte-identical to the gcmtrees magus / es4 alignments.
- magus.jsonl: MAGUS draw 0 (paper flags) from /opt/work/gcmtrees.
