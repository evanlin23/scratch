AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R4: vote-model helper run (t04)

One MAGUS draw (draw 0, paper flags; `cs581/gcmtrees/code/run.sh SIMHIGH_R4`) on a machine with 4 cores and 15 GB RAM.

- `magus.jsonl`: bbtool_bench row for the MAGUS run.
- `aln.jsonl`: gg.py rows (`variant` = linsi, `wsoft0.03:linsi&fftns2#es4` (recipe), linsi#es3,
  linsi&fftns2-op3 (hard), linsi#es4; rep dir /opt/work/gcmtrees/reps/SIMHIGH_R4), followed by
  gcmvote/code/run.py rows (with a `B` key; variants magus, es4, hard, soft, soft2, soft4, hard-bb, soft-bb, hard+mask;
  separate rep dir /opt/work/gcmvote/reps/SIMHIGH_R4 whose inputs/ is linked to the gcmtrees rep).
- `trees.jsonl`: protbench/code/trees.py rows (FastTree -lg -gamma, nRF vs /opt/data/sim/SIMHIGH/R4/tree.nwk).
  Methods: true, magus, es4, recipe, es3, hard (gg.py alignments), vote_hard, vote_soft, vote_soft2, vote_soft4, vote_hard-bb,
  vote_soft-bb, vote_hard_mask (= vote/hard+mask/out.masked.fasta).
- `raw/`: model.json (mixture fit, cutoff per exposure n) for each vote variant.

Notes:
- vote.py `magus` has the same rows as MAGUS's merge (variants/linsi/out.fasta), with a different row order, and the same
  FastSP scores (SPFN 0.286918, SPFP 0.183723).
- None of the vote alignments is byte-identical (md5) to another alignment that was treed, so no tree was reused.
  vote `es4` is not row-identical to gg.py `linsi#es4` (SPFN 0.16151 vs 0.16162).
- hard+mask removed 1 of 9539 columns.
- The tree `true` is built on /opt/data/sim/SIMHIGH/R4/sim.fa; it has the same residues per sequence as the rep's true.fasta.
