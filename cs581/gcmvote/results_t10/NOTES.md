AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R10 (helper t10)

- One MAGUS draw (draw 0, paper flags, 4 threads) via cs581/gcmtrees/code/run.sh SIMHIGH_R10; then gg.py `linsi#es4`.
- `aln.jsonl`: gg.py rows (variant `linsi`, `wsoft0.03:linsi&fftns2#es4` = recipe, `linsi#es3`,
  `linsi&fftns2-op3` = hard, `linsi#es4`) followed by gcmvote run.py rows (field `B`; variants magus, es4, hard,
  soft, soft2, soft4, hard-bb, soft-bb, hard+mask), all FastSP vs the true alignment.
- Vote variants were run in a separate rep dir (/opt/work/gcmvote/reps/SIMHIGH_R10) whose inputs/ and true.fasta
  are symlinks to the gcmtrees rep, so the two results.jsonl files stay separate.
- vote.py `magus` reproduces MAGUS's merge: same rows as gg.py `linsi` (MAGUS's merge re-run on the same subsets
  and backbones) up to sequence order, identical SPFN/SPFP/TC. Both differ from MAGUS's own output in SPFP by
  0.00014 (0.13527 vs 0.13540), as gg.py `linsi` / bbtool_bench merge-mafft do for this draw.
- vote.py `es4` and gg.py `linsi#es4` are not identical alignments (SP error 13.14 vs 13.16 %).
- `hard+mask` out.fasta = `hard` (same rows); the mask removed 2 of 8668 columns; its tree uses out.masked.fasta.
- `magus.jsonl`: bbtool_bench row for the MAGUS draw.
- `trees.jsonl`: FastTree -lg -gamma via protbench/code/trees.py, RF vs /opt/data/sim/SIMHIGH/R10/tree.nwk,
  4 trees at a time. true / magus (gg `linsi`) / recipe / es3 / hard (gg `linsi&fftns2-op3`) / es4 (gg `linsi#es4`)
  and vote_* = vote.py outputs (vote_hard+mask = hard+mask/out.masked.fasta). No vote alignment was byte-identical
  (cmp) to an already-treed one, so no tree was reused.
