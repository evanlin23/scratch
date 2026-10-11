AI-assisted (Claude), exploration code for CS581 project

# Helper t19: SIMHIGH_R19 (vote model)

One MAGUS draw (paper flags) via `cs581/gcmtrees/code/run.sh SIMHIGH_R19`, plus `gg.py run REP 'linsi#es4'`.
Vote variants via `lane.sh` (run.py, B = 10) on a separate rep dir that links the gcmtrees rep's inputs
(run.py and gg.py both append to REP/results.jsonl in different row formats).

- `aln.jsonl`: gg.py rows (variant `linsi` = MAGUS merge, `linsi#es4`, recipe, `linsi#es3`, `linsi&fftns2-op3`)
  followed by run.py rows (vote variants, with `B`).
- `magus.jsonl`: bbtool_bench row of the draw. `trees.jsonl`: trees.py rows (FastTree -lg -gamma).
- Tree names: magus = gg `linsi`, es4 = gg `linsi#es4`, recipe = `wsoft0.03:linsi&fftns2#es4`, es3 = `linsi#es3`,
  hard = `linsi&fftns2-op3`; vote_X = vote.py X; vote_hard_mask_masked = `hard+mask` out.masked.fasta.
- Checks (sequence-content comparison; no pair was byte-identical because row order differs):
  vote `magus` has the same content as MAGUS's merge (gg `linsi`), so it reproduces it;
  vote `hard+mask` out.fasta equals vote `hard` (mask removed 2 of 7663 columns, so the masked file was treed);
  vote `es4` differs from gg `linsi#es4`. No vote alignment duplicates one already treed.
