AI-assisted (Claude), exploration code for CS581 project

SIMHIGH_R8, helper t08 (took over from gcmvote-h1). One MAGUS draw (draw 0, paper flags) via cs581/gcmtrees/code/run.sh.

- `magus.jsonl`: the MAGUS run row (bbtool_bench format).
- `aln.jsonl`: gg.py rows (`linsi`, `wsoft0.03:linsi&fftns2#es4` = recipe, `linsi#es3`, `linsi&fftns2-op3` = hard,
  `linsi#es4`; rep /opt/work/gcmtrees/reps/SIMHIGH_R8), then vote run.py rows (have a `B` field; `magus`, `es4`,
  `hard`, `soft`, `soft2`, `soft4`, `hard-bb`, `soft-bb`, `hard+mask`), run on a separate rep dir
  /opt/work/gcmvote/reps/SIMHIGH_R8 that symlinks the gcmtrees rep's inputs/ and true.fasta (run.py and gg.py
  cannot share one results.jsonl: run.py reads a `B` key gg.py rows lack).
- vote `magus` reproduces MAGUS's merge (same column content as gg `linsi`; only the sequence order in the file differs).
- vote `hard`, `hard-bb` and the unmasked `hard+mask` outputs are byte-identical; `hard+mask` masked 4 of 9915 columns.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, `RF` = nRF vs the true tree). Methods: true, magus,
  es4 (gg `linsi#es4`), recipe, es3, hard (gcmtrees trees.sh set) and vote_* (vote_hard+mask = out.masked.fasta).
  A row with `reused_from` is a vote alignment byte-identical to one already treed; that tree was reused.
