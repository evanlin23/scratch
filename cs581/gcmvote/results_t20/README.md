AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R20: vote-model pilot (helper t20)

One replicate, SIMHIGH_R20 (AliSim, cs581/gcmtrees/code/gen.sh seeds), one MAGUS draw (paper flags) via
`cs581/gcmtrees/code/run.sh SIMHIGH_R20`, plus `gg.py run REP 'linsi#es4'`.

- `magus.jsonl`: MAGUS end-to-end row (bbtool_bench).
- `aln.jsonl`: gg.py rows (`variant` linsi, linsi#es4, recipe = wsoft0.03:linsi&fftns2#es4, linsi#es3,
  hard = linsi&fftns2-op3; no `B` key) followed by vote.py rows from `cs581/gcmvote/code/run.py`
  (`B` = 10; magus, es4, hard, soft, soft2, soft4, hard-bb, soft-bb, hard+mask). FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the true tree). Methods `true`,
  `magus`, `es4`, `recipe`, `es3`, `hard` (gcmtrees set) and `vote_<variant>`; `vote_hard+mask` is the
  masked alignment (out.masked.fasta).

Vote rep dir: /opt/work/gcmvote/reps/SIMHIGH_R20 with `inputs` and `true.fasta` symlinked from the gcmtrees rep
(separate results.jsonl because run.py expects a `B` key on every row).

Checks: vote.py `magus` equals gg.py `linsi` (MAGUS's merge) up to row order (same rows after sorting by name,
`scripts/norm.py`; SPFN/SPFP/TC identical). No vote alignment is byte-identical (cmp/md5) to another treed
alignment, so every listed tree was built separately. vote `hard+mask` out.fasta equals vote `hard` (row-sorted);
its masked alignment drops 9 columns.

`scripts/`: the tree launcher (4 FastTree processes at a time) used here.
