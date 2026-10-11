AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R21, SIMHIGH_R22: vote-model replicates (helper p21)

Same procedure as results_t20. Per replicate (AliSim via cs581/gcmtrees/code/gen.sh seeds), one MAGUS draw
(paper flags) via `cs581/gcmtrees/code/run.sh SIMHIGH_Rk`, then `gg.py run REP 'linsi#es4'`.

- `magus.jsonl`: MAGUS end-to-end row (bbtool_bench).
- `aln.jsonl`: gg.py rows (`variant` linsi = MAGUS's merge, linsi#es4, recipe = wsoft0.03:linsi&fftns2#es4,
  linsi#es3, linsi&fftns2-op3; no `B` key) followed by vote.py rows from `cs581/gcmvote/code/run.py`
  (`B` = 10; magus, hard, hard-bb, soft4). FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the true tree). Methods `true`,
  `magus`, `es4`, `recipe`, `vote_hard`, `vote_hard-bb`; 4 FastTree processes at a time, each its own process.

Vote rep dir: /opt/work/gcmvote/reps/SIMHIGH_Rk with `inputs` and `true.fasta` symlinked from the gcmtrees rep.

`scripts/`: rep.sh (whole per-replicate pipeline), bank.sh (bank tarball + MANIFEST row), collect_p21.sh,
norm.py (row-order-independent md5 of an alignment).

## Checks

- SIMHIGH_R21: vote `magus` equals gg.py `linsi` (MAGUS's merge) up to row order (norm.py 01c15d3630; same
  SP error). vote `hard` and `hard-bb` keep the same 101494 edges and give the same alignment up to row order
  (norm.py 6d8454b9be; files differ byte-wise only in row order), so their two trees come from the same
  alignment; both trees were built separately and both gave nRF 0.1424.
