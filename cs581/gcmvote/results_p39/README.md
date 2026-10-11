AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R39, SIMHIGH_R40: vote-model tree replicates (helper p39)

Each replicate: SIMHIGH_Rk (AliSim, cs581/gcmtrees/code/gen.sh seeds), one MAGUS draw (paper flags) via
`cs581/gcmtrees/code/run.sh SIMHIGH_Rk`, plus `gg.py run REP 'linsi#es4'`; then
`cs581/gcmvote/code/run.py VOTEREP magus hard hard-bb soft4` (B = 10). Pipeline: `scripts/rep.sh`.

- `magus.jsonl`: MAGUS end-to-end row (bbtool_bench).
- `aln.jsonl`: gg.py rows (`variant` linsi = MAGUS's merge, linsi#es4, recipe = wsoft0.03:linsi&fftns2#es4,
  linsi#es3, linsi&fftns2-op3; no `B` key) followed by vote.py rows (`B` = 10; magus, hard, hard-bb, soft4).
  FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the true tree). Methods `true`,
  `magus`, `es4`, `recipe`, `vote_hard`, `vote_hard-bb`, 4 FastTree processes at a time.

Vote rep dir: /opt/work/gcmvote/reps/SIMHIGH_Rk with `inputs`, `true.fasta`, `unaligned.fasta` symlinked from the
gcmtrees rep (separate results.jsonl because run.py expects a `B` key on every row).

Bank: cs581/gcmvote/bank/SIMHIGH_Rk.tar.gz + MANIFEST.tsv (`scripts/bank.sh`). The gcmtrees reps have no
subsets.json, so the tarball holds the same files as helper t20's.

Checks:
- SIMHIGH_R39: vote `magus` equals gg.py `linsi` (MAGUS's merge) up to row order (row-sorted md5 equal,
  results_t20 `norm.py`; SPFN/SPFP/TC identical). vote `hard` and `hard-bb` kept the same 88545 edges and give
  the same alignment up to row order; both trees were built and give the same nRF.

Note: gcmtrees run.sh calls gcmtrees collect.sh, which overwrites cs581/gcmtrees/results/*.jsonl with only this
machine's rows; those files were restored from git and are not part of this branch's changes.
