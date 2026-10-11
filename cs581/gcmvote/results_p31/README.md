AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R31, SIMHIGH_R32: vote-model tree test (helper p31)

New replicates (AliSim, cs581/gcmtrees/code/gen.sh seeds), one MAGUS draw each (paper flags) via
`cs581/gcmtrees/code/run.sh SIMHIGH_Rk`, plus `gg.py run REP 'linsi#es4'`. Same procedure as results_t20.

- `magus.jsonl`: MAGUS end-to-end row (bbtool_bench).
- `aln.jsonl`: gg.py rows (`variant` linsi = MAGUS's merge, linsi#es4 = es4, recipe = wsoft0.03:linsi&fftns2#es4,
  linsi#es3, linsi&fftns2-op3; no `B` key) followed by vote.py rows from `cs581/gcmvote/code/run.py`
  (`B` = 10; magus, hard, hard-bb, soft4). FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the true tree). Methods `true`,
  `magus`, `es4`, `recipe` (gcmtrees variants) and `vote_hard`, `vote_hard-bb`.

Vote rep dir: /opt/work/gcmvote/reps/SIMHIGH_Rk with `inputs` and `true.fasta` symlinked from the gcmtrees rep.

Checks: vote `magus` out.fasta equals gg.py `linsi` (MAGUS's merge) after sorting rows by name (`scripts/norm.py`,
same md5; identical SPFN/SPFP/TC). vote `hard` and `hard-bb` alignments differ (cmp), each tree built separately.

`scripts/`: `rep.sh` (whole lane for one replicate; 6 trees as independent trees.py processes, 4 at a time via
xargs -P 4) and `fin.sh` (bank tarball + MANIFEST row + jsonl collection).
