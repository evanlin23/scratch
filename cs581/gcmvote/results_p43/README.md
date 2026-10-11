AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R43, SIMHIGH_R44: vote-model tree replicates (helper p43)

Same procedure as results_t20. Per replicate (AliSim, cs581/gcmtrees/code/gen.sh seeds): `cs581/gcmtrees/code/run.sh NAME`
(one MAGUS draw, paper flags, plus the gcmtrees merge variants), then `gg.py run REP 'linsi#es4'`.

- `magus.jsonl`: MAGUS end-to-end row (bbtool_bench).
- `aln.jsonl`: gg.py rows (`variant` linsi = MAGUS's merge, linsi#es4, recipe = wsoft0.03:linsi&fftns2#es4, linsi#es3,
  linsi&fftns2-op3; no `B` key) followed by vote.py rows from `cs581/gcmvote/code/run.py` (`B` = 10; magus, hard,
  hard-bb, soft4). FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the true tree). Methods `true`, `magus`,
  `es4`, `recipe` (gcmtrees rep alignments) and `vote_hard`, `vote_hard-bb`; 4 independent trees.py processes at a time.

Vote rep dir: /opt/work/gcmvote/reps/NAME with `inputs`, `true.fasta`, `unaligned.fasta` symlinked from the gcmtrees rep.
Bank: cs581/gcmvote/bank/NAME.tar.gz (merge-only inputs) + MANIFEST.tsv row.

Checks: vote.py `magus` equals gg.py `linsi` (MAGUS's merge): identical SPFN/SPFP/TC and identical row-sorted md5
(results_t20/scripts/norm.py). The six treed alignments have distinct row-sorted md5s.

`scripts/rep.sh` runs one replicate end to end; `scripts/collect.sh` rebuilds these files and the bank.
