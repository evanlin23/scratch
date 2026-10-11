AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R49, SIMHIGH_R50: vote-model tree replicates (helper p49)

Two replicates, SIMHIGH_R49 and SIMHIGH_R50 (AliSim, cs581/gcmtrees/code/gen.sh seeds), one MAGUS draw each (paper
flags) via `cs581/gcmtrees/code/run.sh SIMHIGH_Rk`, plus `gg.py run REP 'linsi#es4'`. Driver: `scripts/rep.sh`.

- `magus.jsonl`: MAGUS end-to-end rows (bbtool_bench).
- `aln.jsonl`: per replicate, gg.py rows (`variant` linsi, linsi#es4, recipe = wsoft0.03:linsi&fftns2#es4, linsi#es3,
  linsi&fftns2-op3; no `B` key) followed by vote.py rows from `cs581/gcmvote/code/run.py` (`B` = 10; magus, hard,
  hard-bb, soft4). FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the true tree). Methods `true`,
  `magus` (gg linsi = MAGUS's merge), `es4` (linsi#es4), `recipe`, `vote_hard`, `vote_hard-bb`. Run as 4 independent
  FastTree processes at a time.

Vote rep dirs: /opt/work/gcmvote/reps/SIMHIGH_Rk with `inputs`, `true.fasta`, `unaligned.fasta` symlinked from the
gcmtrees rep (separate results.jsonl because run.py expects a `B` key on every row).

Checks: on both replicates vote.py `magus` equals gg.py `linsi` (MAGUS's merge) after sorting rows by name
(`scripts/norm.py`); SPFN/SPFP/TC identical. On R49 vote `hard` and `hard-bb` kept the same edge set (79943 edges)
and give the same alignment up to row order (same nRF too); on R50 they differ (102319 vs 100165 kept edges).

Bank: cs581/gcmvote/bank/SIMHIGH_R49.tar.gz, SIMHIGH_R50.tar.gz (`scripts/bank.sh`), rows in bank/MANIFEST.tsv.
