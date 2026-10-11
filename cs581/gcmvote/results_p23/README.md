AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R23, SIMHIGH_R24: vote-model replicates (helper p23)

Two replicates (AliSim, cs581/gcmtrees/code/gen.sh seeds), one MAGUS draw each (paper flags) via
`cs581/gcmtrees/code/run.sh SIMHIGH_Rk`, plus `gg.py run REP 'linsi#es4'` (`scripts/aln_rep.sh`).

- `magus.jsonl`: MAGUS end-to-end rows (bbtool_bench).
- `aln.jsonl`: per replicate, gg.py rows (`variant` linsi, recipe = wsoft0.03:linsi&fftns2#es4, linsi#es3,
  linsi&fftns2-op3, linsi#es4; no `B` key) followed by vote.py rows from `cs581/gcmvote/code/run.py`
  (`B` = 10; magus, hard, hard-bb, soft4). FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the true tree). Methods `true`,
  `magus` (gg.py linsi = MAGUS's merge), `es4`, `recipe`, `vote_hard`, `vote_hard-bb`; 4 FastTree processes at a
  time (`scripts/trees_rep.sh`).

Vote rep dir: /opt/work/gcmvote/reps/SIMHIGH_Rk with `inputs`, `true.fasta`, `unaligned.fasta` symlinked from the
gcmtrees rep. Bank tarballs: `cs581/gcmvote/bank/` (`scripts/bank_rep.sh`).

Checks: vote.py `magus` equals gg.py `linsi` (MAGUS's merge) up to row order (`scripts/norm.py`; same SP scores).
