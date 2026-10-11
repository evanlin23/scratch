AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R35, SIMHIGH_R36: vote-model tree power (helper p35)

New replicates (AliSim, cs581/gcmtrees/code/gen.sh seeds), one MAGUS draw each (paper flags) via
`cs581/gcmtrees/code/run.sh SIMHIGH_Rk`, then `gg.py run REP 'linsi#es4'`. Same formats as results_t20.

- `magus.jsonl`: MAGUS end-to-end row (bbtool_bench).
- `aln.jsonl`: gg.py rows (`variant` linsi, recipe = wsoft0.03:linsi&fftns2#es4, linsi#es3, linsi&fftns2-op3,
  linsi#es4; no `B` key) followed by vote.py rows from `cs581/gcmvote/code/run.py` (`B` = 10; magus, hard,
  hard-bb, soft4). FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the true tree). Methods `true`,
  `magus` (MAGUS end-to-end output), `es4` (linsi#es4), `recipe`, `vote_hard`, `vote_hard-bb`; run as 4
  independent trees.py processes at a time (`scripts/trees.sh`, xargs -P 4).

Vote rep dir: /opt/work/gcmvote/reps/SIMHIGH_Rk with `inputs`, `true.fasta`, `unaligned.fasta` symlinked from the
gcmtrees rep (separate results.jsonl). Check: vote `magus` SP error equals gg.py `linsi` (MAGUS's merge) on each
replicate. The six treed alignments per replicate have distinct md5s.

Bank: cs581/gcmvote/bank/SIMHIGH_Rk.tar.gz (+ MANIFEST.tsv). The gcmtrees rep has no subsets.json, so none is
banked (as in t20); subset membership is in inputs/subalignments.

`scripts/`: the lane scripts used here (aln.sh, trees.sh, bank.sh, collect.sh).
