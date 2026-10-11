AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R33, SIMHIGH_R34: vote-model tree rows (helper p33)

Replicates SIMHIGH_R33 and R34 (AliSim, cs581/gcmtrees/code/gen.sh seeds), one MAGUS draw each (paper flags) via
`cs581/gcmtrees/code/run.sh SIMHIGH_Rk`, plus `gg.py run REP 'linsi#es4'`. Same procedure as results_t20;
driver `scripts/rep.sh`, bank tarball builder `scripts/bank.sh`.

- `magus.jsonl`: MAGUS end-to-end rows (bbtool_bench).
- `aln.jsonl`: per replicate, gg.py rows (`variant` linsi, wsoft0.03:linsi&fftns2#es4 = recipe, linsi#es3,
  linsi&fftns2-op3, linsi#es4; no `B` key) followed by vote.py rows from `cs581/gcmvote/code/run.py`
  (`B` = 10; magus, hard, hard-bb, soft4). FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the true tree). Methods `true`,
  `magus` (= gg linsi), `es4` (= gg linsi#es4), `recipe`, `vote_hard`, `vote_hard-bb`. 4 FastTree processes at a time.

Vote rep dir: /opt/work/gcmvote/reps/SIMHIGH_Rk with `inputs` and `true.fasta` symlinked from the gcmtrees rep.

Checks (SIMHIGH_R33): vote `magus` equals gg `linsi` (MAGUS's merge) up to row order (`scripts/norm.py`; same
SPFN/SPFP/TC). vote `hard` and `hard-bb` give the same alignment up to row order (same cutoff_by_n and
kept_edges in model.json); both were treed separately anyway and gave the same nRF.
