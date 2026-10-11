AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R37, SIMHIGH_R38: vote-model tree rows (helper p37)

Two new replicates (AliSim, cs581/gcmtrees/code/gen.sh seeds), one MAGUS draw each (paper flags) via
`cs581/gcmtrees/code/run.sh SIMHIGH_Rk`, plus `gg.py run REP 'linsi#es4'`. Same procedure as results_t20;
`scripts/rep.sh K` is the whole lane for one replicate, `scripts/bank.sh K` banks it and copies the rows here.

- `magus.jsonl`: MAGUS end-to-end rows (bbtool_bench).
- `aln.jsonl`: per replicate, gg.py rows (`variant` linsi = MAGUS's merge, linsi#es4, recipe =
  wsoft0.03:linsi&fftns2#es4, linsi#es3, hard = linsi&fftns2-op3; no `B` key) then vote.py rows from
  `cs581/gcmvote/code/run.py` (`B` = 10; magus, hard, hard-bb, soft4). FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the true tree). Methods `true`,
  `magus` (gg.py linsi), `es4`, `recipe`, `vote_hard`, `vote_hard-bb`. Four FastTree processes at a time (xargs -P 4);
  R37 trees overlapped with R38's MAGUS run, so `fasttree_wall` is inflated there.

Vote rep dirs: /opt/work/gcmvote/reps/SIMHIGH_Rk with `inputs` and `true.fasta` symlinked from the gcmtrees rep.

Checks (scripts/norm.py: md5 of the alignment with rows sorted by name):
- R37: vote `magus` equals gg.py `linsi` (MAGUS's merge) up to row order; SPFN/SPFP/TC identical.
  vote `hard` and `hard-bb` are the same alignment up to row order (same 99352 kept edges); both were treed
  separately (FastTree input row order differs) and gave the same nRF.
