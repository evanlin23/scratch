AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R27, SIMHIGH_R28: vote-model tree replicates (helper p27)

Two replicates (AliSim, `cs581/gcmtrees/code/gen.sh` seeds: tree 100r+7, sequences 100r+13), one MAGUS draw each
(paper flags) via `cs581/gcmtrees/code/run.sh SIMHIGH_Rk`, then `gg.py run REP 'linsi#es4'`. Same procedure as
results_t20.

- `magus.jsonl`: MAGUS end-to-end rows (bbtool_bench).
- `aln.jsonl`: per replicate, gg.py rows (`variant` linsi = MAGUS's merge, recipe = wsoft0.03:linsi&fftns2#es4,
  linsi#es3, linsi&fftns2-op3, linsi#es4; no `B` key) followed by vote.py rows from `cs581/gcmvote/code/run.py`
  (`B` = 10; magus, hard, hard-bb, soft4). FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the AliSim true tree). Methods `true`,
  `magus` (gg linsi), `es4` (gg linsi#es4), `recipe`, `vote_hard`, `vote_hard-bb`. 4 FastTree processes at a time
  (independent processes, niced).

Vote rep dir: /opt/work/gcmvote/reps/SIMHIGH_Rk with `inputs` and `true.fasta` symlinked from the gcmtrees rep
(separate results.jsonl because run.py expects a `B` key on every row).

Checks: vote `magus` out.fasta equals gg.py `linsi` (MAGUS's merge) after sorting rows by name
(`scripts/norm.py`, same md5) and has identical SPFN/SPFP/TC, on both replicates.

Bank: `cs581/gcmvote/bank/SIMHIGH_Rk.tar.gz` + MANIFEST.tsv row (merge-only inputs; no subsets.json is written by
this pipeline, so none is included). `scripts/`: pipe.sh (whole replicate), bank.sh, collect.sh, norm.py.
