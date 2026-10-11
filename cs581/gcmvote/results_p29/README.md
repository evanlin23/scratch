AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R29, SIMHIGH_R30: vote-model tree rows (helper p29)

Two replicates (AliSim, cs581/gcmtrees/code/gen.sh seeds), one MAGUS draw each (paper flags) via
`cs581/gcmtrees/code/run.sh SIMHIGH_Rk`, plus `gg.py run REP 'linsi#es4'`.

- `magus.jsonl`: MAGUS end-to-end rows (bbtool_bench).
- `aln.jsonl`: per replicate, gg.py rows (`variant` linsi = MAGUS's merge, recipe = wsoft0.03:linsi&fftns2#es4,
  linsi#es3, linsi&fftns2-op3, linsi#es4; no `B` key) followed by vote.py rows from `cs581/gcmvote/code/run.py`
  (`B` = 10; magus, hard, hard-bb, soft4). FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, nRF = `RF` vs the true tree). Methods `true`,
  `magus` (gg linsi), `es4` (gg linsi#es4), `recipe`, `vote_hard`, `vote_hard-bb`.

Vote rep dir: /opt/work/gcmvote/reps/SIMHIGH_Rk with `inputs` and `true.fasta` symlinked from the gcmtrees rep.

Checks: vote.py `magus` equals gg.py `linsi` (MAGUS's merge) after sorting rows by name (`scripts/norm.py`;
same SPFN/SPFP/TC). The six treed alignments of a replicate are pairwise different (norm.py hashes).

`scripts/`: `rep.sh` (vote variants + 6 FastTree runs, 4 at a time via `tree1.sh`), `bank.sh` (bank tarball +
MANIFEST row), `collect.sh` (copies rows here).
