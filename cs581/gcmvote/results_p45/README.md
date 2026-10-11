AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R45, SIMHIGH_R46: vote-model replicates (helper p45)

Same procedure as results_t20. Per replicate (AliSim, cs581/gcmtrees/code/gen.sh seeds): one MAGUS draw (paper
flags) and the gcmtrees merge variants via `cs581/gcmtrees/code/run.sh`, then `gg.py run REP 'linsi#es4'`;
vote.py variants via `cs581/gcmvote/code/run.py` in a separate rep dir /opt/work/gcmvote/reps/<rep> (inputs,
true.fasta, unaligned.fasta symlinked from the gcmtrees rep); FastTree -lg -gamma via protbench trees.py,
4 independent processes at a time. Driver: `scripts/rep.sh <rep>`.

- `magus.jsonl`: MAGUS end-to-end row (bbtool_bench).
- `aln.jsonl`: gg.py rows (`variant` linsi = MAGUS's merge, wsoft0.03:linsi&fftns2#es4 = recipe, linsi#es3,
  linsi&fftns2-op3, linsi#es4; no `B` key) followed by vote.py rows (`B` = 10; magus, hard, hard-bb, soft4).
  FastSP against the true alignment.
- `trees.jsonl`: protbench trees.py rows, nRF = `RF` vs the true tree. Methods true, magus, es4, recipe,
  vote_hard, vote_hard-bb.

Checks: vote `magus` out.fasta equals gg `linsi` out.fasta after sorting rows by name (`scripts/norm.py`), and
their SPFN/SPFP/TC are identical.

Bank: cs581/gcmvote/bank/<rep>.tar.gz (layout as in bank/README.md). No subsets.json is written by this lane
(the subset membership is the 25 files in inputs/subalignments), so none is in the tarball.
