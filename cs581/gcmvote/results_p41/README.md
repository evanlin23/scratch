AI-assisted (Claude), exploration code for CS581 project

# Helper p41: SIMHIGH_R41, SIMHIGH_R42 (vote model, tree accuracy)

New AliSim replicates (cs581/gcmtrees/code/gen.sh seeds: tree 100r+7, sequences 100r+13), one MAGUS draw
(paper flags) each via `cs581/gcmtrees/code/run.sh`, plus `gg.py run REP 'linsi#es4'`. Procedure:
`scripts/pipe.sh NAME` (copied from helper p25), run one replicate at a time on 4 cores.

- `magus.jsonl`: bbtool_bench row of the draw.
- `aln.jsonl`: per replicate, gg.py rows (`linsi` = MAGUS merge, `linsi#es4`, recipe =
  `wsoft0.03:linsi&fftns2#es4`, `linsi#es3`, `linsi&fftns2-op3`; no `B` key) followed by run.py vote rows
  (`B` = 10: magus, hard, hard-bb, soft4). FastSP against the true alignment (bbtool_bench.acc_ref).
- `trees.jsonl`: protbench trees.py rows (FastTree -lg -gamma, `RF` = nRF vs the AliSim tree), 4 trees at a
  time as independent processes. magus = gg `linsi`, es4 = gg `linsi#es4`, recipe = gg recipe,
  vote_X = vote.py X, true = true alignment.
- Vote rep dir /opt/work/gcmvote/reps/NAME links `inputs` and `true.fasta` from the gcmtrees rep (separate
  results.jsonl: run.py needs a `B` key on every row).
- Bank: `cs581/gcmvote/bank/NAME.tar.gz` (`scripts/bank.sh`), same layout as banks t18-t20: inputs/, sets/,
  aligned/, true.fasta, true_tree.nwk, unaligned.fasta, magus.json. No subsets.json: no gcmtrees rep has one
  (vote.py writes it per vote run, with absolute paths), and the t18-t20 banks have none either.

Checks (`scripts/norm.py`, sequence content after sorting rows by name):
- SIMHIGH_R41: vote `magus` has the same content as gg `linsi` (MAGUS's merge); same SPFN/SPFP/TC.
  vote `hard` and `hard-bb` have identical content (same 92062 kept edges); both were treed anyway
  (separate FastTree runs, same nRF). No other treed pair is identical.
