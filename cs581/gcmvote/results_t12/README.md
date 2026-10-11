AI-assisted (Claude), exploration code for CS581 project

SIMHIGH_R12 (helper t12). One MAGUS draw (paper flags) via cs581/gcmtrees/code/run.sh, rep /opt/work/gcmtrees/reps/SIMHIGH_R12.

- aln.jsonl: gg.py rows (variant linsi = MAGUS's merge, linsi#es4, recipe, es3, hard) followed by vote.py rows
  (run.py format, with "B": 10) for magus, es4, hard, soft, soft2, soft4, hard-bb, soft-bb, hard+mask. The vote
  variants ran on /opt/work/gcmvote/reps/SIMHIGH_R12, which symlinks inputs/, sets/ and true.fasta of the gcmtrees rep.
- magus.jsonl: bbtool_bench rows (end-to-end MAGUS and merge-mafft).
- trees.jsonl: protbench trees.py rows (FastTree -lg -gamma, nRF = RF vs the true tree). Scheduling: treesched.py.
- *.model.json: fitted vote models.
- Check: vote.py magus gives the same alignment (every sequence's row identical) as gg.py linsi and bbtool_bench
  merge-mafft; MAGUS's own end-to-end magus.fasta differs from those in 3 of 9417 columns with identical FastSP scores.
- hard+mask removed 5 of 9683 columns.
