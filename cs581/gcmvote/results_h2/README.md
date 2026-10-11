AI-assisted (Claude), exploration code for CS581 project

Helper h2 rows for cs581/gcmvote, SIMHIGH_R9 only (R10–R12 were reassigned to the t10/t11/t12 helpers).
- magus.jsonl: the MAGUS draw (draw 0, paper flags), from gcmtrees/code/run.sh.
- aln.jsonl: gcmgen merge-only variants on that draw (linsi = MAGUS merge, recipe, es3, hard filter, linsi#es4).
- vote_aln.jsonl: cs581/gcmvote/code/run.py --B 10 on the same inputs (magus es4 hard hard-bb soft soft-bb soft2 soft4 hard+mask).
- trees.jsonl: FastTree -lg -gamma (protbench/code/trees.py), RF = nRF vs the true tree. Methods: true, magus, es4
  (gcmgen linsi#es4), recipe, es3, hard (gcmtrees hard filter, linsi&fftns2-op3), vote_<variant>
  (vote_hard_mask = hard+mask's out.masked.fasta). vote magus / vote es4 were not treed; vote magus has the
  same FastSP score as magus and the same rows (sequence order in the file differs).
