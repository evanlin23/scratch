AI-assisted (Claude), exploration code for CS581 project

SIMHIGH_R11 (helper t11). aln.jsonl: gg.py rows (gcmtrees variants + linsi#es4) followed by gcmvote run.py rows
(vote.py variants, B=10). magus.jsonl: the MAGUS draw. trees.jsonl: FastTree -lg -gamma nRF vs true tree
(protbench trees.py); vote methods are named vote_<variant> (vote_hard_mask = hard+mask out.masked.fasta).
vote.py `magus` output = MAGUS's merge (variants/linsi): same 1000 rows, identical sequences, only record order differs.
model_*.json: vote.py fitted mixture per variant.
No vote alignment was byte-identical (cmp) to a gcmtrees alignment, so every listed method has its own FastTree run.
