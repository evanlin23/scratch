AI-assisted (Claude), exploration code for CS581 project

SIMHIGH_R16 (helper t16). One MAGUS draw, paper flags (gcmtrees run.sh).
- aln.jsonl: gg.py rows (no "B" key: linsi = magus, linsi#es4, recipe = wsoft0.03:linsi&fftns2#es4, linsi#es3,
  hard = linsi&fftns2-op3) then vote.py rows via gcmvote/code/run.py (with "B"). FastSP vs true.fasta.
- vote.py `magus` = MAGUS's merge (same rows, only sequence order differs; same FastSP scores).
  vote.py `es4` (16.53 % err) is not byte-identical to gg.py `linsi#es4` (16.58 %).
- trees.jsonl: FastTree -lg -gamma via protbench/code/trees.py, nRF vs the true tree. Methods: true, magus, es4
  (gg linsi#es4), recipe, es3, hard (gcmtrees set), vote_<variant>; vote_hard+mask = its out.masked.fasta.
- magus.jsonl: bbtool_bench row of the MAGUS draw.
