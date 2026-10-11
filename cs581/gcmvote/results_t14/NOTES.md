AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R14 (helper t14)

- Rep: /opt/work/gcmtrees/reps/SIMHIGH_R14 (one MAGUS draw, paper flags, via gcmtrees/code/run.sh), gg.py
  variants + `linsi#es4`. Vote variants run with gcmvote/code/run.py on a separate rep dir
  /opt/work/gcmvote/reps/SIMHIGH_R14 that symlinks inputs/, sets/, true.fasta etc. (run.py's results.jsonl
  cannot be shared with gg.py rows, which have no "B").
- aln.jsonl: first 5 rows gg.py format (linsi = MAGUS, wsoft0.03:linsi&fftns2#es4 = recipe, linsi#es3,
  linsi&fftns2-op3 = hard, linsi#es4), then 9 run.py rows (vote variants, B = 10).
- vote.py `magus` reproduces MAGUS's merge: same rows as variants/linsi/out.fasta (row order differs, so not
  byte-identical; same.py), identical SPFN/SPFP.
- No vote alignment is byte-identical to an already treed one, so every requested tree was run.
  vote/hard+mask/out.fasta has the same rows as vote/hard/out.fasta; the tree uses out.masked.fasta (1 column
  removed).
- Trees: FastTree -lg -gamma via protbench/code/trees.py, at most 4 at once (tree1.sh for the gcmtrees set,
  tree2.sh slot-gated for the vote set). The `true` tree used /opt/data/sim/SIMHIGH/R14/sim.fa, which has the
  same rows as the rep's true.fasta (header padding differs only).
