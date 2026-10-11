AI-assisted (Claude), exploration code for CS581 project

SIMHIGH_R7 (helper t07; took over from gcmvote-h1). One MAGUS draw (paper flags) via cs581/gcmtrees/code/run.sh.

- aln.jsonl: gg.py rows (variant = linsi, linsi#es4, linsi#es3, wsoft0.03:linsi&fftns2#es4 = recipe,
  linsi&fftns2-op3 = hard) then gcmvote/code/run.py rows (variant = magus, es4, hard, soft, soft2, soft4,
  hard-bb, soft-bb, hard+mask; B = 10). Vote rows ran on /opt/work/gcmvote/reps/SIMHIGH_R7, whose inputs are
  symlinks to the gcmtrees rep (MAGUS's own 25 subalignments and 10 L-INS-i backbones).
- vote `magus` = gg `linsi` (MAGUS's merge): same sequences/columns, rows in a different order (not cmp-identical,
  identical after sorting by name; identical FastSP scores).
- vote `es4` vs gg `linsi#es4`: not identical (10530 vs 10526 columns; avgErr 0.194924 vs 0.194928).
- vote `hard` and `hard+mask` out.fasta are identical after sorting; hard+mask is treed from out.masked.fasta.
- magus.jsonl: bbtool_bench row of the MAGUS draw (full MAGUS, and merge-mafft).
- trees.jsonl: protbench/code/trees.py rows, FastTree -lg -gamma, nRF vs the true tree. Methods true, magus
  (= gg linsi), es4 (= gg linsi#es4), recipe, es3, hard (gcmtrees set) and vote_<variant>; vote_hard+mask is
  the masked alignment. No vote alignment was cmp-identical to an already-treed one, so none reuse a tree.
