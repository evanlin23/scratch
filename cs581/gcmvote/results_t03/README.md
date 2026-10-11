AI-assisted (Claude), exploration code for CS581 project

# SIMHIGH_R3 helper run (vote model + gcmtrees set)

One MAGUS draw (draw 0, paper flags, 25 subsets, 10 L-INS-i backbones), AliSim SIMHIGH R3 from gcmtrees/code/gen.sh.

- `magus.jsonl`: the MAGUS draw (bbtool_bench row).
- `aln.jsonl`: FastSP rows. First 5 rows are gg.py/pc.py variants (`linsi` = magus, `wsoft0.03:linsi&fftns2#es4` = recipe,
  `linsi#es3`, `linsi&fftns2-op3` = hard, `linsi#es4` = es4); the rest are vote.py rows (run.py format, B = 10).
- `trees.jsonl`: FastTree -lg -gamma (protbench/code/trees.py) vs the true tree. `vote_hard_mask` = hard+mask's out.masked.fasta.
- `raw/`: vote.py model.json per variant.

Check: vote.py `magus` has the same 1000 rows/strings as MAGUS's merge (`variants/linsi/out.fasta`); the files differ only in
row order, so `cmp` differs but scores are identical.

Trees: no vote alignment was byte-identical to an earlier one, so every tree was built (no reuse). The first magus, recipe
and es3 FastTree runs started on alignment files still being written (es3 tree had 354/1000 leaves); those runs were
discarded and the three trees re-run on the complete files.
