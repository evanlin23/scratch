AI-assisted (Claude), exploration code for CS581 project

Rep bank from helper h9 (ROSE 1000M2 R0, 1000L1 R0). One MAGUS draw per dataset (paper flags, 25 subsets,
20 L-INS-i backbones x 200 sequences); `<DS>_B5` / `_B10` hold backbones 1-5 / 1-10 of the `_B20` draw on the
same subsets. Each tar: `<name>/inputs/{subalignments,backbones}`, `true.fasta` (reference, uppercased),
`unaligned.fasta`, `sets/s0` (unaligned backbone sequence sets), `magus.json` (end-to-end 20-backbone MAGUS
score), `true_tree.tt` (ROSE `rose.tt`). Extract under any directory and pass `<name>` as REP_DIR to
gg.py run / vote.py (bbe.prep is already applied; `aligned/linsi/s0` is rebuilt by copying inputs/backbones).
