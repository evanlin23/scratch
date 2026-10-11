# Held-out plan (written 2026-10-09 ~23:06 UTC, after the training condition finished)

Training condition: n25 FastTree reps 20-39. This is the same condition the paper used to pick t = 0.5.
Held-out conditions: n15 IQ-TREE reps 00-19, n25 IQ-TREE reps 00-19, n50 FastTree reps 00-19.
Caveat: 3-7 n15 replicates had already been aggregated (for monitoring) when this was written.

Primary held-out claims to test, chosen from the training results only:
1. **Adaptive filter:** `z3only` keeps the 2nd quartet topology iff (c2-c3-1)/sqrt(c2+c3) > 3, with no ratio rule.
   It had the lowest training mean error among the z variants (-0.0096 vs default, W/T/L 8/9/3, p=0.45). Secondary: `z5only` (5/15/0, p=0.06).
2. **Alternative base tree:** `tqmc` (TREE-QMC). It had the best training mean among the non-oracle trees (-0.0008, p=0.64).
3. **Score-based base-tree selection** among {default, tqmc, wastral, swap}.

Not claims; reported for context: the `true_major` oracle and the t sweep.
