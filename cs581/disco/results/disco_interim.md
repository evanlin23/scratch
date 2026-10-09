
### Mean species-tree RF, gene trees = g_100

| condition | n | AL-Pro2 | AD-DISCO | AL-DISCO | AD-DISCOR | AL-DISCOR | AD-DISCOR-it2 | AL-DISCOR-it2 | AD-DISCOR-lca | AL-DISCOR-lca | AD-DISCOR-oracle | AL-DISCOR-oracle |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| default | 9 | 0.0782 | 0.0760 | 0.0918 | 0.0737 | 0.0907 | 0.0771 | 0.0930 | 0.1236 | 0.1134 | 0.0726 | 0.0918 |

**Paired comparisons, g_100, all replicates** (diff = mean RF(A) − RF(B); W/T/L = A better/tie/worse, tie = identical RF; two-sided Wilcoxon signed-rank on non-tied pairs)

| A vs B | n | mean A | mean B | diff | W/T/L | p |
|---|---|---|---|---|---|---|
| ASTRID-DISCOR vs ASTRID-DISCO | 9 | 0.0737 | 0.0760 | -0.0023 | 4/3/2 | 0.812 |
| ASTRAL-DISCOR vs ASTRAL-DISCO | 9 | 0.0907 | 0.0918 | -0.0011 | 3/2/4 | 1 |
| ASTRID-DISCOR vs ASTRAL-Pro2 | 9 | 0.0737 | 0.0782 | -0.0045 | 5/2/2 | 0.703 |
| ASTRAL-DISCOR vs ASTRAL-Pro2 | 9 | 0.0907 | 0.0782 | +0.0125 | 2/1/6 | 0.0859 |
| ASTRID-DISCO vs ASTRAL-Pro2 | 9 | 0.0760 | 0.0782 | -0.0023 | 4/3/2 | 0.625 |
| ASTRID-DISCOR-it2 vs ASTRID-DISCO | 9 | 0.0771 | 0.0760 | +0.0011 | 3/3/3 | 0.656 |
| ASTRAL-DISCOR-it2 vs ASTRAL-DISCO | 9 | 0.0930 | 0.0918 | +0.0011 | 3/2/4 | 0.609 |
| ASTRID-DISCOR-oracle vs ASTRID-DISCO | 9 | 0.0726 | 0.0760 | -0.0034 | 4/3/2 | 0.844 |
| ASTRAL-DISCOR-oracle vs ASTRAL-DISCO | 9 | 0.0918 | 0.0918 | +0.0000 | 3/3/3 | 0.969 |
| ASTRID-DISCOR-lca vs ASTRID-DISCO | 9 | 0.1236 | 0.0760 | +0.0476 | 1/0/8 | 0.0117 |
| ASTRAL-DISCOR-lca vs ASTRAL-DISCO | 9 | 0.1134 | 0.0918 | +0.0215 | 3/3/3 | 0.844 |

**Paired comparisons, g_100, held-out replicates 06,07,08,09,10** (diff = mean RF(A) − RF(B); W/T/L = A better/tie/worse, tie = identical RF; two-sided Wilcoxon signed-rank on non-tied pairs)

| A vs B | n | mean A | mean B | diff | W/T/L | p |
|---|---|---|---|---|---|---|
| ASTRID-DISCOR vs ASTRID-DISCO | 4 | 0.0587 | 0.0638 | -0.0051 | 2/1/1 | 0.75 |
| ASTRAL-DISCOR vs ASTRAL-DISCO | 4 | 0.0842 | 0.0893 | -0.0051 | 2/0/2 | 1 |
| ASTRID-DISCOR vs ASTRAL-Pro2 | 4 | 0.0587 | 0.0638 | -0.0051 | 2/1/1 | 1 |
| ASTRAL-DISCOR vs ASTRAL-Pro2 | 4 | 0.0842 | 0.0638 | +0.0204 | 0/0/4 | 0.125 |
| ASTRID-DISCO vs ASTRAL-Pro2 | 4 | 0.0638 | 0.0638 | -0.0000 | 2/0/2 | 1 |
| ASTRID-DISCOR-it2 vs ASTRID-DISCO | 4 | 0.0663 | 0.0638 | +0.0026 | 1/1/2 | 0.5 |
| ASTRAL-DISCOR-it2 vs ASTRAL-DISCO | 4 | 0.0893 | 0.0893 | +0.0000 | 2/0/2 | 1 |
| ASTRID-DISCOR-oracle vs ASTRID-DISCO | 4 | 0.0663 | 0.0638 | +0.0026 | 1/1/2 | 0.5 |
| ASTRAL-DISCOR-oracle vs ASTRAL-DISCO | 4 | 0.0918 | 0.0893 | +0.0026 | 2/0/2 | 1 |
| ASTRID-DISCOR-lca vs ASTRID-DISCO | 4 | 0.0995 | 0.0638 | +0.0357 | 1/0/3 | 0.25 |
| ASTRAL-DISCOR-lca vs ASTRAL-DISCO | 4 | 0.1327 | 0.0893 | +0.0434 | 1/2/1 | 1 |

**Per condition, g_100** (diff | W/T/L | p)

| condition | AD-DISCOR vs AD-DISCO | AL-DISCOR vs AL-DISCO | AD-DISCOR vs Pro2 | AD-oracle vs AD-DISCO |
|---|---|---|---|---|
| default | -0.0023 | 4/3/2 | 0.812 | -0.0011 | 3/2/4 | 1 | -0.0045 | 5/2/2 | 0.703 | -0.0034 | 4/3/2 | 0.844 |

### Tagging accuracy (pooled over replicates; pair acc = fraction of cross-species gene pairs whose ortholog/paralog call matches the locus tree)

| condition | genes | DISCO pair acc | DISCO-R pair acc | DISCO-R-lca | oracle-S | DISCO orth P/R | DISCO-R orth P/R | DISCO root acc* | DISCO-R root acc* |
|---|---|---|---|---|---|---|---|---|---|
| default | g_100 | 0.851 | 0.861 | 0.779 | 0.861 | 0.828/0.612 | 0.825/0.658 | – | – |

\* root acc: fraction of multi-copy true gene trees (input rooting randomised) whose chosen root bipartition equals SimPhy's.

### Runtime (mean seconds per replicate, 1 thread)

| condition | genes | astral-pro2 | DISCO-decomp | DISCO-astrid | DISCO-astral | s0-root | DISCOR-decomp | DISCOR-astrid | DISCOR-astral |
|---|---|---|---|---|---|---|---|---|---|
| default | g_100 | 5.6 | 0.5 | 0.0 | 7.0 | 30.0 | 0.4 | 0.0 | 5.4 |

S0 root correct (min-DL rooting of ASTRAL-Pro2 tree): 6/9 runs.
