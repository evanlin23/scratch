
### Mean species-tree RF, gene trees = g_100

| condition | n | AL-Pro2 | AD-DISCO | AL-DISCO | AD-DISCOR | AL-DISCOR | AD-DISCOR-it2 | AL-DISCOR-it2 | AD-DISCOR-lca | AL-DISCOR-lca | AD-DISCOR-oracle | AL-DISCOR-oracle |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| default | 14 | 0.0700 | 0.0649 | 0.0809 | 0.0641 | 0.0802 | 0.0827 | 0.1020 | 0.1101 | 0.1130 | 0.0620 | 0.0802 |
| gdl_1e-10_1 | 10 | 0.0939 | 0.0786 | 0.0908 | 0.0745 | 0.0949 | 0.0745 | 0.0949 | 0.1429 | 0.1500 | 0.0776 | 0.0949 |
| gdl_1e-9_1 | 10 | 0.0837 | 0.0898 | 0.1061 | 0.0796 | 0.1082 | 0.0765 | 0.1010 | 0.1429 | 0.1408 | 0.0796 | 0.1102 |
| ils_1e4 | 10 | 0.0653 | 0.0592 | 0.0735 | 0.0520 | 0.0663 | 0.0510 | 0.0663 | 0.0898 | 0.0724 | 0.0510 | 0.0663 |
| ils_2e8 | 10 | 0.1735 | 0.1347 | 0.1827 | 0.1337 | 0.1786 | 0.1347 | 0.1796 | 0.2520 | 0.3041 | 0.1327 | 0.1796 |

**Paired comparisons, g_100, all replicates** (diff = mean RF(A) − RF(B); W/T/L = A better/tie/worse, tie = identical RF; two-sided Wilcoxon signed-rank on non-tied pairs)

| A vs B | n | mean A | mean B | diff | W/T/L | p |
|---|---|---|---|---|---|---|
| ASTRID-DISCOR vs ASTRID-DISCO | 54 | 0.0796 | 0.0839 | -0.0043 | 24/19/11 | 0.0593 |
| ASTRAL-DISCOR vs ASTRAL-DISCO | 54 | 0.1037 | 0.1049 | -0.0011 | 22/13/19 | 0.64 |
| ASTRID-DISCOR vs ASTRAL-Pro2 | 54 | 0.0796 | 0.0952 | -0.0157 | 34/6/14 | 0.000289 |
| ASTRAL-DISCOR vs ASTRAL-Pro2 | 54 | 0.1037 | 0.0952 | +0.0085 | 14/13/27 | 0.0252 |
| ASTRID-DISCO vs ASTRAL-Pro2 | 54 | 0.0839 | 0.0952 | -0.0113 | 32/8/14 | 0.00796 |
| ASTRID-DISCOR-it2 vs ASTRID-DISCO | 50 | 0.0839 | 0.0888 | -0.0049 | 22/17/11 | 0.0985 |
| ASTRAL-DISCOR-it2 vs ASTRAL-DISCO | 50 | 0.1088 | 0.1106 | -0.0018 | 22/9/19 | 0.682 |
| ASTRID-DISCOR-oracle vs ASTRID-DISCO | 54 | 0.0792 | 0.0839 | -0.0047 | 25/16/13 | 0.056 |
| ASTRAL-DISCOR-oracle vs ASTRAL-DISCO | 54 | 0.1043 | 0.1049 | -0.0006 | 22/15/17 | 0.552 |
| ASTRID-DISCOR-lca vs ASTRID-DISCO | 54 | 0.1447 | 0.0839 | +0.0608 | 6/4/44 | 1.91e-08 |
| ASTRAL-DISCOR-lca vs ASTRAL-DISCO | 54 | 0.1529 | 0.1049 | +0.0480 | 18/13/23 | 0.0225 |

**Paired comparisons, g_100, held-out replicates 06,07,08,09,10** (diff = mean RF(A) − RF(B); W/T/L = A better/tie/worse, tie = identical RF; two-sided Wilcoxon signed-rank on non-tied pairs)

| A vs B | n | mean A | mean B | diff | W/T/L | p |
|---|---|---|---|---|---|---|
| ASTRID-DISCOR vs ASTRID-DISCO | 25 | 0.0890 | 0.0939 | -0.0049 | 14/5/6 | 0.133 |
| ASTRAL-DISCOR vs ASTRAL-DISCO | 25 | 0.1204 | 0.1176 | +0.0029 | 8/6/11 | 0.364 |
| ASTRID-DISCOR vs ASTRAL-Pro2 | 25 | 0.0890 | 0.1078 | -0.0188 | 16/4/5 | 0.00698 |
| ASTRAL-DISCOR vs ASTRAL-Pro2 | 25 | 0.1204 | 0.1078 | +0.0127 | 4/5/16 | 0.0331 |
| ASTRID-DISCO vs ASTRAL-Pro2 | 25 | 0.0939 | 0.1078 | -0.0139 | 16/2/7 | 0.0413 |
| ASTRID-DISCOR-it2 vs ASTRID-DISCO | 25 | 0.0890 | 0.0939 | -0.0049 | 12/6/7 | 0.343 |
| ASTRAL-DISCOR-it2 vs ASTRAL-DISCO | 25 | 0.1188 | 0.1176 | +0.0012 | 9/5/11 | 0.348 |
| ASTRID-DISCOR-oracle vs ASTRID-DISCO | 25 | 0.0898 | 0.0939 | -0.0041 | 13/4/8 | 0.264 |
| ASTRAL-DISCOR-oracle vs ASTRAL-DISCO | 25 | 0.1220 | 0.1176 | +0.0045 | 8/7/10 | 0.371 |
| ASTRID-DISCOR-lca vs ASTRID-DISCO | 25 | 0.1653 | 0.0939 | +0.0714 | 2/2/21 | 6.72e-05 |
| ASTRAL-DISCOR-lca vs ASTRAL-DISCO | 25 | 0.1849 | 0.1176 | +0.0673 | 4/8/13 | 0.00597 |

**Per condition, g_100** (diff | W/T/L | p)

| condition | AD-DISCOR vs AD-DISCO | AL-DISCOR vs AL-DISCO | AD-DISCOR vs Pro2 | AD-oracle vs AD-DISCO |
|---|---|---|---|---|
| default | -0.0007 | 4/7/3 | 1 | -0.0007 | 4/5/5 | 0.82 | -0.0058 | 8/3/3 | 0.327 | -0.0029 | 5/6/3 | 0.703 |
| gdl_1e-10_1 | -0.0041 | 5/3/2 | 0.469 | +0.0041 | 2/3/5 | 0.344 | -0.0194 | 7/1/2 | 0.0938 | -0.0010 | 4/2/4 | 0.922 |
| gdl_1e-9_1 | -0.0102 | 4/5/1 | 0.188 | +0.0020 | 5/1/4 | 0.883 | -0.0041 | 4/0/6 | 0.736 | -0.0102 | 4/5/1 | 0.188 |
| ils_1e4 | -0.0071 | 7/1/2 | 0.285 | -0.0071 | 6/2/2 | 0.109 | -0.0133 | 6/1/3 | 0.289 | -0.0082 | 7/1/2 | 0.242 |
| ils_2e8 | -0.0010 | 4/3/3 | 0.812 | -0.0041 | 5/2/3 | 0.711 | -0.0398 | 9/1/0 | 0.00391 | -0.0020 | 5/2/3 | 0.672 |

### Tagging accuracy (pooled over replicates; pair acc = fraction of cross-species gene pairs whose ortholog/paralog call matches the locus tree)

| condition | genes | DISCO pair acc | DISCO-R pair acc | DISCO-R-lca | oracle-S | DISCO orth P/R | DISCO-R orth P/R | DISCO root acc* | DISCO-R root acc* |
|---|---|---|---|---|---|---|---|---|---|
| default | g_100 | 0.842 | 0.851 | 0.762 | 0.853 | 0.833/0.605 | 0.829/0.647 | – | – |
| gdl_1e-10_1 | g_100 | 0.899 | 0.914 | 0.604 | 0.920 | 0.969/0.876 | 0.966/0.902 | – | – |
| gdl_1e-9_1 | g_100 | 0.870 | 0.862 | 0.829 | 0.862 | 0.713/0.527 | 0.655/0.598 | – | – |
| ils_1e4 | g_100 | 0.878 | 0.887 | 0.815 | 0.888 | 0.857/0.701 | 0.843/0.759 | – | – |
| ils_2e8 | g_100 | 0.830 | 0.847 | 0.778 | 0.847 | 0.821/0.555 | 0.814/0.627 | – | – |

\* root acc: fraction of multi-copy true gene trees (input rooting randomised) whose chosen root bipartition equals SimPhy's.

### Runtime (mean seconds per replicate, 1 thread)

| condition | genes | astral-pro2 | DISCO-decomp | DISCO-astrid | DISCO-astral | s0-root | DISCOR-decomp | DISCOR-astrid | DISCOR-astral |
|---|---|---|---|---|---|---|---|---|---|
| default | g_100 | 75.8 | 3.1 | 0.0 | 78.4 | 38.6 | 3.5 | 0.1 | 91.7 |
| gdl_1e-10_1 | g_100 | 3.3 | 0.2 | 0.0 | 5.7 | 16.4 | 0.2 | 0.0 | 4.8 |
| gdl_1e-9_1 | g_100 | 5.3 | 0.5 | 0.0 | 5.3 | 32.9 | 0.4 | 0.0 | 5.1 |
| ils_1e4 | g_100 | 3.7 | 0.3 | 0.0 | 4.8 | 22.6 | 0.3 | 0.0 | 4.5 |
| ils_2e8 | g_100 | 4.9 | 0.3 | 0.0 | 5.8 | 25.5 | 0.3 | 0.0 | 5.6 |

S0 root correct (min-DL rooting of ASTRAL-Pro2 tree): 40/54 runs.
