
### Mean species-tree RF, gene trees = g_100

| condition | n | AL-Pro2 | AD-DISCO | AL-DISCO | AD-DISCOR | AL-DISCOR | AD-DISCOR-it2 | AL-DISCOR-it2 | AD-DISCOR-lca | AL-DISCOR-lca | AD-DISCOR-oracle | AL-DISCOR-oracle |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20_gdl_1e-10_0 | 10 | 0.0333 | 0.0333 | 0.0333 | 0.0333 | 0.0333 | 0.0333 | 0.0333 | 0.0333 | 0.0333 | 0.0333 | 0.0333 |
| 20_gdl_1e-10_0_hILS | 5 | 0.0667 | 0.0556 | 0.0667 | 0.0556 | 0.0667 | 0.0556 | 0.0667 | 0.0667 | 0.0667 | 0.0556 | 0.0667 |
| 20_gdl_5e-10_1 | 7 | 0.0317 | 0.0238 | 0.0317 | 0.0238 | 0.0397 | 0.0238 | 0.0397 | 0.0317 | 0.0317 | 0.0238 | 0.0397 |
| 20_gdl_5e-10_1_hILS | 10 | 0.0333 | 0.0389 | 0.0444 | 0.0389 | 0.0389 | 0.0389 | 0.0389 | 0.0333 | 0.0333 | 0.0389 | 0.0389 |

**Paired comparisons, g_100, all replicates** (diff = mean RF(A) − RF(B); W/T/L = A better/tie/worse, tie = identical RF; two-sided Wilcoxon signed-rank on non-tied pairs)

| A vs B | n | mean A | mean B | diff | W/T/L | p |
|---|---|---|---|---|---|---|
| ASTRID-DISCOR vs ASTRID-DISCO | 32 | 0.0365 | 0.0365 | +0.0000 | 1/30/1 | 1 |
| ASTRAL-DISCOR vs ASTRAL-DISCO | 32 | 0.0417 | 0.0417 | +0.0000 | 1/30/1 | 1 |
| ASTRID-DISCOR vs ASTRAL-Pro2 | 32 | 0.0365 | 0.0382 | -0.0017 | 4/25/3 | 1 |
| ASTRAL-DISCOR vs ASTRAL-Pro2 | 32 | 0.0417 | 0.0382 | +0.0035 | 0/30/2 | 0.5 |
| ASTRID-DISCO vs ASTRAL-Pro2 | 32 | 0.0365 | 0.0382 | -0.0017 | 4/25/3 | 1 |
| ASTRID-DISCOR-it2 vs ASTRID-DISCO | 32 | 0.0365 | 0.0365 | +0.0000 | 1/30/1 | 1 |
| ASTRAL-DISCOR-it2 vs ASTRAL-DISCO | 32 | 0.0417 | 0.0417 | +0.0000 | 1/30/1 | 1 |
| ASTRID-DISCOR-oracle vs ASTRID-DISCO | 32 | 0.0365 | 0.0365 | +0.0000 | 1/30/1 | 1 |
| ASTRAL-DISCOR-oracle vs ASTRAL-DISCO | 32 | 0.0417 | 0.0417 | +0.0000 | 1/30/1 | 1 |
| ASTRID-DISCOR-lca vs ASTRID-DISCO | 32 | 0.0382 | 0.0365 | +0.0017 | 3/25/4 | 1 |
| ASTRAL-DISCOR-lca vs ASTRAL-DISCO | 32 | 0.0382 | 0.0417 | -0.0035 | 2/30/0 | 0.5 |

**Paired comparisons, g_100, held-out replicates 06,07,08,09,10** (diff = mean RF(A) − RF(B); W/T/L = A better/tie/worse, tie = identical RF; two-sided Wilcoxon signed-rank on non-tied pairs)

| A vs B | n | mean A | mean B | diff | W/T/L | p |
|---|---|---|---|---|---|---|
| ASTRID-DISCOR vs ASTRID-DISCO | 12 | 0.0231 | 0.0185 | +0.0046 | 0/11/1 | 1 |
| ASTRAL-DISCOR vs ASTRAL-DISCO | 12 | 0.0278 | 0.0278 | +0.0000 | 1/10/1 | 1 |
| ASTRID-DISCOR vs ASTRAL-Pro2 | 12 | 0.0231 | 0.0185 | +0.0046 | 2/7/3 | 1 |
| ASTRAL-DISCOR vs ASTRAL-Pro2 | 12 | 0.0278 | 0.0185 | +0.0093 | 0/10/2 | 0.5 |
| ASTRID-DISCO vs ASTRAL-Pro2 | 12 | 0.0185 | 0.0185 | +0.0000 | 2/8/2 | 1 |
| ASTRID-DISCOR-it2 vs ASTRID-DISCO | 12 | 0.0231 | 0.0185 | +0.0046 | 0/11/1 | 1 |
| ASTRAL-DISCOR-it2 vs ASTRAL-DISCO | 12 | 0.0278 | 0.0278 | +0.0000 | 1/10/1 | 1 |
| ASTRID-DISCOR-oracle vs ASTRID-DISCO | 12 | 0.0231 | 0.0185 | +0.0046 | 0/11/1 | 1 |
| ASTRAL-DISCOR-oracle vs ASTRAL-DISCO | 12 | 0.0278 | 0.0278 | +0.0000 | 1/10/1 | 1 |
| ASTRID-DISCOR-lca vs ASTRID-DISCO | 12 | 0.0185 | 0.0185 | +0.0000 | 2/8/2 | 1 |
| ASTRAL-DISCOR-lca vs ASTRAL-DISCO | 12 | 0.0185 | 0.0278 | -0.0093 | 2/10/0 | 0.5 |

**Per condition, g_100** (diff | W/T/L | p)

| condition | AD-DISCOR vs AD-DISCO | AL-DISCOR vs AL-DISCO | AD-DISCOR vs Pro2 | AD-oracle vs AD-DISCO |
|---|---|---|---|---|
| 20_gdl_1e-10_0 | +0.0000 | 0/10/0 | 1 | +0.0000 | 0/10/0 | 1 | +0.0000 | 0/10/0 | 1 | +0.0000 | 0/10/0 | 1 |
| 20_gdl_1e-10_0_hILS | +0.0000 | 0/5/0 | 1 | +0.0000 | 0/5/0 | 1 | -0.0111 | 1/4/0 | 1 | +0.0000 | 0/5/0 | 1 |
| 20_gdl_5e-10_1 | +0.0000 | 0/7/0 | 1 | +0.0079 | 0/6/1 | 1 | -0.0079 | 1/6/0 | 1 | +0.0000 | 0/7/0 | 1 |
| 20_gdl_5e-10_1_hILS | +0.0000 | 1/8/1 | 1 | -0.0056 | 1/9/0 | 1 | +0.0056 | 2/5/3 | 1 | +0.0000 | 1/8/1 | 1 |

### Mean species-tree RF, gene trees = g_50

| condition | n | AL-Pro2 | AD-DISCO | AL-DISCO | AD-DISCOR | AL-DISCOR | AD-DISCOR-it2 | AL-DISCOR-it2 | AD-DISCOR-lca | AL-DISCOR-lca | AD-DISCOR-oracle | AL-DISCOR-oracle |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20_gdl_1e-10_0 | 10 | 0.0389 | 0.0333 | 0.0389 | 0.0333 | 0.0333 | 0.0333 | 0.0333 | 0.0389 | 0.0389 | 0.0333 | 0.0333 |
| 20_gdl_1e-10_0_hILS | 4 | 0.0833 | 0.1111 | 0.0972 | 0.0694 | 0.0694 | 0.0694 | 0.0694 | 0.0833 | 0.0833 | 0.0694 | 0.0694 |
| 20_gdl_5e-10_1 | 8 | 0.0556 | 0.0486 | 0.0556 | 0.0347 | 0.0556 | 0.0347 | 0.0556 | 0.0556 | 0.0556 | 0.0347 | 0.0556 |
| 20_gdl_5e-10_1_hILS | 10 | 0.0556 | 0.0500 | 0.0611 | 0.0500 | 0.0444 | 0.0500 | 0.0444 | 0.0500 | 0.0556 | 0.0500 | 0.0444 |

**Paired comparisons, g_50, all replicates** (diff = mean RF(A) − RF(B); W/T/L = A better/tie/worse, tie = identical RF; two-sided Wilcoxon signed-rank on non-tied pairs)

| A vs B | n | mean A | mean B | diff | W/T/L | p |
|---|---|---|---|---|---|---|
| ASTRID-DISCOR vs ASTRID-DISCO | 32 | 0.0434 | 0.0521 | -0.0087 | 5/25/2 | 0.25 |
| ASTRAL-DISCOR vs ASTRAL-DISCO | 32 | 0.0469 | 0.0573 | -0.0104 | 7/23/2 | 0.148 |
| ASTRID-DISCOR vs ASTRAL-Pro2 | 32 | 0.0434 | 0.0538 | -0.0104 | 7/22/3 | 0.186 |
| ASTRAL-DISCOR vs ASTRAL-Pro2 | 32 | 0.0469 | 0.0538 | -0.0069 | 5/26/1 | 0.219 |
| ASTRID-DISCO vs ASTRAL-Pro2 | 32 | 0.0521 | 0.0538 | -0.0017 | 5/22/5 | 1 |
| ASTRID-DISCOR-it2 vs ASTRID-DISCO | 32 | 0.0434 | 0.0521 | -0.0087 | 5/25/2 | 0.25 |
| ASTRAL-DISCOR-it2 vs ASTRAL-DISCO | 32 | 0.0469 | 0.0573 | -0.0104 | 7/23/2 | 0.148 |
| ASTRID-DISCOR-oracle vs ASTRID-DISCO | 32 | 0.0434 | 0.0521 | -0.0087 | 5/25/2 | 0.25 |
| ASTRAL-DISCOR-oracle vs ASTRAL-DISCO | 32 | 0.0469 | 0.0573 | -0.0104 | 7/23/2 | 0.148 |
| ASTRID-DISCOR-lca vs ASTRID-DISCO | 32 | 0.0521 | 0.0521 | +0.0000 | 6/21/5 | 1 |
| ASTRAL-DISCOR-lca vs ASTRAL-DISCO | 32 | 0.0538 | 0.0573 | -0.0035 | 4/27/1 | 0.75 |

**Paired comparisons, g_50, held-out replicates 06,07,08,09,10** (diff = mean RF(A) − RF(B); W/T/L = A better/tie/worse, tie = identical RF; two-sided Wilcoxon signed-rank on non-tied pairs)

| A vs B | n | mean A | mean B | diff | W/T/L | p |
|---|---|---|---|---|---|---|
| ASTRID-DISCOR vs ASTRID-DISCO | 13 | 0.0214 | 0.0299 | -0.0085 | 2/10/1 | 0.75 |
| ASTRAL-DISCOR vs ASTRAL-DISCO | 13 | 0.0256 | 0.0385 | -0.0128 | 3/10/0 | 0.25 |
| ASTRID-DISCOR vs ASTRAL-Pro2 | 13 | 0.0214 | 0.0299 | -0.0085 | 4/6/3 | 0.766 |
| ASTRAL-DISCOR vs ASTRAL-Pro2 | 13 | 0.0256 | 0.0299 | -0.0043 | 1/12/0 | 1 |
| ASTRID-DISCO vs ASTRAL-Pro2 | 13 | 0.0299 | 0.0299 | +0.0000 | 2/9/2 | 1 |
| ASTRID-DISCOR-it2 vs ASTRID-DISCO | 13 | 0.0214 | 0.0299 | -0.0085 | 2/10/1 | 0.75 |
| ASTRAL-DISCOR-it2 vs ASTRAL-DISCO | 13 | 0.0256 | 0.0385 | -0.0128 | 3/10/0 | 0.25 |
| ASTRID-DISCOR-oracle vs ASTRID-DISCO | 13 | 0.0214 | 0.0299 | -0.0085 | 2/10/1 | 0.75 |
| ASTRAL-DISCOR-oracle vs ASTRAL-DISCO | 13 | 0.0256 | 0.0385 | -0.0128 | 3/10/0 | 0.25 |
| ASTRID-DISCOR-lca vs ASTRID-DISCO | 13 | 0.0299 | 0.0299 | +0.0000 | 2/9/2 | 1 |
| ASTRAL-DISCOR-lca vs ASTRAL-DISCO | 13 | 0.0299 | 0.0385 | -0.0085 | 2/11/0 | 0.5 |

**Per condition, g_50** (diff | W/T/L | p)

| condition | AD-DISCOR vs AD-DISCO | AL-DISCOR vs AL-DISCO | AD-DISCOR vs Pro2 | AD-oracle vs AD-DISCO |
|---|---|---|---|---|
| 20_gdl_1e-10_0 | +0.0000 | 0/10/0 | 1 | -0.0056 | 1/9/0 | 1 | -0.0056 | 2/7/1 | 1 | +0.0000 | 0/10/0 | 1 |
| 20_gdl_1e-10_0_hILS | -0.0417 | 2/2/0 | 0.5 | -0.0278 | 2/1/1 | 0.75 | -0.0139 | 1/3/0 | 1 | -0.0417 | 2/2/0 | 0.5 |
| 20_gdl_5e-10_1 | -0.0139 | 2/5/1 | 0.75 | +0.0000 | 0/8/0 | 1 | -0.0208 | 2/6/0 | 0.5 | -0.0139 | 2/5/1 | 0.75 |
| 20_gdl_5e-10_1_hILS | +0.0000 | 1/8/1 | 1 | -0.0167 | 4/5/1 | 0.375 | -0.0056 | 2/6/2 | 1 | +0.0000 | 1/8/1 | 1 |

### Mean species-tree RF, gene trees = g_true

| condition | n | AL-Pro2 | AD-DISCO | AL-DISCO | AD-DISCOR | AL-DISCOR | AD-DISCOR-it2 | AL-DISCOR-it2 | AD-DISCOR-lca | AL-DISCOR-lca | AD-DISCOR-oracle | AL-DISCOR-oracle |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20_gdl_1e-10_0 | 10 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 20_gdl_1e-10_0_hILS | 4 | 0.0278 | 0.0139 | 0.0278 | 0.0139 | 0.0278 | 0.0139 | 0.0278 | 0.0278 | 0.0278 | 0.0139 | 0.0278 |
| 20_gdl_5e-10_1 | 5 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 20_gdl_5e-10_1_hILS | 5 | 0.0556 | 0.0111 | 0.0556 | 0.0222 | 0.0333 | 0.0222 | 0.0333 | 0.0556 | 0.0556 | 0.0222 | 0.0333 |

**Paired comparisons, g_true, all replicates** (diff = mean RF(A) − RF(B); W/T/L = A better/tie/worse, tie = identical RF; two-sided Wilcoxon signed-rank on non-tied pairs)

| A vs B | n | mean A | mean B | diff | W/T/L | p |
|---|---|---|---|---|---|---|
| ASTRID-DISCOR vs ASTRID-DISCO | 24 | 0.0069 | 0.0046 | +0.0023 | 0/23/1 | 1 |
| ASTRAL-DISCOR vs ASTRAL-DISCO | 24 | 0.0116 | 0.0162 | -0.0046 | 2/22/0 | 0.5 |
| ASTRID-DISCOR vs ASTRAL-Pro2 | 24 | 0.0069 | 0.0162 | -0.0093 | 4/19/1 | 0.312 |
| ASTRAL-DISCOR vs ASTRAL-Pro2 | 24 | 0.0116 | 0.0162 | -0.0046 | 2/22/0 | 0.5 |
| ASTRID-DISCO vs ASTRAL-Pro2 | 24 | 0.0046 | 0.0162 | -0.0116 | 4/20/0 | 0.125 |
| ASTRID-DISCOR-it2 vs ASTRID-DISCO | 24 | 0.0069 | 0.0046 | +0.0023 | 0/23/1 | 1 |
| ASTRAL-DISCOR-it2 vs ASTRAL-DISCO | 24 | 0.0116 | 0.0162 | -0.0046 | 2/22/0 | 0.5 |
| ASTRID-DISCOR-oracle vs ASTRID-DISCO | 24 | 0.0069 | 0.0046 | +0.0023 | 0/23/1 | 1 |
| ASTRAL-DISCOR-oracle vs ASTRAL-DISCO | 24 | 0.0116 | 0.0162 | -0.0046 | 2/22/0 | 0.5 |
| ASTRID-DISCOR-lca vs ASTRID-DISCO | 24 | 0.0162 | 0.0046 | +0.0116 | 0/20/4 | 0.125 |
| ASTRAL-DISCOR-lca vs ASTRAL-DISCO | 24 | 0.0162 | 0.0162 | +0.0000 | 0/24/0 | 1 |

**Paired comparisons, g_true, held-out replicates 06,07,08,09,10** (diff = mean RF(A) − RF(B); W/T/L = A better/tie/worse, tie = identical RF; two-sided Wilcoxon signed-rank on non-tied pairs)

| A vs B | n | mean A | mean B | diff | W/T/L | p |
|---|---|---|---|---|---|---|
| ASTRID-DISCOR vs ASTRID-DISCO | 5 | 0.0000 | 0.0000 | +0.0000 | 0/5/0 | 1 |
| ASTRAL-DISCOR vs ASTRAL-DISCO | 5 | 0.0000 | 0.0000 | +0.0000 | 0/5/0 | 1 |
| ASTRID-DISCOR vs ASTRAL-Pro2 | 5 | 0.0000 | 0.0000 | +0.0000 | 0/5/0 | 1 |
| ASTRAL-DISCOR vs ASTRAL-Pro2 | 5 | 0.0000 | 0.0000 | +0.0000 | 0/5/0 | 1 |
| ASTRID-DISCO vs ASTRAL-Pro2 | 5 | 0.0000 | 0.0000 | +0.0000 | 0/5/0 | 1 |
| ASTRID-DISCOR-it2 vs ASTRID-DISCO | 5 | 0.0000 | 0.0000 | +0.0000 | 0/5/0 | 1 |
| ASTRAL-DISCOR-it2 vs ASTRAL-DISCO | 5 | 0.0000 | 0.0000 | +0.0000 | 0/5/0 | 1 |
| ASTRID-DISCOR-oracle vs ASTRID-DISCO | 5 | 0.0000 | 0.0000 | +0.0000 | 0/5/0 | 1 |
| ASTRAL-DISCOR-oracle vs ASTRAL-DISCO | 5 | 0.0000 | 0.0000 | +0.0000 | 0/5/0 | 1 |
| ASTRID-DISCOR-lca vs ASTRID-DISCO | 5 | 0.0000 | 0.0000 | +0.0000 | 0/5/0 | 1 |
| ASTRAL-DISCOR-lca vs ASTRAL-DISCO | 5 | 0.0000 | 0.0000 | +0.0000 | 0/5/0 | 1 |

**Per condition, g_true** (diff | W/T/L | p)

| condition | AD-DISCOR vs AD-DISCO | AL-DISCOR vs AL-DISCO | AD-DISCOR vs Pro2 | AD-oracle vs AD-DISCO |
|---|---|---|---|---|
| 20_gdl_1e-10_0 | +0.0000 | 0/10/0 | 1 | +0.0000 | 0/10/0 | 1 | +0.0000 | 0/10/0 | 1 | +0.0000 | 0/10/0 | 1 |
| 20_gdl_1e-10_0_hILS | +0.0000 | 0/4/0 | 1 | +0.0000 | 0/4/0 | 1 | -0.0139 | 1/3/0 | 1 | +0.0000 | 0/4/0 | 1 |
| 20_gdl_5e-10_1 | +0.0000 | 0/5/0 | 1 | +0.0000 | 0/5/0 | 1 | +0.0000 | 0/5/0 | 1 | +0.0000 | 0/5/0 | 1 |
| 20_gdl_5e-10_1_hILS | +0.0111 | 0/4/1 | 1 | -0.0222 | 2/3/0 | 0.5 | -0.0333 | 3/1/1 | 0.5 | +0.0111 | 0/4/1 | 1 |

### Tagging accuracy (pooled over replicates; pair acc = fraction of cross-species gene pairs whose ortholog/paralog call matches the locus tree)

| condition | genes | DISCO pair acc | DISCO-R pair acc | DISCO-R-lca | oracle-S | DISCO orth P/R | DISCO-R orth P/R | DISCO root acc* | DISCO-R root acc* |
|---|---|---|---|---|---|---|---|---|---|
| 20_gdl_1e-10_0 | g_100 | 0.912 | 0.920 | 0.718 | 0.920 | 0.983/0.881 | 0.988/0.889 | – | – |
| 20_gdl_1e-10_0 | g_50 | 0.879 | 0.885 | 0.638 | 0.885 | 0.973/0.840 | 0.980/0.842 | – | – |
| 20_gdl_1e-10_0 | g_true | 0.982 | 0.994 | 0.925 | 0.994 | 1.000/0.973 | 1.000/0.991 | 0.116 | 0.992 |
| 20_gdl_1e-10_0_hILS | g_100 | 0.872 | 0.883 | 0.605 | 0.883 | 0.984/0.821 | 0.989/0.835 | – | – |
| 20_gdl_1e-10_0_hILS | g_50 | 0.843 | 0.849 | 0.556 | 0.849 | 0.975/0.780 | 0.982/0.785 | – | – |
| 20_gdl_1e-10_0_hILS | g_true | 0.937 | 0.950 | 0.740 | 0.950 | 1.000/0.904 | 1.000/0.924 | 0.091 | 0.900 |
| 20_gdl_5e-10_1 | g_100 | 0.867 | 0.884 | 0.820 | 0.884 | 0.812/0.720 | 0.814/0.791 | – | – |
| 20_gdl_5e-10_1 | g_50 | 0.843 | 0.858 | 0.786 | 0.858 | 0.797/0.631 | 0.803/0.693 | – | – |
| 20_gdl_5e-10_1 | g_true | 0.925 | 0.946 | 0.914 | 0.946 | 0.840/0.907 | 0.857/0.967 | 0.341 | 0.687 |
| 20_gdl_5e-10_1_hILS | g_100 | 0.829 | 0.846 | 0.783 | 0.846 | 0.787/0.557 | 0.788/0.638 | – | – |
| 20_gdl_5e-10_1_hILS | g_50 | 0.809 | 0.822 | 0.763 | 0.822 | 0.771/0.476 | 0.772/0.546 | – | – |
| 20_gdl_5e-10_1_hILS | g_true | 0.873 | 0.893 | 0.830 | 0.893 | 0.855/0.705 | 0.865/0.775 | 0.238 | 0.596 |

\* root acc: fraction of multi-copy true gene trees (input rooting randomised) whose chosen root bipartition equals SimPhy's.

### Runtime (mean seconds per replicate, 1 thread)

| condition | genes | astral-pro2 | DISCO-decomp | DISCO-astrid | DISCO-astral | s0-root | DISCOR-decomp | DISCOR-astrid | DISCOR-astral |
|---|---|---|---|---|---|---|---|---|---|
| 20_gdl_1e-10_0 | g_100 | 6.0 | 1.7 | 0.0 | 9.0 | 5.1 | 1.5 | 0.0 | 8.1 |
| 20_gdl_1e-10_0 | g_50 | 6.4 | 1.8 | 0.0 | 9.2 | 5.3 | 1.4 | 0.0 | 8.3 |
| 20_gdl_1e-10_0 | g_true | 6.2 | 1.7 | 0.0 | 9.5 | 5.1 | 1.6 | 0.0 | 8.8 |
| 20_gdl_1e-10_0_hILS | g_100 | 6.1 | 1.4 | 0.0 | 8.3 | 5.3 | 1.8 | 0.0 | 8.8 |
| 20_gdl_1e-10_0_hILS | g_50 | 6.3 | 1.7 | 0.0 | 8.3 | 5.1 | 1.3 | 0.0 | 7.9 |
| 20_gdl_1e-10_0_hILS | g_true | 5.9 | 1.6 | 0.0 | 8.4 | 4.8 | 1.4 | 0.0 | 8.1 |
| 20_gdl_5e-10_1 | g_100 | 4.5 | 1.4 | 0.0 | 5.2 | 4.2 | 1.3 | 0.0 | 5.1 |
| 20_gdl_5e-10_1 | g_50 | 4.8 | 1.3 | 0.0 | 5.0 | 4.3 | 1.3 | 0.0 | 4.7 |
| 20_gdl_5e-10_1 | g_true | 4.2 | 1.4 | 0.0 | 5.7 | 4.3 | 1.4 | 0.0 | 6.0 |
| 20_gdl_5e-10_1_hILS | g_100 | 5.9 | 1.8 | 0.0 | 6.0 | 5.8 | 1.7 | 0.0 | 5.2 |
| 20_gdl_5e-10_1_hILS | g_50 | 6.3 | 1.7 | 0.0 | 5.7 | 5.8 | 1.8 | 0.0 | 5.0 |
| 20_gdl_5e-10_1_hILS | g_true | 6.5 | 2.3 | 0.0 | 7.6 | 6.3 | 2.2 | 0.0 | 6.6 |

S0 root correct (min-DL rooting of ASTRAL-Pro2 tree): 88/88 runs.
