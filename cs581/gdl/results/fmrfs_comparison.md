# FastMulRFS-data comparison (estimated gene trees)

Runs: 360 (dev 108, held-out 252). Metric: species-tree FN rate (= RF rate; all trees binary). W/T/L = ASTRID-Pro better / tie (equal FN count) / worse. Two-sided Wilcoxon signed-rank, paired by (condition, replicate, sequence length, number of genes).


## Dev replicates 01-03: mean FN rate by (sequence length, #genes)

| sqln | ngen | n | astrid-multi | astrid-pro | astrid-pro-min | ortho-allnodes | spec-allpairs | astrid-disco | astral-pro | fastmulrfs | astral-multi(pub) | mulrf(pub) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 25 | 25 | 18 | 0.233 | 0.261 | 0.263 | 0.260 | 0.242 | 0.271 | 0.338 | 0.324 | 0.391 | 0.373 |
| 25 | 100 | 18 | 0.108 | 0.107 | 0.100 | 0.103 | 0.107 | 0.101 | 0.129 | 0.151 | 0.234 | 0.159 |
| 25 | 500 | 18 | 0.052 | 0.044 | 0.045 | 0.041 | 0.049 | 0.045 | 0.053 | 0.066 | 0.121 | 0.064 |
| 100 | 25 | 18 | 0.095 | 0.082 | 0.080 | 0.083 | 0.089 | 0.081 | 0.092 | 0.117 | 0.137 | 0.114 |
| 100 | 100 | 18 | 0.048 | 0.044 | 0.045 | 0.043 | 0.047 | 0.045 | 0.052 | 0.056 | 0.078 | 0.054 |
| 100 | 500 | 18 | 0.020 | 0.019 | 0.017 | 0.022 | 0.020 | 0.021 | 0.022 | 0.022 | 0.038 | 0.022 |
| all | all | 108 | 0.0925 | 0.0925 | 0.0916 | 0.0920 | 0.0922 | 0.0939 | 0.1144 | 0.1227 | 0.1665 | 0.1311 |

Dev choice: mean ASTRID-Pro FN 0.0925 vs closest-copy 0.0916 -> **astrid-pro-min** used as ASTRID-Pro on held-out.


## Held-out replicates 04-10: mean FN rate by (sequence length, #genes)

| sqln | ngen | n | astrid-multi | astrid-pro | astrid-pro-min | ortho-allnodes | spec-allpairs | astrid-disco | astral-pro | fastmulrfs | astral-multi(pub) | mulrf(pub) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 25 | 25 | 42 | 0.229 | 0.244 | 0.252 | 0.243 | 0.235 | 0.254 | 0.284 | 0.326 | 0.357 | 0.359 |
| 25 | 100 | 42 | 0.105 | 0.106 | 0.109 | 0.104 | 0.110 | 0.103 | 0.130 | 0.153 | 0.227 | 0.163 |
| 25 | 500 | 42 | 0.050 | 0.047 | 0.046 | 0.049 | 0.051 | 0.046 | 0.058 | 0.073 | 0.137 | 0.070 |
| 100 | 25 | 42 | 0.101 | 0.094 | 0.092 | 0.092 | 0.097 | 0.089 | 0.100 | 0.121 | 0.143 | 0.116 |
| 100 | 100 | 42 | 0.045 | 0.041 | 0.043 | 0.041 | 0.043 | 0.041 | 0.047 | 0.057 | 0.076 | 0.051 |
| 100 | 500 | 42 | 0.022 | 0.019 | 0.020 | 0.020 | 0.021 | 0.017 | 0.019 | 0.024 | 0.036 | 0.019 |
| all | all | 252 | 0.0920 | 0.0919 | 0.0937 | 0.0914 | 0.0926 | 0.0916 | 0.1064 | 0.1256 | 0.1626 | 0.1298 |

## Held-out (all strata pooled): astrid-pro-min vs baselines

| baseline | n | mean FN (Pro) | mean FN (base) | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|
| astrid-multi | 252 | 0.0937 | 0.0920 | +0.0018 | 90/72/90 | 0.9 |
| astrid-disco | 252 | 0.0937 | 0.0916 | +0.0021 | 67/85/100 | 0.015 |
| astral-pro | 252 | 0.0937 | 0.1064 | -0.0127 | 131/45/76 | 3.6e-07 |
| fastmulrfs | 252 | 0.0937 | 0.1256 | -0.0319 | 174/33/45 | 1.9e-23 |
| astral-multi(pub) | 252 | 0.0937 | 0.1626 | -0.0689 | 221/22/9 | 8.4e-38 |
| mulrf(pub) | 252 | 0.0937 | 0.1298 | -0.0360 | 173/31/48 | 3.5e-23 |
| astrid-pro | 252 | 0.0937 | 0.0919 | +0.0018 | 54/126/72 | 0.1 |

## Held-out, 25 genes: astrid-pro-min vs baselines

| baseline | n | mean FN (Pro) | mean FN (base) | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|
| astrid-multi | 84 | 0.1721 | 0.1648 | +0.0072 | 33/13/38 | 0.35 |
| astrid-disco | 84 | 0.1721 | 0.1715 | +0.0006 | 30/16/38 | 0.45 |
| astral-pro | 84 | 0.1721 | 0.1924 | -0.0204 | 48/6/30 | 0.001 |
| fastmulrfs | 84 | 0.1721 | 0.2234 | -0.0513 | 66/4/14 | 5.5e-11 |
| astral-multi(pub) | 84 | 0.1721 | 0.2501 | -0.0781 | 76/4/4 | 3.3e-14 |
| mulrf(pub) | 84 | 0.1721 | 0.2372 | -0.0652 | 67/6/11 | 9.3e-12 |
| astrid-pro | 84 | 0.1721 | 0.1689 | +0.0032 | 23/31/30 | 0.31 |

## Held-out, 100 genes: astrid-pro-min vs baselines

| baseline | n | mean FN (Pro) | mean FN (base) | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|
| astrid-multi | 84 | 0.0761 | 0.0752 | +0.0009 | 28/24/32 | 0.48 |
| astrid-disco | 84 | 0.0761 | 0.0716 | +0.0045 | 22/22/40 | 0.006 |
| astral-pro | 84 | 0.0761 | 0.0884 | -0.0123 | 47/18/19 | 0.0014 |
| fastmulrfs | 84 | 0.0761 | 0.1052 | -0.0291 | 58/10/16 | 2.8e-09 |
| astral-multi(pub) | 84 | 0.0761 | 0.1512 | -0.0751 | 75/8/1 | 5.6e-14 |
| mulrf(pub) | 84 | 0.0761 | 0.1071 | -0.0311 | 61/7/16 | 5.7e-09 |
| astrid-pro | 84 | 0.0761 | 0.0735 | +0.0026 | 14/42/28 | 0.056 |

## Held-out, 500 genes: astrid-pro-min vs baselines

| baseline | n | mean FN (Pro) | mean FN (base) | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|
| astrid-multi | 84 | 0.0330 | 0.0358 | -0.0028 | 29/35/20 | 0.15 |
| astrid-disco | 84 | 0.0330 | 0.0319 | +0.0011 | 15/47/22 | 0.36 |
| astral-pro | 84 | 0.0330 | 0.0384 | -0.0054 | 36/21/27 | 0.026 |
| fastmulrfs | 84 | 0.0330 | 0.0484 | -0.0153 | 50/19/15 | 2.6e-06 |
| astral-multi(pub) | 84 | 0.0330 | 0.0864 | -0.0534 | 70/10/4 | 7.8e-13 |
| mulrf(pub) | 84 | 0.0330 | 0.0449 | -0.0119 | 45/18/21 | 6.8e-05 |
| astrid-pro | 84 | 0.0330 | 0.0333 | -0.0002 | 17/53/14 | 0.94 |

## Held-out, 25 bp: astrid-pro-min vs baselines

| baseline | n | mean FN (Pro) | mean FN (base) | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|
| astrid-multi | 126 | 0.1359 | 0.1280 | +0.0079 | 41/27/58 | 0.024 |
| astrid-disco | 126 | 0.1359 | 0.1343 | +0.0016 | 38/34/54 | 0.08 |
| astral-pro | 126 | 0.1359 | 0.1574 | -0.0215 | 75/17/34 | 1.1e-06 |
| fastmulrfs | 126 | 0.1359 | 0.1839 | -0.0480 | 100/12/14 | 7.7e-18 |
| astral-multi(pub) | 126 | 0.1359 | 0.2403 | -0.1044 | 126/0/0 | 2e-22 |
| mulrf(pub) | 126 | 0.1359 | 0.1975 | -0.0616 | 108/6/12 | 4.8e-19 |
| astrid-pro | 126 | 0.1359 | 0.1323 | +0.0036 | 28/54/44 | 0.024 |

## Held-out, 100 bp: astrid-pro-min vs baselines

| baseline | n | mean FN (Pro) | mean FN (base) | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|
| astrid-multi | 126 | 0.0515 | 0.0559 | -0.0043 | 49/45/32 | 0.0036 |
| astrid-disco | 126 | 0.0515 | 0.0490 | +0.0025 | 29/51/46 | 0.12 |
| astral-pro | 126 | 0.0515 | 0.0554 | -0.0038 | 56/28/42 | 0.061 |
| fastmulrfs | 126 | 0.0515 | 0.0673 | -0.0158 | 74/21/31 | 2e-06 |
| astral-multi(pub) | 126 | 0.0515 | 0.0848 | -0.0333 | 95/22/9 | 8.5e-16 |
| mulrf(pub) | 126 | 0.0515 | 0.0620 | -0.0105 | 65/25/36 | 0.00021 |
| astrid-pro | 126 | 0.0515 | 0.0515 | +0.0001 | 26/72/28 | 0.65 |

## Held-out, 250 bp: astrid-pro-min vs baselines

| baseline | n | mean FN (Pro) | mean FN (base) | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|

## Runtime (seconds per run, mean; 100 species; single core)

| step | 25 genes | 100 genes | 500 genes |
|---|---|---|---|
| astral-pro | 3.44 | 13.36 | 87.71 |
| astrid-disco | 0.80 | 3.13 | 15.28 |
| astrid-multi | 0.20 | 0.67 | 3.28 |
| astrid-pro | 0.18 | 0.60 | 2.77 |
| astrid-pro-min | 0.16 | 0.52 | 2.33 |
| ortho-allnodes | 0.18 | 0.59 | 2.81 |
| root+tag | 0.22 | 0.91 | 4.45 |
| spec-allpairs | 0.21 | 0.70 | 3.28 |
