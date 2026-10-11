# DISCO-data comparison (doi:10.13012/B2IDB-4050038_V1), estimated gene trees from 100 bp

Mean species-tree FN rate (100 species). `-w` = support-weighted. Conditions: default (dup 5e-10, loss/dup 1, 1000 genes), gdl_1e-9_1 (dup 1e-9, loss/dup 1, 1000 genes), gdl_1e-9_05 (dup 1e-9, loss/dup 0.5, i.e. supercritical; only the first 100 genes because families average ~1000 leaves).

| condition | genes | n | astrid-multi | astrid-multi-w | astrid-pro | astrid-pro-min | astrid-pro-w | ortho-allnodes | astrid-disco | astral-pro |
|---|---|---|---|---|---|---|---|---|---|---|
| default | 1000 | 10 | 0.054 | 0.050 | 0.051 | 0.050 | 0.051 | 0.054 | 0.047 | 0.051 |
| gdl_1e-9_1 | 1000 | 10 | 0.037 | 0.038 | 0.027 | 0.023 | 0.026 | 0.031 | 0.029 | 0.029 |
| gdl_1e-9_05 | 100 | 10 | 0.089 | 0.086 | 0.046 | 0.063 | 0.044 | 0.051 | 0.042 | 0.051 |

## Paired tests, all conditions pooled (W/T/L from the first method's view; tie = equal FN)

| method | baseline | n | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| astrid-pro-min | astrid-multi | 30 | -0.0143 | 16/10/4 | 0.0018 |
| astrid-pro-min | astrid-multi-w | 30 | -0.0122 | 17/7/6 | 0.015 |
| astrid-pro-min | astrid-disco | 30 | +0.0065 | 7/9/14 | 0.087 |
| astrid-pro-min | astral-pro | 30 | +0.0020 | 12/7/11 | 0.84 |
| astrid-pro | astrid-multi | 30 | -0.0187 | 21/6/3 | 0.00022 |
| astrid-pro | astrid-multi-w | 30 | -0.0167 | 20/4/6 | 0.0027 |
| astrid-pro | astrid-disco | 30 | +0.0020 | 7/10/13 | 0.61 |
| astrid-pro | astral-pro | 30 | -0.0024 | 11/11/8 | 0.43 |
| astrid-pro-w | astrid-multi | 30 | -0.0197 | 19/10/1 | 0.00013 |
| astrid-pro-w | astrid-multi-w | 30 | -0.0177 | 17/10/3 | 0.00071 |
| astrid-pro-w | astrid-disco | 30 | +0.0010 | 7/13/10 | 0.98 |
| astrid-pro-w | astral-pro | 30 | -0.0034 | 14/8/8 | 0.22 |

## Runtime (s, mean per run, single core)

| condition | root+tag | astrid-multi | astrid-pro | astrid-disco | astral-pro |
|---|---|---|---|---|---|
| default | 8.0 | 8.2 | 4.8 | 28.0 | 258.8 |
| gdl_1e-9_1 | 14.5 | 14.7 | 6.2 | 38.0 | 330.3 |
| gdl_1e-9_05 | 11.9 | 68.6 | 2.7 | 22.1 | 103.1 |

## True gene trees with GDL + ILS (gtrees_10000_l1): FN rate vs number of genes

| genes | n reps | astrid-multi | astrid-pro | astrid-pro-min | astrid-disco | astral-pro |
|---|---|---|---|---|---|---|
| 100 | 5 | 0.037 | 0.020 | 0.022 | 0.016 | 0.020 |
| 1000 | 5 | 0.010 | 0.000 | 0.000 | 0.000 | 0.002 |
| 10000 | 1 | 0.010 | 0.000 | 0.000 | 0.000 | 0.000 |
