Tie band for tree FN: |diff| <= 0.10 points. Negative diff = variant better than reference.

### Tree FN % (FastTree GTR+G), 21 replicates

Mean FN: MAGUS 9.53, MAGUS(Slow) 9.62, slow-soft-m3 9.53, TRUE 8.90

| comparison | n | mean diff (FN pts) | W/T/L | Wilcoxon p (two-sided) |
|---|---|---|---|---|
| slow-soft-m3 − MAGUS | 21 | -0.00 | 9/2/10 | 0.904 |
| slow-soft-m3 − MAGUS(Slow) | 21 | -0.09 | 12/2/7 | 0.365 |
| MAGUS(Slow) − MAGUS | 21 | +0.09 | 7/3/11 | 0.42 |
| MAGUS − TRUE | 21 | +0.63 | 5/2/14 | 0.0141 |
| slow-soft-m3 − TRUE | 21 | +0.62 | 5/1/15 | 0.0228 |

| replicate | MAGUS | MAGUS(Slow) | slow-soft-m3 | TRUE |
|---|---|---|---|---|
| 1000L1_R0 | 11.58 | 11.28 | 11.38 | 9.16 |
| 1000L1_R1 | 10.55 | 10.85 | 12.16 | 9.55 |
| 1000L1_R2 | 11.29 | 10.58 | 9.88 | 9.27 |
| 1000L2_R0 | 11.87 | 12.27 | 12.98 | 11.27 |
| 1000L2_R1 | 8.35 | 8.25 | 9.05 | 8.25 |
| 1000L3_R0 | 10.56 | 11.67 | 11.57 | 8.95 |
| 1000L3_R3 | 10.84 | 10.94 | 9.64 | 9.54 |
| 1000M2_R0 | 9.68 | 10.58 | 11.29 | 10.18 |
| 1000M3_R0 | 7.55 | 7.45 | 7.65 | 7.85 |
| 1000M3_R1 | 6.24 | 6.94 | 6.34 | 6.84 |
| 1000M3_R3 | 7.26 | 7.46 | 7.26 | 7.06 |
| 1000M4_R0 | 6.29 | 6.29 | 6.19 | 5.98 |
| 1000M4_R1 | 3.79 | 3.79 | 3.89 | 3.79 |
| 1000M4_R4 | 3.82 | 4.03 | 4.03 | 4.34 |
| 1000S1_R0 | 11.78 | 12.29 | 11.68 | 10.27 |
| 1000S1_R1 | 12.89 | 12.79 | 11.98 | 11.68 |
| 1000S1_R2 | 9.77 | 9.16 | 8.96 | 7.55 |
| 1000S2_R0 | 8.04 | 8.14 | 8.14 | 9.05 |
| 1000S2_R1 | 11.78 | 10.67 | 10.57 | 11.18 |
| 1000S3_R0 | 11.90 | 12.20 | 11.90 | 11.90 |
| RNASim_R0 | 14.34 | 14.34 | 13.54 | 13.34 |

### Tree FN % (IQ-TREE --fast GTR+G4), 13 replicates

Mean FN: MAGUS 8.95, MAGUS(Slow) 9.10, slow-soft-m3 8.98, TRUE 8.74

| comparison | n | mean diff (FN pts) | W/T/L | Wilcoxon p (two-sided) |
|---|---|---|---|---|
| slow-soft-m3 − MAGUS | 13 | +0.02 | 5/4/4 | 1 |
| slow-soft-m3 − MAGUS(Slow) | 13 | -0.12 | 4/5/4 | 0.641 |
| MAGUS(Slow) − MAGUS | 13 | +0.15 | 5/4/4 | 0.652 |
| MAGUS − TRUE | 13 | +0.21 | 6/0/7 | 0.787 |
| slow-soft-m3 − TRUE | 13 | +0.23 | 7/1/5 | 0.834 |

| replicate | MAGUS | MAGUS(Slow) | slow-soft-m3 | TRUE |
|---|---|---|---|---|
| 1000L1_R0 | 11.48 | 11.48 | 11.48 | 10.67 |
| 1000L1_R1 | 9.65 | 11.36 | 10.15 | 10.35 |
| 1000L1_R2 | 11.49 | 11.49 | 11.49 | 10.99 |
| 1000L2_R0 | 11.87 | 11.77 | 11.77 | 11.97 |
| 1000L2_R1 | 9.26 | 9.86 | 8.75 | 8.95 |
| 1000L3_R0 | 11.37 | 11.47 | 11.07 | 8.05 |
| 1000L3_R3 | 9.34 | 9.04 | 9.34 | 9.24 |
| 1000M2_R0 | 9.38 | 9.27 | 9.88 | 10.08 |
| 1000M3_R0 | 6.95 | 7.25 | 7.96 | 8.26 |
| 1000M3_R1 | 6.54 | 6.44 | 5.84 | 5.84 |
| 1000M3_R3 | 8.67 | 8.47 | 8.47 | 8.37 |
| 1000M4_R0 | 6.49 | 6.49 | 6.60 | 6.70 |
| 1000M4_R1 | 3.89 | 3.89 | 3.89 | 4.20 |

### Alignment criteria on the same replicates (means)

| variant | n | SPFN % | SPFP % | (SPFN+SPFP)/2 % | compression (est/true length) | merge seconds |
|---|---|---|---|---|---|---|
| MAGUS | 21 | 6.12 | 5.55 | 5.84 | 0.997 | 43 |
| MAGUS(Slow) | 21 | 5.61 | 5.35 | 5.48 | 0.997 | 64 |
| slow-soft-m3 | 21 | 5.25 | 4.91 | 5.08 | 1.185 | 200 |

slow-soft-m3 − MAGUS alignment error: -0.75 pts, W/T/L 19/0/2 (tie 0.05), Wilcoxon p = 2.38e-05

#### ft: within-replicate regressions of tree FN (incl. TRUE, 21 replicates, n=84)

| model | coefs | cluster-robust p | within R² |
|---|---|---|---|
| SPFN | 0.130 | 2.5e-08 | 0.393 |
| SPFP | 0.140 | 1.9e-07 | 0.397 |
| log_comp | 0.864 | 0.2 | 0.018 |
| TC_err | 0.009 | 0.0005 | 0.267 |
| SPFN + SPFP | -0.015, 0.156 | 0.95, 0.54 | 0.397 |
| SPFN + SPFP + log_comp | -0.013, 0.154, 0.016 | 0.95, 0.54, 0.98 | 0.397 |

#### ft: within-replicate regressions of tree FN (estimated only, 21 replicates, n=63)

| model | coefs | cluster-robust p | within R² |
|---|---|---|---|
| SPFN | 0.106 | 0.54 | 0.023 |
| SPFP | 0.310 | 0.022 | 0.122 |
| log_comp | -0.047 | 0.96 | 0.000 |
| TC_err | 0.015 | 0.78 | 0.003 |
| SPFN + SPFP | -0.617, 1.011 | 8.9e-10, 1.2e-09 | 0.265 |
| SPFN + SPFP + log_comp | -0.600, 1.213, 1.658 | 7.8e-15, 5.2e-09, 0.15 | 0.344 |

#### iqfast: within-replicate regressions of tree FN (incl. TRUE, 13 replicates, n=52)

| model | coefs | cluster-robust p | within R² |
|---|---|---|---|
| SPFN | 0.084 | 0.23 | 0.174 |
| SPFP | 0.088 | 0.25 | 0.157 |
| log_comp | -0.737 | 0.38 | 0.013 |
| TC_err | 0.005 | 0.31 | 0.074 |
| SPFN + SPFP | 0.382, -0.334 | 1.2e-05, 0.02 | 0.203 |
| SPFN + SPFP + log_comp | 0.363, -0.313, -0.479 | 0.00091, 0.066, 0.44 | 0.208 |

#### iqfast: within-replicate regressions of tree FN (estimated only, 13 replicates, n=39)

| model | coefs | cluster-robust p | within R² |
|---|---|---|---|
| SPFN | -0.022 | 0.81 | 0.001 |
| SPFP | -0.046 | 0.77 | 0.003 |
| log_comp | -0.312 | 0.62 | 0.008 |
| TC_err | 0.006 | 0.87 | 0.001 |
| SPFN + SPFP | 0.009, -0.057 | 0.76, 0.72 | 0.003 |
| SPFN + SPFP + log_comp | -0.002, -0.121, -0.565 | 0.97, 0.5, 0.49 | 0.022 |

