
### Condition: oracle

Mean (SPFN+SPFP)/2 in %, full alignment; n = replicates

| merger | 1000L2 | 1000M2 | 1000M4 | 1000S2 | RNASim1000 | all |
|---|---|---|---|---|---|---|
| opal | 0.15 (n=2) | 0.28 (n=2) | 0.05 (n=2) | 0.38 (n=2) | 1.25 (n=2) | 0.42 (n=10) |
| muscle3 | 6.66 (n=2) | 19.07 (n=2) | 3.17 (n=2) | 7.65 (n=2) | 28.24 (n=2) | 12.96 (n=10) |
| mafft-merge | 4.88 (n=2) | 10.47 (n=2) | 1.93 (n=2) | 8.68 (n=2) | 10.10 (n=2) | 7.21 (n=10) |
| gcm | 0.63 (n=2) | 0.93 (n=2) | 0.19 (n=2) | 0.63 (n=2) | 2.43 (n=2) | 0.96 (n=10) |
| progdp | 0.57 (n=2) | 0.80 (n=2) | 0.07 (n=2) | 0.64 (n=2) | 2.03 (n=2) | 0.82 (n=10) |

All replicates pooled, in %: full-alignment SPFN/SPFP, cross-half pairs only, runtime

| merger | SPFN | SPFP | cross FN | cross FP | constraints kept | mean seconds |
|---|---|---|---|---|---|---|
| opal | 0.42 | 0.42 | 0.85 | 0.87 | 10/10 | 2.2 |
| muscle3 | 15.65 | 10.27 | 31.99 | 24.87 | 10/10 | 5.2 |
| mafft-merge | 9.36 | 5.07 | 19.05 | 11.00 | 10/10 | 5.7 |
| gcm | 1.14 | 0.79 | 2.32 | 1.62 | 10/10 | 1841.0 |
| progdp | 0.84 | 0.81 | 1.71 | 1.65 | 10/10 | 2.7 |

Paired vs `opal` on avgErr (negative = better than opal; W/T/L = better/tie/worse, tie band 0.1 points)

| merger | n | mean diff (points) | W/T/L | Wilcoxon p |
|---|---|---|---|---|
| muscle3 | 10 | +12.54 | 0/0/10 | 0.002 |
| mafft-merge | 10 | +6.79 | 0/0/10 | 0.002 |
| gcm | 10 | +0.54 | 1/1/8 | 0.02 |
| progdp | 10 | +0.40 | 1/2/7 | 0.027 |

progdp vs gcm (same graph): n=10, mean diff -0.14 points, W/T/L 5/5/0, p=0.084

### Condition: oracle200

Mean (SPFN+SPFP)/2 in %, full alignment; n = replicates

| merger | 1000L2 | 1000M2 | 1000M4 | 1000S2 | RNASim1000 | all |
|---|---|---|---|---|---|---|
| opal | 0.20 (n=2) | 0.28 (n=2) | 0.05 (n=2) | 0.38 (n=2) | 1.28 (n=2) | 0.44 (n=10) |
| muscle3 | 2.84 (n=2) | 14.73 (n=2) | 1.33 (n=2) | 5.93 (n=2) | 27.20 (n=2) | 10.40 (n=10) |
| mafft-merge | 2.86 (n=2) | 9.30 (n=2) | 1.30 (n=2) | 6.79 (n=2) | 9.10 (n=2) | 5.87 (n=10) |
| gcm | 0.60 (n=2) | 1.00 (n=2) | 0.15 (n=2) | 0.54 (n=2) | 2.32 (n=2) | 0.92 (n=10) |
| progdp | 0.55 (n=2) | 0.81 (n=2) | 0.07 (n=2) | 0.38 (n=2) | 2.06 (n=2) | 0.77 (n=10) |

All replicates pooled, in %: full-alignment SPFN/SPFP, cross-half pairs only, runtime

| merger | SPFN | SPFP | cross FN | cross FP | constraints kept | mean seconds |
|---|---|---|---|---|---|---|
| opal | 0.44 | 0.44 | 0.87 | 0.89 | 10/10 | 1.4 |
| muscle3 | 12.32 | 8.49 | 24.68 | 19.34 | 10/10 | 0.7 |
| mafft-merge | 7.38 | 4.36 | 14.77 | 9.14 | 10/10 | 1.4 |
| gcm | 1.11 | 0.73 | 2.22 | 1.47 | 10/10 | 1.5 |
| progdp | 0.79 | 0.76 | 1.58 | 1.52 | 10/10 | 1.6 |

Paired vs `opal` on avgErr (negative = better than opal; W/T/L = better/tie/worse, tie band 0.1 points)

| merger | n | mean diff (points) | W/T/L | Wilcoxon p |
|---|---|---|---|---|
| muscle3 | 10 | +9.96 | 0/0/10 | 0.002 |
| mafft-merge | 10 | +5.43 | 0/0/10 | 0.002 |
| gcm | 10 | +0.48 | 1/2/7 | 0.014 |
| progdp | 10 | +0.33 | 2/2/6 | 0.11 |

progdp vs gcm (same graph): n=10, mean diff -0.15 points, W/T/L 7/2/1, p=0.037

### Condition: linsi200

Mean (SPFN+SPFP)/2 in %, full alignment; n = replicates

| merger | 1000L2 | 1000M2 | 1000M4 | 1000S2 | RNASim1000 | all |
|---|---|---|---|---|---|---|
| opal | 8.31 (n=2) | 24.05 (n=2) | 1.18 (n=2) | 10.50 (n=2) | 10.21 (n=2) | 10.85 (n=10) |
| muscle3 | 9.50 (n=2) | 30.35 (n=2) | 1.61 (n=2) | 13.23 (n=2) | 29.51 (n=2) | 16.84 (n=10) |
| mafft-merge | 10.43 (n=2) | 26.75 (n=2) | 2.21 (n=2) | 13.97 (n=2) | 13.99 (n=2) | 13.47 (n=10) |
| gcm | 8.07 (n=2) | 20.92 (n=2) | 1.25 (n=2) | 10.45 (n=2) | 9.92 (n=2) | 10.12 (n=10) |
| progdp | 8.07 (n=2) | 20.89 (n=2) | 1.20 (n=2) | 10.44 (n=2) | 9.88 (n=2) | 10.10 (n=10) |

All replicates pooled, in %: full-alignment SPFN/SPFP, cross-half pairs only, runtime

| merger | SPFN | SPFP | cross FN | cross FP | constraints kept | mean seconds |
|---|---|---|---|---|---|---|
| opal | 11.56 | 10.14 | 13.60 | 12.10 | 10/10 | 1.3 |
| muscle3 | 18.85 | 14.83 | 28.20 | 22.82 | 10/10 | 0.7 |
| mafft-merge | 15.54 | 11.41 | 21.56 | 14.93 | 10/10 | 1.4 |
| gcm | 11.09 | 9.16 | 12.66 | 10.14 | 10/10 | 1.6 |
| progdp | 10.86 | 9.34 | 12.20 | 10.48 | 10/10 | 1.5 |

Paired vs `opal` on avgErr (negative = better than opal; W/T/L = better/tie/worse, tie band 0.1 points)

| merger | n | mean diff (points) | W/T/L | Wilcoxon p |
|---|---|---|---|---|
| muscle3 | 10 | +5.99 | 0/0/10 | 0.002 |
| mafft-merge | 10 | +2.62 | 0/0/10 | 0.002 |
| gcm | 10 | -0.73 | 4/6/0 | 0.16 |
| progdp | 10 | -0.75 | 6/4/0 | 0.037 |

progdp vs gcm (same graph): n=10, mean diff -0.03 points, W/T/L 1/9/0, p=0.084

### Condition: fftnsi

Mean (SPFN+SPFP)/2 in %, full alignment; n = replicates

| merger | 1000L2 | 1000M2 | 1000M4 | 1000S2 | RNASim1000 | all |
|---|---|---|---|---|---|---|
| opal | 90.52 (n=2) | 89.37 (n=2) | 2.40 (n=2) | 84.99 (n=2) | 20.03 (n=2) | 57.46 (n=10) |
| muscle3 | 89.97 (n=2) | 90.01 (n=2) | 4.29 (n=2) | 84.48 (n=2) | 42.05 (n=2) | 62.16 (n=10) |
| mafft-merge | 89.44 (n=2) | 89.63 (n=2) | 4.17 (n=2) | 84.55 (n=2) | 22.02 (n=2) | 57.96 (n=10) |
| gcm | 86.19 (n=2) | 85.96 (n=2) | 2.48 (n=2) | 81.27 (n=2) | 18.62 (n=2) | 54.90 (n=10) |
| progdp | 88.84 (n=2) | 87.40 (n=2) | 2.40 (n=2) | 83.28 (n=2) | 18.48 (n=2) | 56.08 (n=10) |

All replicates pooled, in %: full-alignment SPFN/SPFP, cross-half pairs only, runtime

| merger | SPFN | SPFP | cross FN | cross FP | constraints kept | mean seconds |
|---|---|---|---|---|---|---|
| opal | 62.03 | 52.89 | 66.07 | 58.14 | 10/10 | 5.7 |
| muscle3 | 66.43 | 57.89 | 75.16 | 68.11 | 10/10 | 8.6 |
| mafft-merge | 61.84 | 54.09 | 65.71 | 57.29 | 10/10 | 7.3 |
| gcm | 60.06 | 49.74 | 62.08 | 39.59 | 10/10 | 61.4 |
| progdp | 59.65 | 52.51 | 61.24 | 53.52 | 10/10 | 16.3 |

Paired vs `opal` on avgErr (negative = better than opal; W/T/L = better/tie/worse, tie band 0.1 points)

| merger | n | mean diff (points) | W/T/L | Wilcoxon p |
|---|---|---|---|---|
| muscle3 | 10 | +4.70 | 4/0/6 | 0.16 |
| mafft-merge | 10 | +0.50 | 3/1/6 | 0.32 |
| gcm | 10 | -2.56 | 8/1/1 | 0.0098 |
| progdp | 10 | -1.38 | 7/2/1 | 0.02 |

progdp vs gcm (same graph): n=10, mean diff +1.18 points, W/T/L 3/1/6, p=0.084
