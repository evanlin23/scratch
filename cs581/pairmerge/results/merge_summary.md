
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

### Condition: linsi200

Mean (SPFN+SPFP)/2 in %, full alignment; n = replicates

| merger | 1000L2 | 1000M2 | 1000M4 | 1000S2 | RNASim1000 | all |
|---|---|---|---|---|---|---|
| opal | 8.31 (n=2) | 24.05 (n=2) | 1.18 (n=2) | 10.50 (n=2) | 10.46 (n=1) | 10.95 (n=9) |
| muscle3 | 9.50 (n=2) | 30.35 (n=2) | 1.61 (n=2) | 13.23 (n=2) | 30.00 (n=1) | 15.49 (n=9) |
| mafft-merge | 10.43 (n=2) | 26.75 (n=2) | 2.21 (n=2) | 13.97 (n=2) | 13.79 (n=1) | 13.39 (n=9) |
| gcm | 8.07 (n=2) | 20.92 (n=2) | 1.25 (n=2) | 10.45 (n=2) | 10.27 (n=1) | 10.18 (n=9) |
| progdp | 8.07 (n=2) | 20.89 (n=2) | 1.20 (n=2) | 10.44 (n=2) | 10.21 (n=1) | 10.16 (n=9) |

All replicates pooled, in %: full-alignment SPFN/SPFP, cross-half pairs only, runtime

| merger | SPFN | SPFP | cross FN | cross FP | constraints kept | mean seconds |
|---|---|---|---|---|---|---|
| opal | 11.73 | 10.17 | 13.77 | 12.11 | 9/9 | 1.3 |
| muscle3 | 17.20 | 13.77 | 24.72 | 20.19 | 9/9 | 0.7 |
| mafft-merge | 15.46 | 11.32 | 21.22 | 14.71 | 9/9 | 1.4 |
| gcm | 11.24 | 9.13 | 12.78 | 10.04 | 9/9 | 1.6 |
| progdp | 10.99 | 9.32 | 12.29 | 10.41 | 9/9 | 1.5 |

Paired vs `opal` on avgErr (negative = better than opal; W/T/L = better/tie/worse, tie band 0.1 points)

| merger | n | mean diff (points) | W/T/L | Wilcoxon p |
|---|---|---|---|---|
| muscle3 | 9 | +4.54 | 0/0/9 | 0.0039 |
| mafft-merge | 9 | +2.44 | 0/0/9 | 0.0039 |
| gcm | 9 | -0.77 | 3/6/0 | 0.3 |
| progdp | 9 | -0.79 | 5/4/0 | 0.074 |

progdp vs gcm (same graph): n=9, mean diff -0.03 points, W/T/L 1/8/0, p=0.13

### Condition: oracle200

Mean (SPFN+SPFP)/2 in %, full alignment; n = replicates

| merger | 1000L2 | 1000M2 | 1000M4 | 1000S2 | RNASim1000 | all |
|---|---|---|---|---|---|---|
| opal | 0.20 (n=2) | 0.28 (n=2) | 0.05 (n=2) | 0.38 (n=2) | 1.78 (n=1) | 0.40 (n=9) |
| muscle3 | 2.84 (n=2) | 14.73 (n=2) | 1.33 (n=2) | 5.93 (n=2) | 26.97 (n=1) | 8.51 (n=9) |
| mafft-merge | 2.86 (n=2) | 9.30 (n=2) | 1.30 (n=2) | 6.79 (n=2) | 8.91 (n=1) | 5.49 (n=9) |
| gcm | 0.60 (n=2) | 1.00 (n=2) | 0.15 (n=2) | 0.54 (n=2) | 2.57 (n=1) | 0.79 (n=9) |
| progdp | 0.55 (n=2) | 0.81 (n=2) | 0.07 (n=2) | 0.38 (n=2) | 2.18 (n=1) | 0.64 (n=9) |

All replicates pooled, in %: full-alignment SPFN/SPFP, cross-half pairs only, runtime

| merger | SPFN | SPFP | cross FN | cross FP | constraints kept | mean seconds |
|---|---|---|---|---|---|---|
| opal | 0.40 | 0.41 | 0.79 | 0.81 | 9/9 | 1.4 |
| muscle3 | 9.91 | 7.12 | 19.83 | 15.68 | 9/9 | 0.7 |
| mafft-merge | 6.84 | 4.15 | 13.67 | 8.65 | 9/9 | 1.4 |
| gcm | 0.98 | 0.61 | 1.95 | 1.23 | 9/9 | 1.5 |
| progdp | 0.66 | 0.63 | 1.32 | 1.26 | 9/9 | 1.6 |

Paired vs `opal` on avgErr (negative = better than opal; W/T/L = better/tie/worse, tie band 0.1 points)

| merger | n | mean diff (points) | W/T/L | Wilcoxon p |
|---|---|---|---|---|
| muscle3 | 9 | +8.11 | 0/0/9 | 0.0039 |
| mafft-merge | 9 | +5.09 | 0/0/9 | 0.0039 |
| gcm | 9 | +0.39 | 1/2/6 | 0.027 |
| progdp | 9 | +0.24 | 2/2/5 | 0.2 |

progdp vs gcm (same graph): n=9, mean diff -0.15 points, W/T/L 6/2/1, p=0.074
