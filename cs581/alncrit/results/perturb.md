Replicates with all alignment types: 21

| alignment | SPFN % | SPFP % | compression | tree FN % | FN excess over TRUE (pts) | Wilcoxon p (excess ≠ 0) | excess per SPFN point |
|---|---|---|---|---|---|---|---|
| TRUE | 0.00 | 0.000 | 1.000 | 8.90 | +0.00 | nan | nan |
| split_0.1 | 4.89 | 0.000 | 1.054 | 8.97 | +0.06 | 0.472 | 0.014 |
| shift_0.1 | 8.39 | 0.187 | 0.997 | 9.25 | +0.35 | 0.136 | 0.066 |
| split_0.2 | 9.75 | 0.000 | 1.107 | 9.40 | +0.49 | 0.0421 | 0.051 |
| MAGUS | 6.12 | 5.547 | 0.997 | 9.53 | +0.63 | 0.0141 | 0.119 |
| MAGUS(Slow) | 5.61 | 5.346 | 0.997 | 9.62 | +0.71 | 0.00549 | 0.144 |
| slow-soft-m3 | 5.25 | 4.911 | 1.185 | 9.53 | +0.62 | 0.0228 | 0.130 |

- split_0.1 − MAGUS: tree FN -0.57 pts, W/T/L 13/3/5 (tie 0.1), Wilcoxon p = 0.0198
- split_0.2 − MAGUS: tree FN -0.14 pts, W/T/L 9/2/10 (tie 0.1), Wilcoxon p = 0.673
- shift_0.1 − MAGUS: tree FN -0.28 pts, W/T/L 12/2/7 (tie 0.1), Wilcoxon p = 0.243
- split_0.1 − slow-soft-m3: tree FN -0.56 pts, W/T/L 15/1/5 (tie 0.1), Wilcoxon p = 0.0304
