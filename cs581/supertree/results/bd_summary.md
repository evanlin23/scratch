### Mean RF error (%) / mean FN (%) / mean wall (s) / mean peak RSS (MB), n replicates

| condition | method | n | RF % | FN % | FP % | wall s | peak MB |
|---|---|---|---|---|---|---|---|
| bd500 | astral3 | 5 | 5.71 | 5.71 | 5.71 | 10.9 | 526 |
| bd500 | mrlft | 5 | 6.84 | 6.84 | 6.84 | 3.9 | 7 |
| bd500 | scs | 5 | 1.83 | 1.85 | 1.81 | 4.8 | 175 |
| bd500 | dc:astral3:100:astral3 | 5 | 7.79 | 7.85 | 7.74 | 6.8 | 114 |
| bd500 | dc:scs:100:astral3 | 5 | 5.57 | 5.63 | 5.52 | 7.4 | 111 |
| bd500 | dc:scs:500:astral3 | 5 | 5.71 | 5.71 | 5.71 | 9.8 | 465 |
| bd1000 | astral3 | 5 | 8.08 | 8.08 | 8.08 | 36.0 | 1640 |
| bd1000 | mrlft | 5 | 8.30 | 8.30 | 8.30 | 14.5 | 18 |
| bd1000 | scs | 5 | 2.80 | 2.89 | 2.71 | 9.2 | 218 |
| bd1000 | dc:astral3:100:astral3 | 5 | 11.03 | 11.19 | 10.91 | 16.4 | 118 |
| bd1000 | dc:scs:100:astral3 | 5 | 8.11 | 8.18 | 8.06 | 16.1 | 116 |
| bd1000 | dc:scs:500:astral3 | 5 | 8.05 | 8.06 | 8.05 | 19.5 | 451 |
| bd2000 | astral3 | 5 | 7.23 | 7.23 | 7.23 | 184.8 | 3699 |
| bd2000 | mrlft | 5 | 7.31 | 7.31 | 7.31 | 43.0 | 62 |
| bd2000 | scs | 5 | 3.20 | 3.32 | 3.08 | 14.9 | 286 |
| bd2000 | dc:astral3:100:astral3 | 5 | 10.02 | 10.14 | 9.93 | 37.6 | 120 |
| bd2000 | dc:scs:100:astral3 | 5 | 7.42 | 7.54 | 7.31 | 42.8 | 120 |
| bd2000 | dc:scs:500:astral3 | 5 | 8.01 | 8.05 | 7.98 | 61.0 | 471 |
| bd2000 | dc:true:100:astral3 | 5 | 6.17 | 6.33 | 6.04 | 16.3 | 121 |
| bd5000 | scs | 5 | 3.39 | 3.68 | 3.11 | 29.2 | 728 |
| bd5000 | dc:scs:100:astral3 | 5 | 7.85 | 7.96 | 7.76 | 58.2 | 121 |
| bd5000 | dc:scs:500:astral3 | 5 | 8.12 | 8.14 | 8.11 | 69.3 | 476 |
| bd10000 | astral4 | 1 | 29.19 | 29.19 | 29.19 | 1595.7 | 92 |
| bd10000 | mrlft | 1 | 8.03 | 8.03 | 8.03 | 1033.1 | 1315 |
| bd10000 | scs | 3 | 3.84 | 4.19 | 3.51 | 64.6 | 1256 |
| bd10000 | dc:mrlft:100:astral3 | 1 | 11.06 | 11.10 | 11.03 | 122.8 | 122 |
| bd10000 | dc:mrlft:500:astral3 | 1 | 9.52 | 9.53 | 9.51 | 155.2 | 517 |
| bd10000 | dc:scs:100:astral3 | 3 | 8.15 | 8.25 | 8.07 | 112.0 | 122 |
| bd10000 | dc:scs:500:astral3 | 3 | 8.32 | 8.33 | 8.31 | 133.5 | 502 |

### Paired comparisons (RF %, new - baseline; negative = new is better). W/T/L = new wins/ties/losses, tie band |diff| <= 0.25 pp; two-sided Wilcoxon signed-rank

| condition | new | baseline | n | mean diff (pp) | W/T/L | p |
|---|---|---|---|---|---|---|
| bd500 | dc:astral3:100:astral3 | astral3 | 5 | +2.07 | 0/0/5 | 0.062 |
| bd500 | mrlft | astral3 | 5 | +1.13 | 3/0/2 | 0.81 |
| bd500 | dc:scs:100:astral3 | scs | 5 | +3.74 | 0/0/5 | 0.062 |
| bd500 | dc:scs:500:astral3 | scs | 5 | +3.88 | 0/0/5 | 0.062 |
| bd500 | dc:scs:100:astral3 | astral3 | 5 | -0.14 | 2/1/2 | 1 |
| bd500 | dc:scs:500:astral3 | astral3 | 5 | +0.00 | 0/5/0 | 1 |
| bd500 | scs | astral3 | 5 | -3.88 | 5/0/0 | 0.062 |
| bd500 | dc:astral3:100:astral3 | scs | 5 | +5.96 | 0/0/5 | 0.062 |
| bd1000 | dc:astral3:100:astral3 | astral3 | 5 | +2.95 | 0/0/5 | 0.062 |
| bd1000 | mrlft | astral3 | 5 | +0.22 | 2/1/2 | 1 |
| bd1000 | dc:scs:100:astral3 | scs | 5 | +5.32 | 0/0/5 | 0.062 |
| bd1000 | dc:scs:500:astral3 | scs | 5 | +5.26 | 0/0/5 | 0.062 |
| bd1000 | dc:scs:100:astral3 | astral3 | 5 | +0.03 | 2/0/3 | 1 |
| bd1000 | dc:scs:500:astral3 | astral3 | 5 | -0.03 | 3/1/1 | 0.62 |
| bd1000 | scs | astral3 | 5 | -5.29 | 5/0/0 | 0.062 |
| bd1000 | dc:astral3:100:astral3 | scs | 5 | +8.23 | 0/0/5 | 0.062 |
| bd2000 | dc:astral3:100:astral3 | astral3 | 5 | +2.79 | 0/0/5 | 0.062 |
| bd2000 | mrlft | astral3 | 5 | +0.08 | 2/0/3 | 0.81 |
| bd2000 | dc:scs:100:astral3 | scs | 5 | +4.22 | 0/0/5 | 0.062 |
| bd2000 | dc:scs:500:astral3 | scs | 5 | +4.81 | 0/0/5 | 0.062 |
| bd2000 | dc:scs:100:astral3 | astral3 | 5 | +0.19 | 1/1/3 | 0.81 |
| bd2000 | dc:scs:500:astral3 | astral3 | 5 | +0.78 | 0/3/2 | 0.12 |
| bd2000 | scs | astral3 | 5 | -4.03 | 5/0/0 | 0.062 |
| bd2000 | dc:astral3:100:astral3 | scs | 5 | +6.82 | 0/0/5 | 0.062 |
| bd5000 | dc:scs:100:astral3 | scs | 5 | +4.47 | 0/0/5 | 0.062 |
| bd5000 | dc:scs:500:astral3 | scs | 5 | +4.73 | 0/0/5 | 0.062 |
| bd10000 | dc:scs:100:astral3 | scs | 3 | +4.31 | 0/0/3 | 0.25 |
| bd10000 | dc:scs:500:astral3 | scs | 3 | +4.48 | 0/0/3 | 0.25 |
| ALL | dc:astral3:100:astral3 | astral3 | 15 | +2.60 | 0/0/15 | 0.00065 |
| ALL | mrlft | astral3 | 15 | +0.48 | 7/1/7 | 0.52 |
| ALL | dc:scs:100:astral3 | scs | 23 | +4.42 | 0/0/23 | 2.4e-07 |
| ALL | dc:scs:500:astral3 | scs | 23 | +4.65 | 0/0/23 | 2.4e-07 |
| ALL | dc:scs:100:astral3 | astral3 | 15 | +0.02 | 5/2/8 | 0.85 |
| ALL | dc:scs:500:astral3 | astral3 | 15 | +0.25 | 3/9/3 | 0.58 |
| ALL | scs | astral3 | 15 | -4.40 | 15/0/0 | 6.1e-05 |
| ALL | dc:astral3:100:astral3 | scs | 15 | +7.00 | 0/0/15 | 6.1e-05 |
