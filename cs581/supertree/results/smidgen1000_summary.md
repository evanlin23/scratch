### Mean RF error (%) / mean FN (%) / mean wall (s) / mean peak RSS (MB), n replicates

| condition | method | n | RF % | FN % | FP % | wall s | peak MB |
|---|---|---|---|---|---|---|---|
| s1000_d20 | astral3 | 10 | 17.19 | 17.19 | 17.19 | 20.8 | 1332 |
| s1000_d20 | astral4 | 8 | 32.07 | 32.07 | 32.07 | 44.9 | 10 |
| s1000_d20 | mrlft | 10 | 19.64 | 19.64 | 19.64 | 16.9 | 37 |
| s1000_d20 | dc:astral3:100:astral3 | 10 | 17.71 | 17.73 | 17.70 | 11.5 | 123 |
| s1000_d20 | dc:mrlft:100:astral3 | 10 | 20.10 | 20.15 | 20.07 | 11.8 | 124 |
| s1000_d20 | dc:mrlft:200:astral3 | 10 | 19.65 | 19.66 | 19.64 | 12.5 | 213 |
| s1000_d50 | astral3 | 10 | 15.76 | 15.76 | 15.76 | 21.4 | 1518 |
| s1000_d50 | mrlft | 10 | 20.88 | 20.88 | 20.88 | 27.2 | 45 |
| s1000_d50 | dc:astral3:100:astral3 | 10 | 15.55 | 15.55 | 15.54 | 11.9 | 119 |
| s1000_d50 | dc:mrlft:100:astral3 | 10 | 19.92 | 19.97 | 19.90 | 11.9 | 121 |
| s1000_d50 | dc:mrlft:200:astral3 | 10 | 19.27 | 19.31 | 19.25 | 13.2 | 215 |
| s1000_d75 | astral3 | 10 | 13.85 | 13.85 | 13.85 | 22.9 | 1544 |
| s1000_d75 | mrlft | 10 | 17.42 | 17.42 | 17.42 | 25.9 | 48 |
| s1000_d75 | dc:astral3:100:astral3 | 10 | 13.74 | 13.75 | 13.73 | 12.2 | 118 |
| s1000_d75 | dc:mrlft:100:astral3 | 10 | 16.76 | 16.78 | 16.74 | 12.5 | 120 |
| s1000_d75 | dc:mrlft:200:astral3 | 10 | 16.52 | 16.53 | 16.52 | 13.7 | 217 |
| s1000_d100 | astral3 | 10 | 11.62 | 11.62 | 11.62 | 23.9 | 1662 |
| s1000_d100 | mrlft | 10 | 11.93 | 11.93 | 11.93 | 24.0 | 48 |
| s1000_d100 | dc:astral3:100:astral3 | 10 | 11.35 | 11.36 | 11.35 | 13.3 | 119 |
| s1000_d100 | dc:mrlft:100:astral3 | 10 | 11.69 | 11.69 | 11.69 | 13.3 | 117 |
| s1000_d100 | dc:mrlft:200:astral3 | 10 | 11.67 | 11.68 | 11.67 | 13.8 | 212 |

### Paired comparisons (RF %, new - baseline; negative = new is better). W/T/L = new wins/ties/losses, tie band |diff| <= 0.25 pp; two-sided Wilcoxon signed-rank

| condition | new | baseline | n | mean diff (pp) | W/T/L | p |
|---|---|---|---|---|---|---|
| s1000_d20 | dc:astral3:100:astral3 | astral3 | 10 | +0.53 | 4/0/6 | 0.23 |
| s1000_d20 | dc:mrlft:100:astral3 | astral3 | 10 | +2.91 | 0/0/10 | 0.002 |
| s1000_d20 | dc:mrlft:200:astral3 | astral3 | 10 | +2.46 | 0/1/9 | 0.0039 |
| s1000_d20 | dc:mrlft:100:astral3 | mrlft | 10 | +0.46 | 2/2/6 | 0.16 |
| s1000_d20 | dc:mrlft:200:astral3 | mrlft | 10 | +0.01 | 4/1/5 | 0.92 |
| s1000_d20 | mrlft | astral3 | 10 | +2.45 | 0/0/10 | 0.002 |
| s1000_d20 | astral4 | astral3 | 8 | +13.99 | 0/0/8 | 0.0078 |
| s1000_d50 | dc:astral3:100:astral3 | astral3 | 10 | -0.21 | 4/3/3 | 0.57 |
| s1000_d50 | dc:mrlft:100:astral3 | astral3 | 10 | +4.17 | 0/0/10 | 0.002 |
| s1000_d50 | dc:mrlft:200:astral3 | astral3 | 10 | +3.51 | 0/0/10 | 0.002 |
| s1000_d50 | dc:mrlft:100:astral3 | mrlft | 10 | -0.96 | 7/0/3 | 0.064 |
| s1000_d50 | dc:mrlft:200:astral3 | mrlft | 10 | -1.61 | 8/2/0 | 0.002 |
| s1000_d50 | mrlft | astral3 | 10 | +5.12 | 0/0/10 | 0.002 |
| s1000_d75 | dc:astral3:100:astral3 | astral3 | 10 | -0.11 | 6/1/3 | 0.62 |
| s1000_d75 | dc:mrlft:100:astral3 | astral3 | 10 | +2.91 | 0/0/10 | 0.002 |
| s1000_d75 | dc:mrlft:200:astral3 | astral3 | 10 | +2.68 | 0/1/9 | 0.002 |
| s1000_d75 | dc:mrlft:100:astral3 | mrlft | 10 | -0.66 | 4/5/1 | 0.098 |
| s1000_d75 | dc:mrlft:200:astral3 | mrlft | 10 | -0.90 | 5/4/1 | 0.16 |
| s1000_d75 | mrlft | astral3 | 10 | +3.57 | 0/0/10 | 0.002 |
| s1000_d100 | dc:astral3:100:astral3 | astral3 | 10 | -0.27 | 7/1/2 | 0.23 |
| s1000_d100 | dc:mrlft:100:astral3 | astral3 | 10 | +0.06 | 5/3/2 | 0.67 |
| s1000_d100 | dc:mrlft:200:astral3 | astral3 | 10 | +0.05 | 2/5/3 | 0.94 |
| s1000_d100 | dc:mrlft:100:astral3 | mrlft | 10 | -0.24 | 5/3/2 | 0.17 |
| s1000_d100 | dc:mrlft:200:astral3 | mrlft | 10 | -0.26 | 7/0/3 | 0.28 |
| s1000_d100 | mrlft | astral3 | 10 | +0.30 | 4/1/5 | 0.68 |
| ALL | dc:astral3:100:astral3 | astral3 | 40 | -0.02 | 21/5/14 | 0.65 |
| ALL | dc:mrlft:100:astral3 | astral3 | 40 | +2.51 | 5/3/32 | 4.5e-07 |
| ALL | dc:mrlft:200:astral3 | astral3 | 40 | +2.17 | 2/7/31 | 8.4e-07 |
| ALL | dc:mrlft:100:astral3 | mrlft | 40 | -0.35 | 18/10/12 | 0.11 |
| ALL | dc:mrlft:200:astral3 | mrlft | 40 | -0.69 | 24/7/9 | 0.0016 |
| ALL | mrlft | astral3 | 40 | +2.86 | 4/1/35 | 2.8e-07 |
| ALL | astral4 | astral3 | 8 | +13.99 | 0/0/8 | 0.0078 |
