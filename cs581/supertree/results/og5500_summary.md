### Mean RF error (%) / mean FN (%) / mean wall (s) / mean peak RSS (MB), n replicates

| condition | method | n | RF % | FN % | FP % | wall s | peak MB |
|---|---|---|---|---|---|---|---|
| og5500 | scs | 5 | 32.41 | 39.31 | 29.59 | 310.8 | 2723 |
| og5500 | dc:scs:200:astral3 | 5 | 45.04 | 45.13 | 45.03 | 107.1 | 788 |
| og5500 | dc:scs:500:astral3 | 1 | 46.35 | 46.38 | 46.35 | 172.9 | 1609 |

### Paired comparisons (RF %, new - baseline; negative = new is better). W/T/L = new wins/ties/losses, tie band |diff| <= 0.25 pp; two-sided Wilcoxon signed-rank

| condition | new | baseline | n | mean dRF (pp) | W/T/L (RF) | p (RF) | mean dFN (pp) | W/T/L (FN) | p (FN) |
|---|---|---|---|---|---|---|---|---|---|
| og5500 | dc:scs:200:astral3 | scs | 5 | +12.63 | 0/0/5 | 0.062 | +5.82 | 0/0/5 | 0.062 |
| ALL | dc:scs:200:astral3 | scs | 5 | +12.63 | 0/0/5 | 0.062 | +5.82 | 0/0/5 | 0.062 |
