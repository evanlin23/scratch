# Fast MAGUS pilot: paired comparison with MAGUS (paper flags), same machine

Error = (SPFN+SPFP)/2 in %. Delta < 0: lower error than MAGUS. W/T/L: replicates where the variant is better / within 0.05 pts / worse. Runtime ratio = MAGUS wall / variant wall (> 1: faster), mean over replicates.

## Per replicate (error %, wall s)

| variant | 1000L1_R0 | 1000M2_R0 | 1000S2_R0 | BBA0101_R0 |
|---|---|---|---|---|
| magus | 7.88 (1552) | 9.74 (1323) | 3.89 (1049) | 28.03 (989) |
| magus+ss | 7.09 (1675) | 8.76 (1470) | 3.08 (1116) | 27.60 (1185) |
| pt-r4 | - | 44.75 (652) | - | - |
| pt-r4+ss | - | 39.69 (1272) | - | - |
| cl-r4 | - | 26.59 (577) | - | - |
| r4 | - | 11.87 (740) | - | - |
| r4+ss | - | 10.04 (892) | - | - |
| r4+ss2 | - | 9.69 (1063) | - | - |
| sk100-r4 | 8.93 (732) | 11.37 (581) | 4.71 (453) | 28.83 (385) |
| sk100-r4+ss | 8.01 (867) | 9.02 (734) | 3.90 (523) | 27.95 (593) |
| sk100-r4+ss2 | 7.83 (996) | 8.79 (903) | 3.81 (589) | 28.10 (736) |
| sk100-r2 | - | 15.97 (398) | - | - |
| sk100-r2+ss | - | 13.18 (573) | - | - |
| sk100-r2+ss2 | - | 12.49 (762) | - | - |

## Paired summary

| variant | n | mean err | mean delta (pts) | W/T/L | Wilcoxon p | mean wall s | speedup vs MAGUS (mean of ratios) | CPU ratio |
|---|---|---|---|---|---|---|---|---|
| magus | 4 | 12.38 | +0.00 | 0/4/0 | - | 1228 | 1.00x | 1.00x |
| magus+ss | 4 | 11.63 | -0.75 | 4/0/0 | 0.125 | 1361 | 0.90x | 0.95x |
| pt-r4 | 1 | 44.75 | +35.01 | 0/0/1 | - | 652 | 2.03x | 2.05x |
| pt-r4+ss | 1 | 39.69 | +29.95 | 0/0/1 | - | 1272 | 1.04x | 1.38x |
| cl-r4 | 1 | 26.59 | +16.85 | 0/0/1 | - | 577 | 2.29x | 2.28x |
| r4 | 1 | 11.87 | +2.13 | 0/0/1 | - | 740 | 1.79x | 1.82x |
| r4+ss | 1 | 10.04 | +0.30 | 0/0/1 | - | 892 | 1.48x | 1.63x |
| r4+ss2 | 1 | 9.69 | -0.05 | 0/1/0 | - | 1063 | 1.24x | 1.44x |
| sk100-r4 | 4 | 13.46 | +1.07 | 0/0/4 | 0.125 | 538 | 2.32x | 2.39x |
| sk100-r4+ss | 4 | 12.22 | -0.16 | 2/1/1 | 0.875 | 679 | 1.82x | 2.12x |
| sk100-r4+ss2 | 4 | 12.13 | -0.25 | 2/1/1 | 0.375 | 806 | 1.54x | 1.92x |
| sk100-r2 | 1 | 15.97 | +6.23 | 0/0/1 | - | 398 | 3.33x | 3.52x |
| sk100-r2+ss | 1 | 13.18 | +3.44 | 0/0/1 | - | 573 | 2.31x | 2.80x |
| sk100-r2+ss2 | 1 | 12.49 | +2.75 | 0/0/1 | - | 762 | 1.74x | 2.26x |

## Against the published MAGUS(Fast) alignment of the same replicate

| variant | n | mean delta vs published MAGUS(Fast) (pts) | W/T/L |
|---|---|---|---|
| magus | 4 | +0.31 | 2/0/2 |
| magus+ss | 4 | -0.44 | 3/0/1 |
| pt-r4 | 1 | +36.50 | 0/0/1 |
| pt-r4+ss | 1 | +31.43 | 0/0/1 |
| cl-r4 | 1 | +18.33 | 0/0/1 |
| r4 | 1 | +3.62 | 0/0/1 |
| r4+ss | 1 | +1.78 | 0/0/1 |
| r4+ss2 | 1 | +1.44 | 0/0/1 |
| sk100-r4 | 4 | +1.38 | 0/0/4 |
| sk100-r4+ss | 4 | +0.15 | 2/0/2 |
| sk100-r4+ss2 | 4 | +0.05 | 2/0/2 |
| sk100-r2 | 1 | +7.71 | 0/0/1 |
| sk100-r2+ss | 1 | +4.93 | 0/0/1 |
| sk100-r2+ss2 | 1 | +4.24 | 0/0/1 |

## MAGUS (paper flags) stage breakdown on this machine

| replicate | wall s | guide tree s (%) | subsets+backbones until graph built s (%) | cluster+trace s (%) |
|---|---|---|---|---|
| 1000L1_R0 | 1552 | 288 (19%) | 1245 (80%) | 18 (1%) |
| 1000M2_R0 | 1323 | 253 (19%) | 1052 (79%) | 19 (1%) |
| 1000S2_R0 | 1049 | 205 (20%) | 839 (80%) | 5 (0%) |
| BBA0101_R0 | 989 | 175 (18%) | 776 (78%) | 38 (4%) |

PASTA(3) (published alignments and logs, same 4 replicates): 2.17x MAGUS(Fast)'s runtime, +2.14 pts error vs MAGUS(Fast).
