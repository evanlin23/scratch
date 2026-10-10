## Per dataset: MAGUS error and Δ error (points) of each variant

| dataset | type | MAGUS err | `linsi&fftns2-op3` | `soft0.5:linsi&fftns2-op3` | `linsi\|cons0.7` | `ss:linsi` | `ss:linsi&fftns2-op3` | `ss:soft0.5:linsi&fftns2-op3` | `ss:linsi\|cons0.7` | `soft0.2:linsi&fftns2-op3` | `sm0.5:linsi\|cons0.7` | `ss:soft0.2:linsi&fftns2-op3` | `wsoft0.03:linsi&fftns2` | `ss:wsoft0.03:linsi&fftns2` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1000L1_R0 | DNA/RNA | 7.47 | +7.22 | -0.01 | +9.16 | -0.70 | +4.56 | -0.75 | +10.25 | -0.01 | -0.13 | -0.65 | -0.04 | -0.62 |
| 1000M2_R0 | DNA/RNA | 8.23 | +15.65 | +0.00 | +9.49 | -0.73 |  | -0.68 |  | +0.00 | +0.10 | -0.49 | +0.01 | -0.67 |
| 1000S1_R0 | DNA/RNA | 9.78 | +3.89 | -0.01 | +3.77 |  |  |  |  |  |  |  |  |  |
| BBA0039_R0 | protein | 4.69 | -0.16 | -0.04 | -0.04 | -0.01 | -0.12 | -0.07 | -0.02 | -0.11 | -0.04 | -0.05 | -0.06 | +0.01 |
| BBA0067_R0 | protein | 26.28 | -0.96 | +0.07 | -2.47 | -0.31 | -0.97 | -0.58 | -2.53 | +0.13 | -0.21 | -0.61 | -0.68 | -1.04 |
| BBA0101_R0 | protein | 29.19 | -1.59 | -0.11 | -2.94 | -1.22 | -1.85 | -0.76 | -3.25 | -0.76 | -0.35 | -1.02 | -1.89 | -1.83 |
| BBA0134_R0 | protein | 18.34 | -1.26 | +5.08 | -1.50 | +11.07 | timeout | timeout |  |  |  |  | -0.20 | -0.20 |
| RNASim_R0 | DNA/RNA | 9.92 | +0.58 | -0.34 | -0.08 | -0.31 | -0.17 | -0.65 | -0.57 | -0.23 | -0.04 | -0.71 | +0.42 | -0.71 |

## Summary per variant (Δ = variant − MAGUS, points; W/T/L with |Δ| < 0.05 a tie; two-sided Wilcoxon)

| variant | group | n | Δ error | W/T/L | p | Δ SPFN | Δ SPFP | extra wall s over MAGUS (median) |
|---|---|---|---|---|---|---|---|---|
| `linsi&fftns2-op3` | DNA/RNA | 4 | +6.84 | 0/0/4 | 0.125 | +13.70 | -0.02 | 75 |
| `linsi&fftns2-op3` | protein | 4 | -0.99 | 4/0/0 | 0.125 | +0.42 | -2.40 | 1 |
| `linsi&fftns2-op3` | pooled | 8 | +2.92 | 4/0/4 | 0.547 | +7.06 | -1.21 | 50 |
| `soft0.5:linsi&fftns2-op3` | DNA/RNA | 4 | -0.09 | 1/3/0 | 0.25 | -0.04 | -0.14 | 71 |
| `soft0.5:linsi&fftns2-op3` | protein | 4 | +1.25 | 1/1/2 | 0.875 | +2.56 | -0.06 | 10 |
| `soft0.5:linsi&fftns2-op3` | pooled | 8 | +0.58 | 2/4/2 | 0.641 | +1.26 | -0.10 | 43 |
| `linsi\|cons0.7` | DNA/RNA | 4 | +5.59 | 1/0/3 | 0.25 | +12.68 | -1.51 | 10 |
| `linsi\|cons0.7` | protein | 4 | -1.74 | 3/1/0 | 0.125 | +0.82 | -4.30 | -23 |
| `linsi\|cons0.7` | pooled | 8 | +1.92 | 4/1/3 | 0.742 | +6.75 | -2.91 | 5 |
| `ss:linsi` | DNA/RNA | 3 | -0.58 | 3/0/0 | 0.25 | -0.73 | -0.43 | 169 |
| `ss:linsi` | protein | 4 | +2.38 | 2/1/1 | 0.875 | +5.79 | -1.02 | 212 |
| `ss:linsi` | pooled | 7 | +1.11 | 5/1/1 | 0.297 | +3.00 | -0.77 | 187 |
| `ss:linsi&fftns2-op3` | DNA/RNA | 2 | +2.20 | 1/0/1 | 1 | +3.84 | +0.56 | 406 |
| `ss:linsi&fftns2-op3` | protein | 4 | -0.98 | 3/0/1 | 0.25 | +0.45 | -2.41 | 208 |
| `ss:linsi&fftns2-op3` | pooled | 6 | +0.29 | 4/0/2 | 0.625 | +1.81 | -1.23 | 314 |
| `ss:soft0.5:linsi&fftns2-op3` | DNA/RNA | 3 | -0.69 | 3/0/0 | 0.25 | -0.81 | -0.57 | 242 |
| `ss:soft0.5:linsi&fftns2-op3` | protein | 4 | -0.47 | 3/0/1 | 0.25 | +0.08 | -1.03 | 208 |
| `ss:soft0.5:linsi&fftns2-op3` | pooled | 7 | -0.58 | 6/0/1 | 0.0312 | -0.36 | -0.80 | 239 |
| `ss:linsi\|cons0.7` | DNA/RNA | 2 | +4.84 | 1/0/1 | 1 | +10.93 | -1.26 | 479 |
| `ss:linsi\|cons0.7` | protein | 3 | -1.93 | 2/1/0 | 0.25 | +1.59 | -5.46 | 168 |
| `ss:linsi\|cons0.7` | pooled | 5 | +0.78 | 3/1/1 | 0.625 | +5.33 | -3.78 | 172 |
| `soft0.2:linsi&fftns2-op3` | DNA/RNA | 3 | -0.08 | 1/2/0 | 0.5 | +0.24 | -0.39 | 209 |
| `soft0.2:linsi&fftns2-op3` | protein | 3 | -0.25 | 2/0/1 | 0.75 | +0.23 | -0.72 | 32 |
| `soft0.2:linsi&fftns2-op3` | pooled | 6 | -0.16 | 3/2/1 | 0.312 | +0.23 | -0.56 | 54 |
| `sm0.5:linsi\|cons0.7` | DNA/RNA | 3 | -0.02 | 1/1/1 | 0.75 | +0.06 | -0.10 | 23 |
| `sm0.5:linsi\|cons0.7` | protein | 3 | -0.20 | 2/1/0 | 0.25 | +0.07 | -0.47 | 12 |
| `sm0.5:linsi\|cons0.7` | pooled | 6 | -0.11 | 3/2/1 | 0.156 | +0.07 | -0.29 | 18 |
| `ss:soft0.2:linsi&fftns2-op3` | DNA/RNA | 3 | -0.62 | 3/0/0 | 0.25 | -0.72 | -0.52 | 372 |
| `ss:soft0.2:linsi&fftns2-op3` | protein | 3 | -0.56 | 2/1/0 | 0.25 | +0.37 | -1.48 | 221 |
| `ss:soft0.2:linsi&fftns2-op3` | pooled | 6 | -0.59 | 5/1/0 | 0.0312 | -0.17 | -1.00 | 343 |
| `wsoft0.03:linsi&fftns2` | DNA/RNA | 3 | +0.13 | 0/2/1 | 0.75 | +0.99 | -0.73 | 124 |
| `wsoft0.03:linsi&fftns2` | protein | 4 | -0.71 | 4/0/0 | 0.125 | +1.26 | -2.68 | 21 |
| `wsoft0.03:linsi&fftns2` | pooled | 7 | -0.35 | 4/2/1 | 0.219 | +1.14 | -1.84 | 47 |
| `ss:wsoft0.03:linsi&fftns2` | DNA/RNA | 3 | -0.67 | 3/0/0 | 0.25 | -0.68 | -0.66 | 304 |
| `ss:wsoft0.03:linsi&fftns2` | protein | 4 | -0.77 | 3/1/0 | 0.25 | +1.26 | -2.79 | 242 |
| `ss:wsoft0.03:linsi&fftns2` | pooled | 7 | -0.72 | 6/1/0 | 0.0312 | +0.43 | -1.88 | 288 |
