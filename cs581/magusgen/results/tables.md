## Per dataset: MAGUS error and Δ error (points) of each variant

| dataset | type | MAGUS err | `linsi&fftns2-op3` | `soft0.5:linsi&fftns2-op3` | `linsi\|cons0.7` | `ss:linsi` | `ss:linsi&fftns2-op3` | `ss:soft0.5:linsi&fftns2-op3` | `ss:linsi\|cons0.7` | `soft0.2:linsi&fftns2-op3` | `sm0.5:linsi\|cons0.7` | `ss:soft0.2:linsi&fftns2-op3` | `wsoft0.03:linsi&fftns2` | `ss:wsoft0.03:linsi&fftns2` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1000L1_R0 | DNA/RNA | 7.47 | +7.22 | -0.01 | +9.16 | -0.70 | +4.56 | -0.75 | +10.25 | -0.01 | -0.13 | -0.65 | -0.04 | -0.62 |
| 1000M2_R0 | DNA/RNA | 8.23 | +15.65 | +0.00 | +9.49 | -0.73 |  | -0.68 |  | +0.00 | +0.10 | -0.49 | +0.01 | -0.67 |
| 1000S1_R0 | DNA/RNA | 9.78 | +3.89 | -0.01 | +3.77 | -2.30 |  | -2.29 |  |  |  |  | +0.01 | -2.28 |
| 1000S2_R0 | DNA/RNA | 4.74 | +0.04 | -0.00 | +0.14 | -0.36 |  |  |  |  |  |  | -0.21 | -0.65 |
| 16S.M_R0 | DNA/RNA | 12.86 | -0.40 | -0.07 | -0.37 | -0.06 |  |  |  |  |  |  | -0.64 | -0.34 |
| BBA0039_R0 | protein | 4.69 | -0.16 | -0.04 | -0.04 | -0.01 | -0.12 | -0.07 | -0.02 | -0.11 | -0.04 | -0.05 | -0.06 | +0.01 |
| BBA0067_R0 | protein | 26.28 | -0.96 | +0.07 | -2.47 | -0.31 | -0.97 | -0.58 | -2.53 | +0.13 | -0.21 | -0.61 | -0.68 | -1.04 |
| BBA0101_R0 | protein | 29.19 | -1.59 | -0.11 | -2.94 | -1.22 | -1.85 | -0.76 | -3.25 | -0.76 | -0.35 | -1.02 | -1.89 | -1.83 |
| BBA0117 | protein | 12.61 | +0.37 | -0.00 | -0.04 | -0.21 |  |  |  |  |  |  | +1.07 | +1.02 |
| BBA0134_R0 | protein | 18.34 | -1.26 | +5.08 | -1.50 | +11.07 | timeout | timeout |  |  |  |  | -0.20 | -0.20 |
| BBA0154_R0 | protein | 21.66 | -1.76 | -0.23 | -1.05 | -0.45 |  |  |  |  |  |  | -1.64 | -1.34 |
| BBA0190_R0 | protein | 23.40 | -1.96 | -1.05 | -0.23 | +0.13 |  |  |  |  |  |  | -1.84 | -1.61 |
| RNASim_R0 | DNA/RNA | 9.92 | +0.58 | -0.34 | -0.08 | -0.31 | -0.17 | -0.65 | -0.57 | -0.23 | -0.04 | -0.71 | +0.42 | -0.71 |
| SIMHIGH_R1 | protein | 25.41 | -2.51 | -0.76 | +0.76 | timeout |  |  |  |  |  |  | +1.08 | timeout |
| SIMMOD_R1 | protein | 13.62 | -1.97 | -0.19 | -1.44 | +1.23 |  |  |  |  |  |  | +1.33 | +2.83 |

## Summary per variant (Δ = variant − MAGUS, points; W/T/L with |Δ| < 0.05 a tie; two-sided Wilcoxon)

| variant | group | n | Δ error | W/T/L | p | Δ SPFN | Δ SPFP | extra wall s over MAGUS (median) | extra % of MAGUS wall (median) |
|---|---|---|---|---|---|---|---|---|---|
| `linsi&fftns2-op3` | DNA/RNA | 6 | +4.50 | 1/1/4 | 0.0938 | +9.56 | -0.57 | 74 | 5% |
| `linsi&fftns2-op3` | protein | 9 | -1.31 | 8/0/1 | 0.0117 | -0.15 | -2.48 | -3 | -0% |
| `linsi&fftns2-op3` | pooled | 15 | +1.01 | 9/1/5 | 0.639 | +3.74 | -1.71 | 16 | 2% |
| `soft0.5:linsi&fftns2-op3` | DNA/RNA | 6 | -0.07 | 2/4/0 | 0.0625 | -0.01 | -0.14 | 65 | 5% |
| `soft0.5:linsi&fftns2-op3` | protein | 9 | +0.31 | 5/2/2 | 0.25 | +1.00 | -0.38 | 10 | 1% |
| `soft0.5:linsi&fftns2-op3` | pooled | 15 | +0.16 | 7/6/2 | 0.0479 | +0.60 | -0.29 | 19 | 3% |
| `linsi\|cons0.7` | DNA/RNA | 6 | +3.69 | 2/0/4 | 0.219 | +8.58 | -1.20 | 14 | 1% |
| `linsi\|cons0.7` | protein | 9 | -0.99 | 6/2/1 | 0.0273 | +0.42 | -2.41 | -10 | -1% |
| `linsi\|cons0.7` | pooled | 15 | +0.88 | 8/2/5 | 0.72 | +3.68 | -1.93 | 2 | 0% |
| `ss:linsi` | DNA/RNA | 6 | -0.74 | 6/0/0 | 0.0312 | -0.86 | -0.63 | 140 | 8% |
| `ss:linsi` | protein | 9 | +1.28 | 4/1/4 | 0.945 | +3.12 | -0.56 | 238 | 25% |
| `ss:linsi` | pooled | 15 | +0.41 | 10/1/4 | 0.153 | +1.41 | -0.59 | 169 | 23% |
| `ss:linsi&fftns2-op3` | DNA/RNA | 2 | +2.20 | 1/0/1 | 1 | +3.84 | +0.56 | 406 | 34% |
| `ss:linsi&fftns2-op3` | protein | 4 | -0.98 | 3/0/1 | 0.25 | +0.45 | -2.41 | 208 | 22% |
| `ss:linsi&fftns2-op3` | pooled | 6 | +0.29 | 4/0/2 | 0.625 | +1.81 | -1.23 | 314 | 22% |
| `ss:soft0.5:linsi&fftns2-op3` | DNA/RNA | 4 | -1.09 | 4/0/0 | 0.125 | -1.22 | -0.97 | 240 | 14% |
| `ss:soft0.5:linsi&fftns2-op3` | protein | 4 | -0.47 | 3/0/1 | 0.25 | +0.08 | -1.03 | 208 | 22% |
| `ss:soft0.5:linsi&fftns2-op3` | pooled | 8 | -0.83 | 7/0/1 | 0.0156 | -0.66 | -0.99 | 238 | 20% |
| `ss:linsi\|cons0.7` | DNA/RNA | 2 | +4.84 | 1/0/1 | 1 | +10.93 | -1.26 | 479 | 36% |
| `ss:linsi\|cons0.7` | protein | 3 | -1.93 | 2/1/0 | 0.25 | +1.59 | -5.46 | 168 | 21% |
| `ss:linsi\|cons0.7` | pooled | 5 | +0.78 | 3/1/1 | 0.625 | +5.33 | -3.78 | 172 | 22% |
| `soft0.2:linsi&fftns2-op3` | DNA/RNA | 3 | -0.08 | 1/2/0 | 0.5 | +0.24 | -0.39 | 209 | 12% |
| `soft0.2:linsi&fftns2-op3` | protein | 3 | -0.25 | 2/0/1 | 0.75 | +0.23 | -0.72 | 32 | 4% |
| `soft0.2:linsi&fftns2-op3` | pooled | 6 | -0.16 | 3/2/1 | 0.312 | +0.23 | -0.56 | 54 | 6% |
| `sm0.5:linsi\|cons0.7` | DNA/RNA | 3 | -0.02 | 1/1/1 | 0.75 | +0.06 | -0.10 | 23 | 1% |
| `sm0.5:linsi\|cons0.7` | protein | 3 | -0.20 | 2/1/0 | 0.25 | +0.07 | -0.47 | 12 | 1% |
| `sm0.5:linsi\|cons0.7` | pooled | 6 | -0.11 | 3/2/1 | 0.156 | +0.07 | -0.29 | 18 | 1% |
| `ss:soft0.2:linsi&fftns2-op3` | DNA/RNA | 3 | -0.62 | 3/0/0 | 0.25 | -0.72 | -0.52 | 372 | 25% |
| `ss:soft0.2:linsi&fftns2-op3` | protein | 3 | -0.56 | 2/1/0 | 0.25 | +0.37 | -1.48 | 221 | 27% |
| `ss:soft0.2:linsi&fftns2-op3` | pooled | 6 | -0.59 | 5/1/0 | 0.0312 | -0.17 | -1.00 | 343 | 26% |
| `wsoft0.03:linsi&fftns2` | DNA/RNA | 6 | -0.08 | 2/3/1 | 0.688 | +0.57 | -0.72 | 121 | 8% |
| `wsoft0.03:linsi&fftns2` | protein | 9 | -0.32 | 6/0/3 | 0.426 | +1.98 | -2.61 | 22 | 3% |
| `wsoft0.03:linsi&fftns2` | pooled | 15 | -0.22 | 8/3/4 | 0.359 | +1.42 | -1.85 | 42 | 4% |
| `ss:wsoft0.03:linsi&fftns2` | DNA/RNA | 6 | -0.88 | 6/0/0 | 0.0312 | -0.84 | -0.92 | 265 | 18% |
| `ss:wsoft0.03:linsi&fftns2` | protein | 9 | -0.27 | 5/1/3 | 0.461 | +1.94 | -2.49 | 250 | 28% |
| `ss:wsoft0.03:linsi&fftns2` | pooled | 15 | -0.53 | 11/1/3 | 0.0676 | +0.75 | -1.81 | 250 | 20% |

## Runtime on fresh draws (MAGUS measured on this machine, 4 threads; extra = measured wall of the added steps)

| dataset | MAGUS wall s | `linsi&fftns2-op3` extra s (%) | `soft0.5:linsi&fftns2-op3` extra s (%) | `linsi\|cons0.7` extra s (%) | `ss:linsi` extra s (%) | `ss:linsi&fftns2-op3` extra s (%) | `ss:soft0.5:linsi&fftns2-op3` extra s (%) | `ss:linsi\|cons0.7` extra s (%) | `soft0.2:linsi&fftns2-op3` extra s (%) | `sm0.5:linsi\|cons0.7` extra s (%) | `ss:soft0.2:linsi&fftns2-op3` extra s (%) | `wsoft0.03:linsi&fftns2` extra s (%) | `ss:wsoft0.03:linsi&fftns2` extra s (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BBA0117 | 14 | 2 (+13%) | 1 (+8%) | 1 (+7%) | 4 (+26%) |  |  |  |  |  |  | 2 (+17%) | 7 (+51%) |
| SIMHIGH_R1 | 980 | -15 (-1%) | 24 (+3%) | -8 (-1%) | 942 (+96%) |  |  |  |  |  |  | 42 (+4%) | 985 (+101%) |
| SIMMOD_R1 | 982 | -21 (-2%) | 5 (+1%) | -16 (-2%) | 243 (+25%) |  |  |  |  |  |  | 13 (+1%) | 289 (+29%) |
