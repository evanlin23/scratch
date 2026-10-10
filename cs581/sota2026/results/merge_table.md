## MAGUS with a different base method (merge pilot)

Same MAGUS decomposition (25 subsets) and the same 10 MAFFT L-INS-i backbones as MAGUS's own run; only the subset aligner changes, then MAGUS's GCM merge (paper flags). Control = MAGUS's own L-INS-i subsets (reproduces MAGUS(4c)). Errors in %, Δ in points (negative = better than MAGUS).

| rep | base method | subset err (mean of 25) | MAGUS final err | Δ final vs control | subset align s (1 core) |
|---|---|---|---|---|---|
| 1000L3_R0 | L-INS-i (MAGUS's own) | 6.22 | 11.28 | – | – |
| 1000M2_R0 | L-INS-i (MAGUS's own) | 3.88 | 8.23 | – | – |
| 16S.M_R0 | L-INS-i (MAGUS's own) | 6.97 | 12.86 | – | – |
| BBA0067_R0 | L-INS-i (MAGUS's own) | 19.44 | 26.28 | – | – |
| BBA0067_R0 | clustalo | 18.93 | 25.97 | -0.31 | 4 |
| BBA0067_R0 | mafft-ginsi | 19.40 | 26.20 | -0.09 | 16 |
| BBA0067_R0 | muscle5 | 18.58 | 25.37 | -0.91 | 43 |
| BBA0101_R0 | L-INS-i (MAGUS's own) | 16.37 | 29.19 | – | – |
| BBA0101_R0 | clustalo | 16.00 | 29.00 | -0.19 | 7 |
| BBA0101_R0 | mafft-ginsi | 16.28 | 29.32 | +0.14 | 19 |
| BBA0101_R0 | muscle5 | 16.11 | 28.11 | -1.08 | 61 |
| BBA0154_R0 | L-INS-i (MAGUS's own) | 7.95 | 21.66 | – | – |
| BBA0154_R0 | clustalo | 7.47 | 21.16 | -0.50 | 6 |
| BBA0154_R0 | mafft-ginsi | 8.06 | 21.65 | -0.01 | 16 |
| BBA0154_R0 | muscle5 | 7.78 | 21.86 | +0.20 | 54 |
| BBA0190_R0 | L-INS-i (MAGUS's own) | 7.26 | 23.40 | – | – |
| BBA0190_R0 | mafft-ginsi | 7.40 | 23.55 | +0.15 | 34 |
| RNASim_R0 | L-INS-i (MAGUS's own) | 4.92 | 9.92 | – | – |

| data | base method | n | mean Δ | better/tie/worse |
|---|---|---|---|---|
| BAliBASE | clustalo | 3 | -0.33 | 3/0/0 |
| BAliBASE | mafft-ginsi | 4 | +0.05 | 1/1/2 |
| BAliBASE | muscle5 | 3 | -0.60 | 2/0/1 |
