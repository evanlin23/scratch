## MAGUS with a different base method (merge pilot)

Same MAGUS decomposition (25 subsets) and the same 10 MAFFT L-INS-i backbones as MAGUS's own run; only the subset aligner changes, then MAGUS's GCM merge (paper flags). Control = MAGUS's own L-INS-i subsets (reproduces MAGUS(4c)). Errors in %, Δ in points (negative = better than MAGUS).

| rep | base method | subset err (mean of 25) | MAGUS final err | Δ final vs control | subset align s (1 core) |
|---|---|---|---|---|---|
| 1000L3_R0 | L-INS-i (MAGUS's own) | 6.22 | 11.28 | – | – |
| 1000L3_R0 | mafft-ginsi | 5.96 | 11.23 | -0.05 | 310 |
| 1000M2_R0 | L-INS-i (MAGUS's own) | 3.88 | 8.23 | – | – |
| 1000M2_R0 | mafft-ginsi | 3.59 | 7.97 | -0.25 | 288 |
| 1000S1_R0 | L-INS-i (MAGUS's own) | 5.08 | 9.78 | – | – |
| 1000S1_R0 | mafft-ginsi | 5.03 | 9.93 | +0.16 | 274 |
| 16S.M_R0 | L-INS-i (MAGUS's own) | 6.97 | 12.86 | – | – |
| 16S.M_R0 | mafft-ginsi | 6.97 | 12.89 | +0.03 | 113 |
| BBA0039_R0 | L-INS-i (MAGUS's own) | 2.72 | 4.69 | – | – |
| BBA0039_R0 | clustalo | 2.50 | 4.45 | -0.24 | 9 |
| BBA0039_R0 | clustalo (subsets+backbones) | 2.50 | 4.24 | -0.45 | 9 |
| BBA0039_R0 | mafft-ginsi | 2.75 | 4.64 | -0.05 | 28 |
| BBA0039_R0 | muscle5 | 2.58 | 4.48 | -0.21 | 127 |
| BBA0067_R0 | L-INS-i (MAGUS's own) | 19.44 | 26.28 | – | – |
| BBA0067_R0 | clustalo | 18.93 | 25.97 | -0.31 | 4 |
| BBA0067_R0 | clustalo (backbones only) | 19.44 | 24.72 | -1.56 | 0 |
| BBA0067_R0 | clustalo (subsets+backbones) | 18.93 | 24.86 | -1.43 | 5 |
| BBA0067_R0 | mafft-ginsi | 19.40 | 26.20 | -0.09 | 16 |
| BBA0067_R0 | muscle5 | 18.58 | 25.37 | -0.91 | 43 |
| BBA0101_R0 | L-INS-i (MAGUS's own) | 16.37 | 29.19 | – | – |
| BBA0101_R0 | clustalo | 16.00 | 29.00 | -0.19 | 7 |
| BBA0101_R0 | clustalo (backbones only) | 16.37 | 26.35 | -2.83 | 0 |
| BBA0101_R0 | clustalo (subsets+backbones) | 16.00 | 26.39 | -2.79 | 6 |
| BBA0101_R0 | mafft-ginsi | 16.28 | 29.32 | +0.14 | 19 |
| BBA0101_R0 | muscle5 | 16.11 | 28.11 | -1.08 | 61 |
| BBA0154_R0 | L-INS-i (MAGUS's own) | 7.95 | 21.66 | – | – |
| BBA0154_R0 | clustalo | 7.47 | 21.16 | -0.50 | 6 |
| BBA0154_R0 | clustalo (backbones only) | 7.95 | 19.40 | -2.27 | 0 |
| BBA0154_R0 | clustalo (subsets+backbones) | 7.47 | 18.70 | -2.97 | 6 |
| BBA0154_R0 | mafft-ginsi | 8.06 | 21.65 | -0.01 | 16 |
| BBA0154_R0 | muscle5 | 7.78 | 21.86 | +0.20 | 54 |
| BBA0190_R0 | L-INS-i (MAGUS's own) | 7.26 | 23.40 | – | – |
| BBA0190_R0 | clustalo | 6.96 | 22.98 | -0.42 | 20 |
| BBA0190_R0 | clustalo (subsets+backbones) | 6.96 | 21.59 | -1.81 | 19 |
| BBA0190_R0 | mafft-ginsi | 7.40 | 23.55 | +0.15 | 34 |
| BBA0190_R0 | muscle5 | 7.04 | 22.98 | -0.42 | 130 |
| RNASim_R0 | L-INS-i (MAGUS's own) | 4.92 | 9.92 | – | – |
| RNASim_R0 | mafft-ginsi | 4.67 | 9.72 | -0.20 | 621 |

| data | base method | n | mean Δ | better/tie/worse |
|---|---|---|---|---|
| 16S | mafft-ginsi | 1 | +0.03 | 0/1/0 |
| BAliBASE | clustalo | 5 | -0.33 | 5/0/0 |
| BAliBASE | clustalo (backbones only) | 3 | -2.22 | 3/0/0 |
| BAliBASE | clustalo (subsets+backbones) | 5 | -1.89 | 5/0/0 |
| BAliBASE | mafft-ginsi | 5 | +0.03 | 1/2/2 |
| BAliBASE | muscle5 | 5 | -0.48 | 4/0/1 |
| RNASim | mafft-ginsi | 1 | -0.20 | 1/0/0 |
| ROSE | mafft-ginsi | 3 | -0.05 | 1/1/1 |
