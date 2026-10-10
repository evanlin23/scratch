
### 16S_b1000

WITCH default wall 322 s = decomposition 16 s + all-vs-all hmmsearch 243 s + alignment/merge/other 63 s (search = 75% of wall)

| method | SPFN % | SPFP % | dSPFN vs WITCH (pp) | W/T/L vs WITCH | Wilcoxon p | HMMs scored/query | mean k | time (s) | % of WITCH time | % excl. decomposition |
|---|---|---|---|---|---|---|---|---|---|---|
| WITCH (default, own run) | 0.54 | 0.39 | +0.00 | 0/1000/0 | 1 | all | 10 | 322 | 100% | 100% |
| all HMMs, k=10 (WITCH re-run via weights) | 0.54 | 0.39 | +0.00 | 0/1000/0 | 1 | 277.0 | 9.9 | 318 | 99% | 99% |
| all HMMs, k=3 | 0.58 | 0.37 | +0.04 | 10/930/60 | 1.8e-07 | 277.0 | 3.0 | 294 | 91% | 91% |
| all HMMs, k=1 (UPP-like, adj. bit-score) | 0.69 | 0.33 | +0.15 | 22/804/174 | 3.1e-19 | 277.0 | 1.0 | 281 | 87% | 87% |
| all HMMs, adaptive k (tau=.99) | 0.57 | 0.38 | +0.04 | 0/962/38 | 7.7e-08 | 277.0 | 4.9 | 309 | 96% | 96% |
| all HMMs, adaptive k (tau=.95) | 0.60 | 0.38 | +0.06 | 4/920/76 | 7.2e-12 | 277.0 | 3.5 | 299 | 93% | 92% |
| hier descent (UPP2-style), k=10 | 0.73 | 0.48 | +0.19 | 13/871/116 | 1.8e-13 | 15.4 | 9.9 | 88 | 27% | 23% |
| hier descent, adaptive k | 0.80 | 0.48 | +0.26 | 11/841/148 | 4.6e-20 | 15.4 | 4.3 | 77 | 24% | 20% |
| hier EarlyStop (UPP2), k<=10 | 1.90 | 0.90 | +1.36 | 18/508/474 | 5.9e-73 | 7.0 | 5.7 | 62 | 19% | 15% |
| beam-2 descent, k=10 | 0.58 | 0.41 | +0.04 | 6/958/36 | 5.8e-05 | 27.8 | 9.9 | 102 | 32% | 28% |
| beam-2 descent, adaptive k | 0.62 | 0.41 | +0.08 | 5/925/70 | 1.2e-10 | 27.8 | 4.8 | 92 | 29% | 25% |
| BLAST path, k<=10 | 0.56 | 0.36 | +0.02 | 12/935/53 | 0.00031 | 8.2 | 8.0 | 73 | 23% | 19% |
| BLAST path, adaptive k | 0.59 | 0.36 | +0.05 | 12/913/75 | 1.2e-07 | 8.2 | 3.2 | 66 | 21% | 16% |
| BLAST path+siblings, k=10 | 0.54 | 0.37 | -0.00 | 8/970/22 | 0.37 | 15.4 | 9.8 | 86 | 27% | 23% |
| BLAST path+siblings, adaptive k (tau=.99) | 0.58 | 0.37 | +0.04 | 7/937/56 | 5.3e-06 | 15.4 | 3.9 | 79 | 24% | 21% |
| BLAST path+siblings, adaptive k (tau=.95) | 0.62 | 0.37 | +0.08 | 7/899/94 | 3.2e-12 | 15.4 | 2.8 | 71 | 22% | 18% |
| beam-3 descent, k=10 | 0.55 | 0.39 | +0.01 | 4/983/13 | 0.1 | 38.1 | 9.9 | 113 | 35% | 32% |
| beam-3 descent, adaptive k | 0.59 | 0.39 | +0.05 | 3/952/45 | 2.1e-07 | 38.1 | 4.8 | 103 | 32% | 29% |
| beam-4 descent, k=10 | 0.54 | 0.39 | +0.00 | 2/996/2 | 0.72 | 48.3 | 9.9 | 122 | 38% | 35% |
| hybrid: BLAST paths if hit>=100 bits else beam-3 (sel. time est.) | 0.53 | 0.37 | -0.00 | 6/982/12 | 0.81 | 17.0 | 10.0 | 95 | 30% | 26% |
| BLASTN only (-task blastn) | 0.57 | 0.25 | +0.03 | 131/773/96 | 0.23 | - | - | 4 | 1% | 1% |
| BLASTN only (megablast, TIPP3-fast) | 1.22 | 0.34 | +0.68 | 128/748/124 | 0.017 | - | - | 2 | 0% | 1% |

### 16S_b5000

WITCH default wall 1419 s = decomposition 94 s + all-vs-all hmmsearch 1217 s + alignment/merge/other 107 s (search = 86% of wall)

| method | SPFN % | SPFP % | dSPFN vs WITCH (pp) | W/T/L vs WITCH | Wilcoxon p | HMMs scored/query | mean k | time (s) | % of WITCH time | % excl. decomposition |
|---|---|---|---|---|---|---|---|---|---|---|
| WITCH (default, own run) | 1.29 | 0.68 | +0.00 | 0/1000/0 | 1 | all | 10 | 1419 | 100% | 100% |
| all HMMs, k=10 (WITCH re-run via weights) | 1.29 | 0.68 | +0.00 | 0/1000/0 | 1 | 1387.0 | 10.0 | 1394 | 98% | 98% |
| all HMMs, k=3 | 1.27 | 0.65 | -0.02 | 15/971/14 | 0.29 | 1387.0 | 3.0 | 1378 | 97% | 97% |
| all HMMs, k=1 (UPP-like, adj. bit-score) | 1.45 | 0.59 | +0.16 | 30/857/113 | 1.2e-06 | 1387.0 | 1.0 | 1373 | 97% | 97% |
| all HMMs, adaptive k (tau=.99) | 1.34 | 0.66 | +0.05 | 4/962/34 | 3.5e-05 | 1387.0 | 2.9 | 1388 | 98% | 98% |
| all HMMs, adaptive k (tau=.95) | 1.39 | 0.65 | +0.10 | 8/928/64 | 1.5e-07 | 1387.0 | 2.3 | 1383 | 97% | 97% |
| hier descent (UPP2-style), k=10 | 1.59 | 0.79 | +0.30 | 5/902/93 | 7.6e-16 | 19.9 | 9.6 | 201 | 14% | 8% |
| hier descent, adaptive k | 1.65 | 0.77 | +0.36 | 7/872/121 | 1.5e-19 | 19.9 | 2.2 | 203 | 14% | 8% |
| hier EarlyStop (UPP2), k<=10 | 2.34 | 1.04 | +1.05 | 8/728/264 | 1.7e-44 | 3.1 | 2.9 | 177 | 12% | 6% |
| beam-2 descent, k=10 | 1.36 | 0.70 | +0.07 | 3/971/26 | 0.00012 | 37.0 | 9.9 | 231 | 16% | 10% |
| beam-2 descent, adaptive k | 1.41 | 0.68 | +0.12 | 5/938/57 | 1.6e-08 | 37.0 | 2.8 | 226 | 16% | 10% |
| BLAST path, k<=10 | 1.34 | 0.69 | +0.05 | 5/966/29 | 5.7e-05 | 10.4 | 9.6 | 207 | 15% | 9% |
| BLAST path, adaptive k | 1.39 | 0.67 | +0.10 | 8/932/60 | 1.2e-08 | 10.4 | 1.9 | 203 | 14% | 8% |
| BLAST path+siblings, k=10 | 1.31 | 0.68 | +0.02 | 3/983/14 | 0.022 | 19.9 | 9.9 | 220 | 15% | 9% |
| BLAST path+siblings, adaptive k (tau=.99) | 1.37 | 0.66 | +0.08 | 7/943/50 | 4e-06 | 19.9 | 2.3 | 221 | 16% | 10% |
| BLAST path+siblings, adaptive k (tau=.95) | 1.40 | 0.66 | +0.11 | 11/919/70 | 4e-07 | 19.9 | 1.8 | 218 | 15% | 9% |
| beam-3 descent, k=10 | 1.32 | 0.69 | +0.03 | 2/985/13 | 0.009 | 52.1 | 9.9 | 244 | 17% | 11% |
| beam-3 descent, adaptive k | 1.36 | 0.67 | +0.07 | 5/951/44 | 5.1e-06 | 52.1 | 2.9 | 237 | 17% | 11% |
| beam-4 descent, k=10 | 1.32 | 0.69 | +0.03 | 2/989/9 | 0.041 | 66.9 | 9.9 | 253 | 18% | 12% |
| hybrid: BLAST paths if hit>=100 bits else beam-3 (sel. time est.) | 1.30 | 0.68 | +0.01 | 1/993/6 | 0.043 | 24.4 | 10.0 | 231 | 16% | 10% |
| BLASTN only (-task blastn) | 0.48 | 0.21 | -0.81 | 448/495/57 | 5.9e-48 | - | - | 20 | 1% | 2% |
| BLASTN only (megablast, TIPP3-fast) | 0.94 | 0.22 | -0.35 | 443/491/66 | 5.7e-40 | - | - | 6 | 0% | 0% |

### rna_b1000

WITCH default wall 368 s = decomposition 19 s + all-vs-all hmmsearch 280 s + alignment/merge/other 69 s (search = 76% of wall)

| method | SPFN % | SPFP % | dSPFN vs WITCH (pp) | W/T/L vs WITCH | Wilcoxon p | HMMs scored/query | mean k | time (s) | % of WITCH time | % excl. decomposition |
|---|---|---|---|---|---|---|---|---|---|---|
| WITCH (default, own run) | 2.72 | 1.42 | +0.00 | 0/1000/0 | 1 | all | 10 | 368 | 100% | 100% |
| all HMMs, k=10 (WITCH re-run via weights) | 2.72 | 1.42 | +0.00 | 0/1000/0 | 1 | 283.0 | 10.0 | 361 | 98% | 98% |
| all HMMs, k=3 | 2.86 | 1.41 | +0.14 | 6/875/119 | 1.6e-19 | 283.0 | 3.0 | 350 | 95% | 95% |
| all HMMs, k=1 (UPP-like, adj. bit-score) | 3.31 | 1.32 | +0.60 | 41/575/384 | 3.2e-48 | 283.0 | 1.0 | 334 | 91% | 90% |
| all HMMs, adaptive k (tau=.99) | 2.84 | 1.40 | +0.13 | 5/877/118 | 4.5e-18 | 283.0 | 2.8 | 351 | 95% | 95% |
| all HMMs, adaptive k (tau=.95) | 3.01 | 1.39 | +0.30 | 3/773/224 | 1.7e-37 | 283.0 | 2.0 | 345 | 94% | 94% |
| hier descent (UPP2-style), k=10 | 3.22 | 1.62 | +0.51 | 9/900/91 | 1.3e-16 | 15.5 | 10.0 | 108 | 29% | 26% |
| hier descent, adaptive k | 3.35 | 1.61 | +0.63 | 8/802/190 | 1.6e-32 | 15.5 | 2.7 | 96 | 26% | 22% |
| hier EarlyStop (UPP2), k<=10 | 3.91 | 1.69 | +1.19 | 29/614/357 | 4e-51 | 6.6 | 6.0 | 89 | 24% | 20% |
| beam-2 descent, k=10 | 2.76 | 1.44 | +0.04 | 3/981/16 | 0.0015 | 27.9 | 10.0 | 121 | 33% | 29% |
| beam-2 descent, adaptive k | 2.89 | 1.42 | +0.17 | 6/861/133 | 1.9e-21 | 27.9 | 2.8 | 112 | 31% | 27% |
| BLAST path, k<=10 | 2.99 | 1.41 | +0.27 | 16/914/70 | 1.7e-06 | 8.2 | 8.1 | 97 | 26% | 22% |
| BLAST path, adaptive k | 3.11 | 1.39 | +0.39 | 17/823/160 | 3e-20 | 8.2 | 2.5 | 87 | 24% | 20% |
| BLAST path+siblings, k=10 | 2.80 | 1.42 | +0.08 | 8/951/41 | 0.00015 | 15.3 | 10.0 | 108 | 29% | 25% |
| BLAST path+siblings, adaptive k (tau=.99) | 2.92 | 1.40 | +0.21 | 10/847/143 | 2.3e-19 | 15.3 | 2.7 | 98 | 27% | 23% |
| BLAST path+siblings, adaptive k (tau=.95) | 3.07 | 1.39 | +0.35 | 8/752/240 | 7.3e-37 | 15.3 | 2.0 | 93 | 25% | 21% |
| beam-3 descent, k=10 | 2.73 | 1.44 | +0.02 | 0/993/7 | 0.018 | 38.3 | 10.0 | 138 | 37% | 34% |
| beam-3 descent, adaptive k | 2.86 | 1.41 | +0.15 | 5/869/126 | 2e-19 | 38.3 | 2.8 | 125 | 34% | 31% |
| beam-4 descent, k=10 | 2.72 | 1.43 | +0.00 | 0/997/3 | 0.11 | 48.8 | 10.0 | 146 | 40% | 36% |
| hybrid: BLAST paths if hit>=100 bits else beam-3 (sel. time est.) | 2.73 | 1.43 | +0.01 | 1/990/9 | 0.022 | 23.8 | 10.0 | 118 | 32% | 28% |
| BLASTN only (-task blastn) | 11.12 | 1.95 | +8.40 | 175/78/747 | 4.2e-104 | - | - | 2 | 1% | 1% |
| BLASTN only (megablast, TIPP3-fast) | 78.71 | 2.58 | +75.99 | 25/12/963 | 3.3e-162 | - | - | 1 | 0% | 0% |

### roseM1_R0

WITCH default wall 61 s = decomposition 7 s + all-vs-all hmmsearch 40 s + alignment/merge/other 15 s (search = 65% of wall)

| method | SPFN % | SPFP % | dSPFN vs WITCH (pp) | W/T/L vs WITCH | Wilcoxon p | HMMs scored/query | mean k | time (s) | % of WITCH time | % excl. decomposition |
|---|---|---|---|---|---|---|---|---|---|---|
| WITCH (default, own run) | 7.85 | 2.30 | +0.00 | 0/300/0 | 1 | all | 10 | 61 | 100% | 100% |
| all HMMs, k=10 (WITCH re-run via weights) | 7.85 | 2.30 | +0.00 | 0/300/0 | 1 | 203.0 | 8.9 | 60 | 99% | 99% |
| all HMMs, k=3 | 8.44 | 1.70 | +0.59 | 4/244/52 | 7.4e-07 | 203.0 | 3.0 | 57 | 93% | 92% |
| all HMMs, k=1 (UPP-like, adj. bit-score) | 10.33 | 1.39 | +2.48 | 4/164/132 | 2.5e-21 | 203.0 | 1.0 | 55 | 91% | 90% |
| all HMMs, adaptive k (tau=.99) | 8.12 | 2.29 | +0.28 | 0/265/35 | 2.5e-07 | 203.0 | 3.2 | 59 | 97% | 96% |
| all HMMs, adaptive k (tau=.95) | 8.47 | 2.20 | +0.63 | 2/223/75 | 5.9e-12 | 203.0 | 2.5 | 58 | 95% | 95% |
| hier descent (UPP2-style), k=10 | 10.96 | 3.24 | +3.11 | 4/255/41 | 1e-07 | 14.4 | 9.0 | 26 | 42% | 35% |
| hier descent, adaptive k | 11.25 | 3.20 | +3.41 | 3/228/69 | 2.3e-12 | 14.4 | 3.0 | 24 | 40% | 32% |
| hier EarlyStop (UPP2), k<=10 | 19.13 | 4.63 | +11.29 | 12/121/167 | 7.6e-29 | 7.8 | 5.8 | 23 | 37% | 30% |
| beam-2 descent, k=10 | 8.95 | 2.59 | +1.11 | 1/285/14 | 0.0015 | 26.1 | 9.0 | 28 | 46% | 40% |
| beam-2 descent, adaptive k | 9.21 | 2.58 | +1.37 | 1/253/46 | 1.3e-08 | 26.1 | 3.2 | 27 | 44% | 38% |
| BLAST path, k<=10 | 17.72 | 4.25 | +9.87 | 16/143/141 | 3.7e-24 | 20.0 | 6.2 | 26 | 43% | 36% |
| BLAST path, adaptive k | 17.84 | 4.06 | +10.00 | 16/134/150 | 1.5e-25 | 20.0 | 2.2 | 25 | 40% | 33% |
| BLAST path+siblings, k=10 | 14.20 | 4.11 | +6.35 | 16/151/133 | 6e-22 | 26.2 | 8.7 | 29 | 47% | 41% |
| BLAST path+siblings, adaptive k (tau=.99) | 14.34 | 4.03 | +6.50 | 16/141/143 | 1.7e-23 | 26.2 | 2.8 | 27 | 45% | 38% |
| BLAST path+siblings, adaptive k (tau=.95) | 14.67 | 3.96 | +6.82 | 13/126/161 | 5.6e-26 | 26.2 | 2.1 | 25 | 42% | 34% |
| beam-3 descent, k=10 | 8.18 | 2.18 | +0.34 | 0/290/10 | 0.0051 | 35.7 | 9.0 | 30 | 49% | 43% |
| beam-3 descent, adaptive k | 8.45 | 2.16 | +0.60 | 0/257/43 | 1.1e-08 | 35.7 | 3.2 | 28 | 47% | 40% |
| beam-4 descent, k=10 | 8.04 | 2.22 | +0.19 | 0/294/6 | 0.028 | 45.5 | 9.0 | 32 | 53% | 47% |
| hybrid: BLAST paths if hit>=100 bits else beam-3 (sel. time est.) | 8.18 | 2.18 | +0.34 | 0/290/10 | 0.0051 | 33.0 | 10.0 | 30 | 50% | 44% |
| BLASTN only (-task blastn) | 67.65 | 10.01 | +59.80 | 53/11/236 | 1.3e-40 | - | - | 0 | 1% | 1% |
| BLASTN only (megablast, TIPP3-fast) | 89.16 | 0.10 | +81.31 | 27/7/266 | 4.8e-48 | - | - | 0 | 1% | 1% |

### roseM1_R1

WITCH default wall 58 s = decomposition 6 s + all-vs-all hmmsearch 38 s + alignment/merge/other 14 s (search = 66% of wall)

| method | SPFN % | SPFP % | dSPFN vs WITCH (pp) | W/T/L vs WITCH | Wilcoxon p | HMMs scored/query | mean k | time (s) | % of WITCH time | % excl. decomposition |
|---|---|---|---|---|---|---|---|---|---|---|
| WITCH (default, own run) | 10.03 | 4.42 | +0.00 | 0/300/0 | 1 | all | 10 | 58 | 100% | 100% |
| all HMMs, k=10 (WITCH re-run via weights) | 10.03 | 4.42 | +0.00 | 0/300/0 | 1 | 199.0 | 8.7 | 58 | 99% | 98% |
| all HMMs, k=3 | 10.66 | 3.92 | +0.63 | 0/259/41 | 2.4e-08 | 199.0 | 3.0 | 54 | 92% | 92% |
| all HMMs, k=1 (UPP-like, adj. bit-score) | 11.47 | 2.67 | +1.44 | 6/174/120 | 5e-16 | 199.0 | 1.0 | 53 | 90% | 89% |
| all HMMs, adaptive k (tau=.99) | 10.38 | 4.39 | +0.36 | 0/268/32 | 8e-07 | 199.0 | 3.1 | 56 | 96% | 96% |
| all HMMs, adaptive k (tau=.95) | 10.68 | 4.30 | +0.65 | 0/244/56 | 7.5e-11 | 199.0 | 2.5 | 55 | 95% | 94% |
| hier descent (UPP2-style), k=10 | 16.38 | 6.46 | +6.35 | 3/247/50 | 9.8e-10 | 14.5 | 8.6 | 25 | 42% | 35% |
| hier descent, adaptive k | 16.54 | 6.36 | +6.52 | 3/225/72 | 3.4e-13 | 14.5 | 3.1 | 23 | 40% | 33% |
| hier EarlyStop (UPP2), k<=10 | 24.93 | 7.22 | +14.90 | 12/111/177 | 5.4e-30 | 7.4 | 5.4 | 22 | 37% | 30% |
| beam-2 descent, k=10 | 11.83 | 4.91 | +1.80 | 1/281/18 | 0.0004 | 25.9 | 8.7 | 26 | 45% | 39% |
| beam-2 descent, adaptive k | 12.18 | 4.85 | +2.15 | 1/253/46 | 9.1e-09 | 25.9 | 3.1 | 25 | 44% | 37% |
| BLAST path, k<=10 | 18.80 | 6.07 | +8.77 | 7/174/119 | 9e-20 | 17.2 | 6.2 | 24 | 41% | 33% |
| BLAST path, adaptive k | 18.84 | 5.90 | +8.81 | 7/167/126 | 7.9e-21 | 17.2 | 2.1 | 22 | 38% | 31% |
| BLAST path+siblings, k=10 | 16.70 | 6.35 | +6.67 | 9/181/110 | 5.3e-18 | 23.5 | 8.3 | 27 | 47% | 40% |
| BLAST path+siblings, adaptive k (tau=.99) | 16.88 | 6.20 | +6.85 | 7/174/119 | 1.1e-19 | 23.5 | 2.8 | 26 | 44% | 37% |
| BLAST path+siblings, adaptive k (tau=.95) | 17.08 | 5.80 | +7.05 | 5/167/128 | 2.9e-21 | 23.5 | 2.1 | 24 | 41% | 34% |
| beam-3 descent, k=10 | 11.51 | 4.86 | +1.49 | 0/291/9 | 0.0077 | 35.6 | 8.7 | 29 | 50% | 44% |
| beam-3 descent, adaptive k | 11.87 | 4.85 | +1.84 | 0/259/41 | 2.4e-08 | 35.6 | 3.2 | 28 | 48% | 41% |
| beam-4 descent, k=10 | 11.18 | 4.85 | +1.15 | 0/293/7 | 0.018 | 45.3 | 8.7 | 31 | 53% | 47% |
| hybrid: BLAST paths if hit>=100 bits else beam-3 (sel. time est.) | 11.51 | 4.86 | +1.49 | 0/291/9 | 0.0077 | 32.1 | 10.0 | 28 | 48% | 42% |
| BLASTN only (-task blastn) | 61.42 | 7.04 | +51.39 | 64/15/221 | 8.9e-37 | - | - | 0 | 1% | 1% |
| BLASTN only (megablast, TIPP3-fast) | 88.83 | 0.10 | +78.80 | 26/10/264 | 2.1e-47 | - | - | 0 | 1% | 1% |

### Placement delta error (EPA-ng, best-LWR edge; mean added missing true bipartitions per query)

| instance | WITCH | all HMMs k=1 | hier (UPP2-style) k=10 | beam-3 k=10 | hybrid | BLAST path+sib k=10 | BLASTN -task blastn |
|---|---|---|---|---|---|---|---|
| rna_b1000 | 1.117 | 1.130 (+0.013, p=0.46) | 1.157 (+0.040, p=0.6) | 1.118 (+0.001, p=0.32) | 1.118 (+0.001, p=0.32) | 1.114 (-0.003, p=0.33) | 1.903 (+0.789, p=3e-15, unplaced 1) |
| roseM1_R0 | 0.830 | 0.790 (-0.040, p=0.34) | 1.123 (+0.293, p=0.087) | 0.823 (-0.007, p=0.53) | 0.823 (-0.007, p=0.53) | 1.193 (+0.363, p=0.066) | 6.796 (+6.303, p=1.5e-17, unplaced 99) |
| roseM1_R1 | 0.797 | 0.760 (-0.037, p=0.99) | 1.393 (+0.597, p=0.00025) | 0.947 (+0.150, p=0.063) | 0.947 (+0.150, p=0.063) | 1.437 (+0.640, p=1.3e-05) | 4.714 (+4.102, p=8.3e-15, unplaced 94) |
