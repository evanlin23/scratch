
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
| BLASTN only (-task blastn) | 0.57 | 0.25 | +0.03 | 131/773/96 | 0.23 | - | - | 4 | 1% | 1% |
| BLASTN only (megablast, TIPP3-fast) | 1.22 | 0.34 | +0.68 | 128/748/124 | 0.017 | - | - | 2 | 0% | 1% |

### rna_b1000

WITCH default wall 407 s = decomposition 19 s + all-vs-all hmmsearch 320 s + alignment/merge/other 68 s (search = 79% of wall)

| method | SPFN % | SPFP % | dSPFN vs WITCH (pp) | W/T/L vs WITCH | Wilcoxon p | HMMs scored/query | mean k | time (s) | % of WITCH time | % excl. decomposition |
|---|---|---|---|---|---|---|---|---|---|---|
| WITCH (default, own run) | 2.72 | 1.42 | +0.00 | 0/1000/0 | 1 | all | 10 | 407 | 100% | 100% |
| all HMMs, k=10 (WITCH re-run via weights) | 2.72 | 1.42 | +0.00 | 0/1000/0 | 1 | 283.0 | 10.0 | 361 | 89% | 88% |
| all HMMs, k=3 | 2.86 | 1.41 | +0.14 | 6/875/119 | 1.6e-19 | 283.0 | 3.0 | 350 | 86% | 85% |
| all HMMs, k=1 (UPP-like, adj. bit-score) | 3.31 | 1.32 | +0.60 | 41/575/384 | 3.2e-48 | 283.0 | 1.0 | 334 | 82% | 81% |
| all HMMs, adaptive k (tau=.99) | 2.84 | 1.40 | +0.13 | 5/877/118 | 4.5e-18 | 283.0 | 2.8 | 351 | 86% | 86% |
| all HMMs, adaptive k (tau=.95) | 3.01 | 1.39 | +0.30 | 3/773/224 | 1.7e-37 | 283.0 | 2.0 | 345 | 85% | 84% |
| hier descent (UPP2-style), k=10 | 3.22 | 1.62 | +0.51 | 9/900/91 | 1.3e-16 | 15.5 | 10.0 | 108 | 27% | 23% |
| hier descent, adaptive k | 3.35 | 1.61 | +0.63 | 8/802/190 | 1.6e-32 | 15.5 | 2.7 | 96 | 24% | 20% |
| hier EarlyStop (UPP2), k<=10 | 3.91 | 1.69 | +1.19 | 29/614/357 | 4e-51 | 6.6 | 6.0 | 89 | 22% | 18% |
| beam-2 descent, k=10 | 2.76 | 1.44 | +0.04 | 3/981/16 | 0.0015 | 27.9 | 10.0 | 121 | 30% | 26% |
| beam-2 descent, adaptive k | 2.89 | 1.42 | +0.17 | 6/861/133 | 1.9e-21 | 27.9 | 2.8 | 112 | 28% | 24% |
| BLAST path, k<=10 | 2.99 | 1.41 | +0.27 | 16/914/70 | 1.7e-06 | 8.2 | 8.1 | 97 | 24% | 20% |
| BLAST path, adaptive k | 3.11 | 1.39 | +0.39 | 17/823/160 | 3e-20 | 8.2 | 2.5 | 87 | 21% | 18% |
| BLAST path+siblings, k=10 | 2.80 | 1.42 | +0.08 | 8/951/41 | 0.00015 | 15.3 | 10.0 | 107 | 26% | 23% |
| BLAST path+siblings, adaptive k (tau=.99) | 2.92 | 1.40 | +0.21 | 10/847/143 | 2.3e-19 | 15.3 | 2.7 | 98 | 24% | 20% |
| BLAST path+siblings, adaptive k (tau=.95) | 3.07 | 1.39 | +0.35 | 8/752/240 | 7.3e-37 | 15.3 | 2.0 | 93 | 23% | 19% |
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
| BLASTN only (-task blastn) | 67.65 | 10.01 | +59.80 | 53/11/236 | 1.3e-40 | - | - | 0 | 1% | 1% |
| BLASTN only (megablast, TIPP3-fast) | 89.16 | 0.10 | +81.31 | 27/7/266 | 4.8e-48 | - | - | 0 | 1% | 1% |

### roseM1_R1

WITCH default wall 65 s = decomposition 6 s + all-vs-all hmmsearch 44 s + alignment/merge/other 14 s (search = 68% of wall)

| method | SPFN % | SPFP % | dSPFN vs WITCH (pp) | W/T/L vs WITCH | Wilcoxon p | HMMs scored/query | mean k | time (s) | % of WITCH time | % excl. decomposition |
|---|---|---|---|---|---|---|---|---|---|---|
| WITCH (default, own run) | 10.03 | 4.42 | +0.00 | 0/300/0 | 1 | all | 10 | 65 | 100% | 100% |
| BLASTN only (-task blastn) | 61.42 | 7.04 | +51.39 | 64/15/221 | 8.9e-37 | - | - | 0 | 1% | 1% |
| BLASTN only (megablast, TIPP3-fast) | 88.83 | 0.10 | +78.80 | 26/10/264 | 2.1e-47 | - | - | 0 | 1% | 1% |
