
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
| BLASTN only (-task blastn) | 11.12 | 1.95 | +8.40 | 175/78/747 | 4.2e-104 | - | - | 2 | 1% | 1% |
| BLASTN only (megablast, TIPP3-fast) | 78.71 | 2.58 | +75.99 | 25/12/963 | 3.3e-162 | - | - | 1 | 0% | 0% |
