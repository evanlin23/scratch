
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
| BLAST path, k<=10 | 0.56 | 0.36 | +0.02 | 13/934/53 | 0.00045 | 8.7 | 8.0 | 70 | 22% | 18% |
| BLAST path, adaptive k | 0.60 | 0.36 | +0.06 | 12/913/75 | 4.3e-08 | 8.7 | 3.2 | 64 | 20% | 16% |
| BLAST path+siblings, k=10 | 0.55 | 0.37 | +0.01 | 8/966/26 | 0.13 | 15.9 | 9.8 | 85 | 26% | 22% |
| BLAST path+siblings, adaptive k (tau=.99) | 0.58 | 0.37 | +0.04 | 7/936/57 | 3.8e-06 | 15.9 | 3.9 | 76 | 24% | 20% |
| BLAST path+siblings, adaptive k (tau=.95) | 0.63 | 0.37 | +0.09 | 7/897/96 | 1.3e-12 | 15.9 | 2.8 | 69 | 21% | 17% |
| BLASTN only (TIPP3-fast) | 1.22 | 0.34 | +0.68 | 128/748/124 | 0.017 | - | - | 2 | 0% | 1% |
