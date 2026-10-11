## Final MAGUS alignment SP error (%) by subset aligner (merge-only, same subsets & backbones)

| dataset | linsi | ginsi | muscle5 | famsa | probcons | clustalo | kalign | prank |
|---|---|---|---|---|---|---|---|---|
| BBA0039 | 4.66 | 4.64 | 4.48 | 4.92 | 4.61 | **4.45** | 5.03 |  |
| BBA0067 | 26.40 | 26.20 | 25.37 | 26.53 | **25.33** | 25.97 | 27.27 |  |
| BBA0101 | 28.88 | 29.32 | **28.11** | 29.63 | 28.65 | 29.00 | 29.21 |  |
| BBA0154 | 21.76 | 21.65 | 21.86 | 21.56 | 21.72 | **21.16** | 21.31 |  |
| BBA0190 | 23.43 | 23.55 | 22.98 | 23.77 | 22.98 | **22.98** | 23.63 |  |
| HF_Acetyltransf | **26.03** | 34.07 | 27.90 | 40.67 | 44.79 | 36.75 | 41.41 |  |
| HF_PDZ | 16.61 | **15.39** | 17.79 | 17.53 | 18.05 | 16.84 | 26.17 |  |
| HF_aat | 15.07 | **14.23** | 15.46 | 15.67 | 15.56 | 16.41 | 17.49 |  |
| HF_rrm | 18.41 | 19.61 | 20.81 | 20.09 | 20.51 | 20.48 | **16.98** |  |
| HF_sdr | 23.35 | **21.93** | 26.48 | 26.73 | 27.55 | 27.50 | 35.08 |  |
| HF_zf-CCHH | 12.25 | 12.33 | 14.83 | 12.13 | 15.51 | 15.51 | **9.32** |  |
| SIMMOD_R1 | 12.79 | 12.83 | 13.82 | 16.41 | 13.73 | 20.38 | 14.99 |  |
| SIMMOD_R2 | 11.06 | **10.86** | 11.41 | 13.88 | 11.46 | 17.06 | 12.07 |  |
| SIMHIGH_R1 | 25.00 | **24.19** | 28.45 | 30.21 | 28.06 | 34.53 | 30.68 |  |
| SIMHIGH_R2 | 23.63 | **22.92** | 23.57 | 29.32 | 24.29 | 29.36 | 26.51 |  |
| 1000M2 | 8.24 | **7.97** |  | 23.01 |  | 13.26 | 11.52 |  |

## Paired differences vs L-INS-i subsets (SP-error points, negative = better than MAGUS default)

| method | set | n | mean d | median d | W/T/L | p (Wilcoxon) | d SPFN | d SPFP | d TC |
|---|---|---|---|---|---|---|---|---|---|
| ginsi | all protein | 15 | +0.29 | -0.11 | 8/2/5 | 0.421 | +0.22 | +0.37 | +0.1 |
| ginsi | BAliBASE | 5 | +0.05 | -0.02 | 2/1/2 | 1 | +0.22 | -0.13 | +0.3 |
| ginsi | HomFam | 6 | +0.97 | -0.39 | 3/0/3 | 1 | +0.75 | +1.19 | -0.7 |
| ginsi | SIMMOD | 2 | -0.08 | -0.08 | 1/1/0 | 1 | +0.09 | -0.24 | +1.0 |
| ginsi | SIMHIGH | 2 | -0.76 | -0.76 | 2/0/0 | 0.5 | -1.27 | -0.25 | +1.3 |
| ginsi | nucleotide | 1 | -0.26 | -0.26 | 1/0/0 | nan | -0.25 | -0.27 | +0.7 |
| muscle5 | all protein | 15 | +0.93 | +0.39 | 5/0/10 | 0.0554 | +1.22 | +0.64 | -0.8 |
| muscle5 | BAliBASE | 5 | -0.47 | -0.45 | 4/0/1 | 0.125 | -0.25 | -0.68 | +1.4 |
| muscle5 | HomFam | 6 | +1.92 | +2.13 | 0/0/6 | 0.0312 | +2.58 | +1.26 | -0.4 |
| muscle5 | SIMMOD | 2 | +0.69 | +0.69 | 0/0/2 | 0.5 | +0.54 | +0.84 | -5.3 |
| muscle5 | SIMHIGH | 2 | +1.69 | +1.69 | 1/0/1 | 1 | +1.48 | +1.90 | -3.0 |
| famsa | all protein | 15 | +2.65 | +0.93 | 2/0/13 | 0.000427 | +2.64 | +2.66 | -4.1 |
| famsa | BAliBASE | 5 | +0.26 | +0.26 | 1/0/4 | 0.188 | +0.44 | +0.07 | -0.5 |
| famsa | HomFam | 6 | +3.52 | +1.30 | 1/0/5 | 0.0625 | +3.22 | +3.82 | -3.0 |
| famsa | SIMMOD | 2 | +3.22 | +3.22 | 0/0/2 | 0.5 | +3.66 | +2.77 | -10.6 |
| famsa | SIMHIGH | 2 | +5.45 | +5.45 | 0/0/2 | 0.5 | +5.39 | +5.51 | -9.8 |
| famsa | nucleotide | 1 | +14.78 | +14.78 | 0/0/1 | nan | +15.33 | +14.22 | -12.4 |
| probcons | all protein | 15 | +2.23 | +0.67 | 4/1/10 | 0.0215 | +2.53 | +1.94 | -1.7 |
| probcons | BAliBASE | 5 | -0.37 | -0.23 | 4/1/0 | 0.0625 | -0.01 | -0.72 | +0.8 |
| probcons | HomFam | 6 | +5.04 | +2.68 | 0/0/6 | 0.0312 | +5.55 | +4.53 | -0.6 |
| probcons | SIMMOD | 2 | +0.67 | +0.67 | 0/0/2 | 0.5 | +0.46 | +0.88 | -7.6 |
| probcons | SIMHIGH | 2 | +1.86 | +1.86 | 0/0/2 | 0.5 | +1.86 | +1.87 | -4.9 |
| clustalo | all protein | 15 | +3.27 | +2.07 | 4/0/11 | 0.0125 | +3.13 | +3.41 | -4.6 |
| clustalo | BAliBASE | 5 | -0.31 | -0.43 | 4/0/1 | 0.125 | -0.15 | -0.48 | +1.3 |
| clustalo | HomFam | 6 | +3.63 | +2.66 | 0/0/6 | 0.0312 | +3.36 | +3.89 | -1.9 |
| clustalo | SIMMOD | 2 | +6.80 | +6.80 | 0/0/2 | 0.5 | +7.20 | +6.40 | -16.0 |
| clustalo | SIMHIGH | 2 | +7.63 | +7.63 | 0/0/2 | 0.5 | +6.55 | +8.71 | -15.9 |
| clustalo | nucleotide | 1 | +5.03 | +5.03 | 0/0/1 | nan | +5.38 | +4.67 | -13.5 |
| kalign | all protein | 15 | +3.19 | +1.01 | 3/0/12 | 0.0302 | +6.61 | -0.23 | -5.5 |
| kalign | BAliBASE | 5 | +0.26 | +0.33 | 1/0/4 | 0.438 | +3.39 | -2.86 | +0.6 |
| kalign | HomFam | 6 | +5.79 | +5.99 | 2/0/4 | 0.219 | +11.09 | +0.48 | -11.6 |
| kalign | SIMMOD | 2 | +1.61 | +1.61 | 0/0/2 | 0.5 | +2.23 | +0.98 | -1.9 |
| kalign | SIMHIGH | 2 | +4.28 | +4.28 | 0/0/2 | 0.5 | +5.59 | +2.97 | -6.1 |
| kalign | nucleotide | 1 | +3.28 | +3.28 | 0/0/1 | nan | +7.32 | -0.76 | +3.3 |

## Pre-registered primary test (protein sets, Holm over 4 candidates)

| candidate | n | mean d vs L-INS-i | W/T/L | raw p | Holm p |
|---|---|---|---|---|---|
| muscle5 | 15 | +0.93 | 5/0/10 | 0.0554 | 0.0554 |
| probcons | 15 | +2.23 | 4/1/10 | 0.0215 | 0.0431 |
| famsa | 15 | +2.65 | 2/0/13 | 0.000427 | 0.00171 |
| clustalo | 15 | +3.27 | 4/0/11 | 0.0125 | 0.0374 |

## Subset-level accuracy (pooled over subsets; HomFam only subsets with >= 2 seed sequences)

| set | linsi | ginsi | muscle5 | famsa | probcons | clustalo | kalign | prank |
|---|---|---|---|---|---|---|---|---|
| BAliBASE (n=5) | 10.43 | 10.48 | 10.14 | 10.77 | 10.21 | 10.06 | 10.91 |  |
| HomFam (n=6) | 8.86 | 9.40 | 12.30 | 11.19 | 17.60 | 11.84 | 14.17 |  |
| SIMMOD (n=2) | 3.35 | 3.12 | 4.26 | 6.01 | 4.30 | 10.18 | 5.02 |  |
| SIMHIGH (n=2) | 7.69 | 7.43 | 9.65 | 12.88 | 9.87 | 17.77 | 12.23 |  |
| nucleotide (n=1) | 3.91 | 3.66 |  | 17.40 |  | 8.37 | 6.94 |  |

## Subset-alignment cost: total CPU seconds over the 25 subsets (1 thread each); ratio to L-INS-i

| set | linsi | ginsi | muscle5 | famsa | probcons | clustalo | kalign | prank |
|---|---|---|---|---|---|---|---|---|
| BAliBASE | 81 (1x) | 77 (0.95x) | 317 (3.9x) | 24 (0.3x) | 627 (7.7x) | 45 (0.55x) | 17 (0.21x) |  |
| HomFam | 215 (1x) | 202 (0.94x) | 1478 (6.9x) | 6 (0.03x) | 1995 (9.3x) | 38 (0.17x) | 6 (0.028x) |  |
| SIMMOD | 164 (1x) | 178 (1.1x) | 1100 (6.7x) | 27 (0.17x) | 1560 (9.5x) | 53 (0.32x) | 12 (0.071x) |  |
| SIMHIGH | 234 (1x) | 251 (1.1x) | 1267 (5.4x) | 26 (0.11x) | 1632 (7x) | 64 (0.27x) | 12 (0.053x) |  |
| nucleotide | 1184 (1x) | 1306 (1.1x) |  | 71 (0.06x) |  | 313 (0.26x) | 568 (0.48x) |  |

GCM merge wall time (4 threads), L-INS-i subsets: median 27 s, range 1-54 s.

## Control: GCM merge on re-aligned L-INS-i subsets vs the original MAGUS run (SP error %)

| dataset | MAGUS (original) | merge-linsi |
|---|---|---|
| BBA0039 | 4.69 | 4.66 |
| BBA0067 | 26.28 | 26.40 |
| BBA0101 | 29.18 | 28.88 |
| BBA0154 | 21.66 | 21.76 |
| BBA0190 | 23.40 | 23.43 |
| HF_Acetyltransf | 26.76 | 26.03 |
| HF_PDZ | 16.14 | 16.61 |
| HF_aat | 15.16 | 15.07 |
| HF_rrm | 18.57 | 18.41 |
| HF_sdr | 23.65 | 23.35 |
| HF_zf-CCHH | 12.25 | 12.25 |
| SIMMOD_R1 | 12.71 | 12.79 |
| SIMMOD_R2 | 10.99 | 11.06 |
| SIMHIGH_R1 | 25.04 | 25.00 |
| SIMHIGH_R2 | 23.51 | 23.63 |
| 1000M2 | 8.23 | 8.24 |

## Whole-dataset baselines and PASTA variants (SP error %, wall s, 4 threads)

| dataset | MAGUS (L-INS-i subsets) | best MAGUS variant | famsa | muscle5 | pasta-mafft-it1 | pasta-prank-it1 | pasta-probcons-it1 | regressive |
|---|---|---|---|---|---|---|---|---|
| BBA0039 | 4.66 | 4.45 (clustalo) | 5.49 (268 s) | fail/timeout | 4.48 (542 s) | fail/timeout | 4.56 (1752 s) | 18.30 (33 s) |
| BBA0067 | 26.40 | 25.33 (probcons) | 32.04 (343 s) | fail/timeout |  |  |  | 36.14 (22 s) |
| BBA0101 | 28.88 | 28.11 (muscle5) | 35.24 (284 s) |  |  |  |  | 38.30 (25 s) |
| BBA0154 | 21.76 | 21.16 (clustalo) | 25.26 (255 s) |  |  |  |  | 28.19 (25 s) |
| BBA0190 | 23.43 | 22.98 (clustalo) | 33.67 (436 s) |  |  |  |  | 31.52 (77 s) |
| HF_Acetyltransf | 26.03 | 26.03 (linsi) | 51.32 (20 s) |  |  |  |  |  |
| HF_PDZ | 16.61 | 15.39 (ginsi) | 12.32 (21 s) |  |  |  |  |  |
| HF_aat | 15.07 | 14.23 (ginsi) | 15.89 (124 s) |  |  |  |  |  |
| HF_rrm | 18.41 | 16.98 (kalign) | 18.60 (24 s) | 24.71 (155 s) |  |  |  | 24.71 (41 s) |
| HF_sdr | 23.35 | 21.93 (ginsi) | 24.51 (128 s) |  |  |  |  |  |
| HF_zf-CCHH | 12.25 | 9.32 (kalign) | 10.06 (24 s) |  |  |  |  | 15.27 (32 s) |
| SIMMOD_R1 | 12.79 | 12.71 (orig) | 19.99 (282 s) | fail/timeout |  |  |  | 54.00 (60 s) |
| SIMMOD_R2 | 11.06 | 10.86 (ginsi) | 16.09 (295 s) |  |  |  |  |  |
| SIMHIGH_R1 | 25.00 | 24.19 (ginsi) | 38.22 (374 s) |  |  |  |  | 71.50 (68 s) |
| SIMHIGH_R2 | 23.63 | 22.92 (ginsi) | 29.55 (372 s) |  |  |  |  |  |
| 1000L1 |  |  | 33.08 (1558 s) |  |  |  |  |  |
| 1000M2 | 8.24 | 7.97 (ginsi) | 43.59 (922 s) |  | 12.96 (2029 s) | fail/timeout | fail/timeout | 81.16 (275 s) |
| 16S.M |  |  | 16.50 (555 s) |  |  |  |  |  |
| RNASim |  |  | fail/timeout |  |  |  |  |  |

## Trees: FastTree -lg -gamma, normalized RF (%) vs the true tree

| dataset | merge_ginsi | merge_linsi | merge_muscle5 | merge_probcons | true |
|---|---|---|---|---|---|
| SIMHIGH_R1 | 6.82 | 10.13 | 10.63 | 8.83 | 6.12 |
| SIMHIGH_R2 | 8.53 | 6.72 | 8.32 | 7.22 | 5.22 |
| SIMMOD_R1 |  | 6.32 | 6.12 | 6.22 | 6.52 |
| SIMMOD_R2 | 7.02 | 7.52 | 7.82 | 7.82 | 6.52 |
