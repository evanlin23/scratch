## Final MAGUS alignment SP error (%) by subset aligner (merge-only, same subsets & backbones)

| dataset | linsi | ginsi | muscle5 | famsa | probcons | clustalo | kalign | prank |
|---|---|---|---|---|---|---|---|---|
| BBA0039 | 4.66 | 4.64 | 4.48 | 4.92 | 4.61 | **4.45** | 5.03 |  |
| BBA0067 | 26.40 | 26.20 | 25.37 | 26.53 | **25.33** | 25.97 | 27.27 |  |
| BBA0101 | 28.88 | 29.32 | **28.11** | 29.63 | 28.65 | 29.00 | 29.21 |  |
| BBA0154 | 21.76 | 21.65 | 21.86 | 21.56 | 21.72 | **21.16** | 21.31 |  |
| BBA0190 | 23.43 | 23.55 | 22.98 | 23.77 | 22.98 | **22.98** | 23.63 |  |
| HF_PDZ | 16.61 | **15.39** | 17.79 | 17.53 | 18.05 | 16.84 | 26.17 |  |
| HF_rrm | 18.41 | 19.61 | 20.81 | 20.09 | 20.51 | 20.48 | **16.98** |  |
| HF_zf-CCHH | 12.25 | 12.33 | 14.83 | 12.13 | 15.51 | 15.51 | **9.32** |  |
| SIMMOD_R1 | 12.79 | 12.83 | 13.82 | 16.41 | 13.73 | 20.38 | 14.99 |  |
| SIMHIGH_R1 | 25.00 | **24.19** | 28.45 | 30.21 | 28.06 | 34.53 | 30.68 |  |

## Paired differences vs L-INS-i subsets (SP-error points, negative = better than MAGUS default)

| method | set | n | mean d | median d | W/T/L | p (Wilcoxon) | d SPFN | d SPFP | d TC |
|---|---|---|---|---|---|---|---|---|---|
| ginsi | all protein | 10 | -0.05 | +0.01 | 4/2/4 | 0.922 | -0.00 | -0.09 | +0.0 |
| ginsi | BAliBASE | 5 | +0.05 | -0.02 | 2/1/2 | 1 | +0.22 | -0.13 | +0.3 |
| ginsi | HomFam | 3 | +0.02 | +0.07 | 1/0/2 | 1 | +0.00 | +0.03 | -1.0 |
| ginsi | SIMMOD | 1 | +0.04 | +0.04 | 0/1/0 | nan | +0.25 | -0.16 | +0.3 |
| ginsi | SIMHIGH | 1 | -0.81 | -0.81 | 1/0/0 | nan | -1.42 | -0.21 | +1.5 |
| muscle5 | all protein | 10 | +0.83 | +0.57 | 4/0/6 | 0.232 | +1.31 | +0.35 | -0.1 |
| muscle5 | BAliBASE | 5 | -0.47 | -0.45 | 4/0/1 | 0.125 | -0.25 | -0.68 | +1.4 |
| muscle5 | HomFam | 3 | +2.05 | +2.40 | 0/0/3 | 0.25 | +2.88 | +1.22 | -0.1 |
| muscle5 | SIMMOD | 1 | +1.03 | +1.03 | 0/0/1 | nan | +0.96 | +1.10 | -4.9 |
| muscle5 | SIMHIGH | 1 | +3.44 | +3.44 | 0/0/1 | nan | +4.71 | +2.17 | -2.5 |
| famsa | all protein | 10 | +1.26 | +0.55 | 2/0/8 | 0.0137 | +1.38 | +1.14 | -2.2 |
| famsa | BAliBASE | 5 | +0.26 | +0.26 | 1/0/4 | 0.188 | +0.44 | +0.07 | -0.5 |
| famsa | HomFam | 3 | +0.83 | +0.93 | 1/0/2 | 0.5 | +0.85 | +0.80 | +0.0 |
| famsa | SIMMOD | 1 | +3.62 | +3.62 | 0/0/1 | nan | +4.18 | +3.06 | -11.2 |
| famsa | SIMHIGH | 1 | +5.21 | +5.21 | 0/0/1 | nan | +4.89 | +5.52 | -8.8 |
| probcons | all protein | 10 | +0.90 | +0.45 | 4/1/5 | 0.275 | +1.28 | +0.51 | +0.4 |
| probcons | BAliBASE | 5 | -0.37 | -0.23 | 4/1/0 | 0.0625 | -0.01 | -0.72 | +0.8 |
| probcons | HomFam | 3 | +2.27 | +2.10 | 0/0/3 | 0.25 | +2.66 | +1.88 | +3.1 |
| probcons | SIMMOD | 1 | +0.95 | +0.95 | 0/0/1 | nan | +0.67 | +1.22 | -6.0 |
| probcons | SIMHIGH | 1 | +3.06 | +3.06 | 0/0/1 | nan | +4.24 | +1.87 | -3.6 |
| clustalo | all protein | 10 | +2.11 | +0.18 | 4/0/6 | 0.322 | +2.17 | +2.05 | -1.9 |
| clustalo | BAliBASE | 5 | -0.31 | -0.43 | 4/0/1 | 0.125 | -0.15 | -0.48 | +1.3 |
| clustalo | HomFam | 3 | +1.85 | +2.07 | 0/0/3 | 0.25 | +1.35 | +2.36 | +2.6 |
| clustalo | SIMMOD | 1 | +7.59 | +7.59 | 0/0/1 | nan | +8.55 | +6.63 | -16.7 |
| clustalo | SIMHIGH | 1 | +9.52 | +9.52 | 0/0/1 | nan | +9.84 | +9.21 | -16.6 |
| kalign | all protein | 10 | +1.44 | +0.35 | 3/0/7 | 0.375 | +4.42 | -1.54 | -2.4 |
| kalign | BAliBASE | 5 | +0.26 | +0.33 | 1/0/4 | 0.438 | +3.39 | -2.86 | +0.6 |
| kalign | HomFam | 3 | +1.74 | -1.43 | 2/0/1 | 1 | +5.51 | -2.04 | -5.4 |
| kalign | SIMMOD | 1 | +2.20 | +2.20 | 0/0/1 | nan | +2.97 | +1.44 | -3.1 |
| kalign | SIMHIGH | 1 | +5.67 | +5.67 | 0/0/1 | nan | +7.75 | +3.60 | -8.1 |

## Pre-registered primary test (protein sets, Holm over 4 candidates)

| candidate | n | mean d vs L-INS-i | W/T/L | raw p | Holm p |
|---|---|---|---|---|---|
| muscle5 | 10 | +0.83 | 4/0/6 | 0.232 | 0.697 |
| probcons | 10 | +0.90 | 4/1/5 | 0.275 | 0.697 |
| famsa | 10 | +1.26 | 2/0/8 | 0.0137 | 0.0547 |
| clustalo | 10 | +2.11 | 4/0/6 | 0.322 | 0.697 |

## Subset-level accuracy (pooled over subsets; HomFam only subsets with >= 2 seed sequences)

| set | linsi | ginsi | muscle5 | famsa | probcons | clustalo | kalign | prank |
|---|---|---|---|---|---|---|---|---|
| BAliBASE (n=5) | 10.43 | 10.48 | 10.14 | 10.77 | 10.21 | 10.06 | 10.91 |  |
| HomFam (n=3) | 6.04 | 5.53 | 8.31 | 4.31 | 8.54 | 7.63 | 3.41 |  |
| SIMMOD (n=1) | 3.50 | 3.27 | 4.47 | 6.40 | 4.56 | 10.33 | 5.59 |  |
| SIMHIGH (n=1) | 8.10 | 7.93 | 10.20 | 13.26 | 10.31 | 18.35 | 13.12 |  |
| nucleotide (n=0) |  |  |  |  |  |  |  |  |

## Subset-alignment cost: total CPU seconds over the 25 subsets (1 thread each); ratio to L-INS-i

| set | linsi | ginsi | muscle5 | famsa | probcons | clustalo | kalign | prank |
|---|---|---|---|---|---|---|---|---|
| BAliBASE | 81 (1x) | 77 (0.95x) | 317 (3.9x) | 24 (0.3x) | 627 (7.7x) | 45 (0.55x) | 17 (0.21x) |  |
| HomFam | 59 (1x) | 62 (1.1x) | 309 (5.2x) | 10 (0.17x) | 338 (5.7x) | 9 (0.15x) | 3 (0.045x) |  |
| SIMMOD | 165 (1x) | 181 (1.1x) | 1125 (6.8x) | 4 (0.027x) | 1626 (9.9x) | 28 (0.17x) | 10 (0.061x) |  |
| SIMHIGH | 245 (1x) | 253 (1x) | 1326 (5.4x) | 45 (0.18x) | 1706 (7x) | 51 (0.21x) | 10 (0.04x) |  |
| nucleotide |  |  |  |  |  |  |  |  |

GCM merge wall time (4 threads), L-INS-i subsets: median 28 s, range 1-54 s.

## Control: GCM merge on re-aligned L-INS-i subsets vs the original MAGUS run (SP error %)

| dataset | MAGUS (original) | merge-linsi |
|---|---|---|
| BBA0039 | 4.69 | 4.66 |
| BBA0067 | 26.28 | 26.40 |
| BBA0101 | 29.18 | 28.88 |
| BBA0154 | 21.66 | 21.76 |
| BBA0190 | 23.40 | 23.43 |
| HF_PDZ | 16.14 | 16.61 |
| HF_rrm | 18.57 | 18.41 |
| HF_zf-CCHH | 12.25 | 12.25 |
| SIMMOD_R1 | 12.71 | 12.79 |
| SIMHIGH_R1 | 25.04 | 25.00 |

## Whole-dataset baselines and PASTA variants (SP error %, wall s, 4 threads)

| dataset | MAGUS (L-INS-i subsets) | best MAGUS variant | pasta-mafft-it1 | pasta-prank-it1 | pasta-probcons-it1 |
|---|---|---|---|---|---|
| BBA0039 | 4.66 | 4.45 (clustalo) | 4.48 (542 s) | fail/timeout | 4.56 (1752 s) |

## Trees: FastTree -lg -gamma, normalized RF (%) vs the true tree

| dataset | merge_linsi | true |
|---|---|---|
| SIMMOD_R1 | 6.32 | 6.52 |
