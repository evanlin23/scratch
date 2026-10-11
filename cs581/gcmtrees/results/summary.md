## Per dataset: FastTree nRF (%) to the true tree and SP error

| dataset | RF true | RF MAGUS | room | RF recipe | RF es3 | RF hard | Δ recipe | Δ es3 | Δ hard | SP err MAGUS | Δ SP recipe | Δ SPFN / ΔSPFP recipe | Δ SP es3 | Δ SP hard |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SIMHIGH_R1 | 6.12 | 11.94 | +5.82 | 12.34 | 11.63 | 12.34 | +0.40 | -0.31 | +0.40 | 25.32 | -5.32 | -9.66 / -0.98 | -4.49 | -0.39 |
| SIMHIGH_R2 | 5.22 | 9.33 | +4.11 | 7.52 | 8.73 | 7.82 | -1.81 | -0.60 | -1.51 | 22.91 | -8.63 | -17.65 / +0.40 | -7.52 | -3.16 |
| SIMHIGH_R3 | 6.42 | 10.03 | +3.61 | 9.23 | 9.83 | 8.93 | -0.80 | -0.20 | -1.10 | 27.55 | -9.47 | -17.97 / -0.97 | -7.95 | -5.60 |
| SIMHIGH_R4 | 6.62 | 9.63 | +3.01 | 10.53 | 11.94 | 8.63 | +0.90 | +2.31 | -1.00 | 22.95 | -5.52 | -10.29 / -0.74 | -5.14 | -2.68 |
| SIMHIGH_R5 | 6.52 | 8.63 | +2.11 | 9.73 | 9.53 | 9.13 | +1.10 | +0.90 | +0.50 | 20.62 | -6.63 | -13.27 / +0.01 | -6.04 | -1.72 |
| SIMHIGH_R6 | 9.73 | 12.14 | +2.41 | 11.43 | 11.74 | 11.84 | -0.71 | -0.40 | -0.30 | 20.63 | -7.83 | -15.56 / -0.10 | -6.93 | -3.68 |
| SIMHIGH_R7 | 6.62 | 10.63 | +4.01 | 8.32 | 8.53 | 10.63 | -2.31 | -2.10 | +0.00 | 27.28 | -8.14 | -14.51 / -1.77 | -6.14 | -3.53 |
| SIMHIGH_R8 | 8.22 | 11.23 | +3.01 | 10.53 | 11.03 | 11.03 | -0.70 | -0.20 | -0.20 | 26.19 | -10.14 | -19.18 / -1.10 | -9.42 | -4.02 |
| SIMMOD_R1 | 6.52 | 5.92 | -0.60 | 6.62 | 5.82 | 6.22 | +0.70 | -0.10 | +0.30 | 12.50 | -2.39 | -4.11 / -0.68 | -3.17 | -0.36 |
| SIMMOD_R2 | 6.52 | 7.02 | +0.50 | 7.02 | 6.52 | 6.92 | +0.00 | -0.50 | -0.10 | 10.15 | -4.52 | -8.74 / -0.30 | -4.41 | -2.82 |
| SIMMOD_R3 | 6.12 | 7.32 | +1.20 | 7.12 | 7.22 | 7.02 | -0.20 | -0.10 | -0.30 | 11.32 | -4.22 | -7.80 / -0.63 | -3.40 | -2.25 |
| SIMMOD_R4 | 5.52 | 6.22 | +0.70 | 6.92 | 6.82 | 6.62 | +0.70 | +0.60 | +0.40 | 10.69 | -3.82 | -7.52 / -0.12 | -3.10 | -1.02 |
| SIMMOD_R5 | 8.53 | 8.02 | -0.51 | 8.43 | 7.92 | 7.52 | +0.41 | -0.10 | -0.50 | 9.88 | -4.41 | -8.65 / -0.18 | -3.87 | -2.62 |
| SIMMOD_R6 | 7.92 | 8.32 | +0.40 | 7.82 | 8.02 | 8.63 | -0.50 | -0.30 | +0.31 | 9.81 | -2.52 | -5.07 / +0.03 | -2.80 | -1.02 |
| SIMMOD_R7 | 6.82 | 6.92 | +0.10 | 6.62 | 6.82 | 6.52 | -0.30 | -0.10 | -0.40 | 11.49 | -4.72 | -8.19 / -1.25 | -3.71 | -2.62 |
| SIMMOD_R8 | 8.53 | 8.02 | -0.51 | 8.12 | 8.02 | 8.22 | +0.10 | +0.00 | +0.20 | 9.55 | -2.71 | -4.49 / -0.93 | -2.60 | -1.94 |

## Paired summaries (method − MAGUS; RF in nRF points, SP in error points; W/T/L tie band 0.1 RF point)

| subset | method | n | mean Δ RF | median Δ RF | W/T/L | Wilcoxon p | mean Δ FN | mean Δ FP | mean Δ SP err | mean ΔSPFN | mean ΔSPFP | mean room |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| all | recipe | 16 | -0.19 | -0.10 | 8/2/6 | 0.629 | -0.19 | -0.19 | -5.69 | -10.79 | -0.58 | +1.84 |
| all | es3 | 16 | -0.07 | -0.15 | 8/5/3 | 0.256 | -0.07 | -0.07 | -5.04 | -10.53 | +0.44 | +1.84 |
| all | hard | 16 | -0.21 | -0.15 | 8/2/6 | 0.349 | -0.21 | -0.21 | -2.46 | -4.48 | -0.44 | +1.84 |
| all | true aln | 16 | -1.84 | -1.66 | 12/1/3 | 0.00612 | | | | | | |
| SIMHIGH | recipe | 8 | -0.49 | -0.70 | 5/0/3 | 0.461 | -0.49 | -0.49 | -7.71 | -14.76 | -0.66 | +3.51 |
| SIMHIGH | es3 | 8 | -0.07 | -0.26 | 6/0/2 | 0.617 | -0.07 | -0.07 | -6.70 | -14.05 | +0.64 | +3.51 |
| SIMHIGH | hard | 8 | -0.40 | -0.25 | 5/1/2 | 0.297 | -0.40 | -0.40 | -3.10 | -5.57 | -0.62 | +3.51 |
| SIMHIGH | true aln | 8 | -3.51 | -3.31 | 8/0/0 | 0.00781 | | | | | | |
| SIMMOD | recipe | 8 | +0.11 | +0.05 | 3/2/3 | 0.578 | +0.11 | +0.11 | -3.66 | -6.82 | -0.51 | +0.16 |
| SIMMOD | es3 | 8 | -0.07 | -0.10 | 2/5/1 | 0.297 | -0.07 | -0.07 | -3.38 | -7.01 | +0.25 | +0.16 |
| SIMMOD | hard | 8 | -0.01 | +0.05 | 3/1/4 | 0.844 | -0.01 | -0.01 | -1.83 | -3.40 | -0.27 | +0.16 |
| SIMMOD | true aln | 8 | -0.16 | -0.25 | 4/1/3 | 0.711 | | | | | | |

## Does Δ RF track Δ SPFN or Δ SPFP? (Spearman over dataset × method, n = 48)

- Δ RF vs Δ SPFN: rho = +0.32 (p = 0.0267)
- Δ RF vs Δ SPFP: rho = +0.03 (p = 0.847)
- Δ RF vs Δ SP error: rho = +0.34 (p = 0.0188)
- room (RF MAGUS − RF true) vs MAGUS SPFN: rho = +0.83 (p = 8e-05, n = 16)
- room (RF MAGUS − RF true) vs MAGUS SPFP: rho = +0.82 (p = 0.000116, n = 16)

## IQ-TREE 3 (-m LG+G4 --fast) nRF (%)

| dataset | true | magus | recipe | es3 | hard | Δ recipe | Δ es3 | Δ hard |
|---|---|---|---|---|---|---|---|---|
| SIMHIGH_R1 | 6.02 | 13.34 | 12.84 |  |  | -0.50 |  |  |
| SIMHIGH_R2 | 5.12 | 9.83 | 9.03 |  |  | -0.80 |  |  |
| SIMHIGH_R3 |  | 10.13 | 10.03 |  |  | -0.10 |  |  |
| SIMHIGH_R6 | 7.82 | 10.53 |  |  |  |  |  |  |
| SIMHIGH_R7 | 6.82 | 10.23 |  |  |  |  |  |  |

- IQ-TREE recipe − MAGUS: n = 3, mean -0.47, W/T/L 2/1/0, p = 0.25

## Why-not diagnostic: zero-SPFP refinements (FastTree nRF %)

split(X) = common refinement of the true alignment and X: X's true-positive pairs only (X's SPFN, SPFP = 0).

| dataset | true | split(MAGUS) | MAGUS | split(recipe) | recipe | FP cost MAGUS | FP cost recipe | split(recipe) − split(MAGUS) |
|---|---|---|---|---|---|---|---|---|
| SIMHIGH_R1 | 6.12 | 12.24 | 11.94 | 8.32 | 12.34 | -0.30 | +4.02 | -3.92 |
| SIMHIGH_R2 | 5.22 | 6.32 | 9.33 | 8.93 | 7.52 | +3.01 | -1.41 | +2.61 |
| SIMHIGH_R3 | 6.42 | 8.93 | 10.03 | 7.32 | 9.23 | +1.10 | +1.91 | -1.61 |
| SIMHIGH_R4 | 6.62 | 8.43 | 9.63 | 8.22 | 10.53 | +1.20 | +2.31 | -0.21 |
| SIMHIGH_R5 | 6.52 | 8.12 | 8.63 | 7.32 | 9.73 | +0.51 | +2.41 | -0.80 |
| mean | | | | | | +1.10 | +1.85 | -0.79 |

## Where the alignment changes land (true-alignment columns)

Gains = recovered true pairs as % of all true pairs. 'dense' = true columns holding >= 50% of taxa; 'informative' = parsimony-informative columns. Residue level: misplaced = not in the largest true-column group of its estimated column; split = outside the largest estimated fragment of its true column (% of residues).

| dataset | pairs in dense cols | MAGUS recall dense / gappy | recipe gain dense / gappy | es3 gain dense / gappy | gain in informative cols (recipe) | misplaced res MAGUS / recipe / es3 | of which in blocks >= 10 (MAGUS / recipe) | split res MAGUS / recipe / es3 |
|---|---|---|---|---|---|---|---|---|
| SIMHIGH_R1 | 97.4% | 68.7% / 73.6% | +9.32 / +0.34 | +9.17 / +0.32 | +9.66 | 12.2% / 10.9% / 12.1% | 7.2% / 6.5% | 18.8% / 12.2% / 12.3% |
| SIMHIGH_R2 | 98.6% | 68.7% / 70.7% | +17.43 / +0.22 | +15.63 / +0.15 | +17.65 | 9.6% / 8.9% / 9.4% | 4.7% / 4.6% | 19.2% / 7.4% / 8.9% |
| SIMHIGH_R3 | 97.7% | 63.3% / 73.7% | +17.64 / +0.33 | +16.88 / +0.28 | +17.97 | 11.2% / 10.2% / 11.5% | 6.8% / 6.3% | 21.4% / 9.7% / 10.3% |
| SIMHIGH_R4 | 98.8% | 72.3% / 75.9% | +10.18 / +0.11 | +10.29 / +0.05 | +10.29 | 10.9% / 10.1% / 10.8% | 6.7% / 6.3% | 16.2% / 9.2% / 9.4% |
| SIMHIGH_R5 | 98.0% | 72.8% / 75.3% | +13.04 / +0.23 | +12.25 / +0.18 | +13.27 | 8.9% / 8.3% / 8.8% | 4.6% / 4.5% | 16.5% / 7.4% / 8.1% |

## MAGUS runtime

- SIMHIGH_R1 draw 0: MAGUS wall 1089.6 s
- SIMHIGH_R2 draw 0: MAGUS wall 1299.9 s
- SIMHIGH_R3 draw 0: MAGUS wall 2236.7 s
- SIMHIGH_R4 draw 0: MAGUS wall 2416.1 s
- SIMHIGH_R5 draw 0: MAGUS wall 3081.6 s
- SIMHIGH_R6 draw 0: MAGUS wall 1051.6 s
- SIMHIGH_R7 draw 0: MAGUS wall 1218.3 s
- SIMHIGH_R8 draw 0: MAGUS wall 1272.9 s
- SIMMOD_R1 draw 0: MAGUS wall 794.3 s
- SIMMOD_R2 draw 0: MAGUS wall 842.4 s
- SIMMOD_R3 draw 0: MAGUS wall 730.3 s
- SIMMOD_R4 draw 0: MAGUS wall 752.6 s
- SIMMOD_R5 draw 0: MAGUS wall 557.4 s
- SIMMOD_R6 draw 0: MAGUS wall 543.8 s
- SIMMOD_R7 draw 0: MAGUS wall 671.4 s
- SIMMOD_R8 draw 0: MAGUS wall 631.1 s
