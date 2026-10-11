## Part A tables (protein: 268 estimated alignments with trees, 33 draws)

### Within-draw (fixed effects) and paired-vs-MAGUS association with FastTree nRF

| measure | within r [95% CI] | within R² | within ρ | slope (RF pts per pt) | paired ΔRF~Δx ρ (p) | paired r | paired LODO-CV R² | across-draw ρ(RF) | ρ(RF − RF_true) |
|---|---|---|---|---|---|---|---|---|---|
| gap openings est/true | -0.21 [-0.38, -0.04] | 0.04 | -0.18 | -13.613 | -0.09 (0.1) | -0.11 | -0.00 | -0.21 | -0.17 |
| SP error (SPFN+SPFP)/2 | +0.21 [+0.03, +0.40] | 0.04 | +0.21 | +0.053 | +0.25 (9e-05) | +0.22 | +0.04 | +0.44 | +0.47 |
| 1 - column precision | +0.21 [+0.01, +0.39] | 0.04 | +0.23 | +0.032 | +0.19 (0.003) | +0.15 | +0.01 | +0.30 | +0.31 |
| PI columns est/true | +0.21 [+0.01, +0.40] | 0.04 | +0.20 | +1.394 | +0.26 (5e-05) | +0.22 | +0.04 | +0.03 | +0.05 |
| p-distance distortion, random pairs | +0.20 [+0.11, +0.29] | 0.04 | +0.21 | +0.800 | +0.20 (0.002) | +0.19 | +0.01 | +0.49 | +0.48 |
| error at indel boundaries | +0.20 [+0.01, +0.40] | 0.04 | +0.19 | +0.036 | +0.23 (0.0003) | +0.21 | +0.03 | +0.40 | +0.43 |
| residue error | +0.20 [+0.02, +0.40] | 0.04 | +0.20 | +0.038 | +0.23 (0.0003) | +0.21 | +0.03 | +0.39 | +0.43 |
| over-merged est. columns | +0.20 [+0.01, +0.40] | 0.04 | +0.20 | +0.030 | +0.18 (0.007) | +0.15 | +0.01 | +0.24 | +0.27 |
| Fitch excess on true tree | +0.20 [+0.00, +0.40] | 0.04 | +0.21 | +0.169 | +0.10 (0.1) | +0.11 | -0.00 | +0.36 | +0.43 |
| 1 - TC (column recall error) | +0.19 [+0.01, +0.38] | 0.04 | +0.21 | +0.035 | +0.16 (0.01) | +0.13 | +0.00 | +0.30 | +0.31 |
| split residues | +0.19 [-0.00, +0.39] | 0.04 | +0.18 | +0.038 | +0.23 (0.0004) | +0.20 | +0.02 | +0.37 | +0.40 |
| split residues, PI | +0.19 [+0.01, +0.39] | 0.04 | +0.18 | +0.034 | +0.23 (0.0005) | +0.20 | +0.02 | +0.37 | +0.40 |
| SPFN | +0.19 [+0.01, +0.39] | 0.04 | +0.19 | +0.024 | +0.22 (0.0005) | +0.20 | +0.02 | +0.37 | +0.40 |
| SPFN, PI columns | +0.19 [-0.01, +0.39] | 0.04 | +0.19 | +0.024 | +0.22 (0.0005) | +0.20 | +0.02 | +0.37 | +0.40 |
| Fitch excess, PI cols | +0.19 [+0.00, +0.39] | 0.04 | +0.20 | +0.211 | +0.11 (0.1) | +0.11 | -0.00 | +0.36 | +0.42 |
| per-taxon error, 90th pct | +0.19 [-0.00, +0.38] | 0.04 | +0.19 | +0.035 | +0.22 (0.0006) | +0.21 | +0.02 | +0.36 | +0.41 |
| length est/true | -0.19 [-0.39, +0.01] | 0.03 | -0.19 | -4.637 | -0.09 (0.2) | -0.08 | -0.01 | -0.23 | -0.29 |
| misplaced residues | +0.18 [+0.08, +0.28] | 0.03 | +0.19 | +0.306 | +0.12 (0.07) | +0.15 | -0.00 | +0.52 | +0.56 |
| over-merged PI est. columns | +0.17 [-0.03, +0.40] | 0.03 | +0.18 | +0.043 | +0.14 (0.03) | +0.10 | -0.01 | +0.17 | +0.20 |
| compared sites est/true | -0.17 [-0.38, +0.03] | 0.03 | -0.16 | -1.744 | -0.21 (0.001) | -0.18 | +0.01 | -0.14 | -0.14 |
| error away from indels | +0.16 [+0.01, +0.31] | 0.02 | +0.18 | +0.047 | +0.12 (0.07) | +0.09 | -0.00 | +0.28 | +0.28 |
| misplaced residues, PI | +0.16 [+0.05, +0.27] | 0.02 | +0.17 | +0.253 | +0.10 (0.1) | +0.14 | -0.01 | +0.51 | +0.56 |
| 1 - rank corr. of p-distances | +0.12 [-0.05, +0.31] | 0.01 | +0.11 | +8.238 | +0.17 (0.01) | +0.16 | -0.01 | +0.53 | +0.56 |
| SPFN between cherry taxa | +0.10 [-0.08, +0.29] | 0.01 | +0.11 | +1.984 | +0.09 (0.2) | +0.01 | -0.06 | +0.57 | +0.61 |
| SPFP, PI est. columns | +0.10 [-0.06, +0.24] | 0.01 | +0.08 | +0.106 | +0.10 (0.1) | +0.11 | -0.02 | +0.49 | +0.53 |
| SPFP | +0.10 [-0.06, +0.23] | 0.01 | +0.08 | +0.106 | +0.10 (0.1) | +0.11 | -0.02 | +0.49 | +0.53 |
| SPFP, random taxon pairs | +0.09 [-0.06, +0.23] | 0.01 | +0.09 | +0.104 | +0.10 (0.1) | +0.10 | -0.02 | +0.49 | +0.53 |
| p-distance distortion, cherries | +0.07 [-0.13, +0.22] | 0.00 | +0.01 | +3.325 | -0.26 (4e-05) | -0.17 | -0.04 | +0.57 | +0.66 |
| over-split true columns | +0.06 [-0.08, +0.24] | 0.00 | +0.15 | +0.312 | -0.18 (0.006) | -0.28 | +0.01 | +0.46 | +0.46 |
| SPFP between cherry taxa | -0.02 [-0.22, +0.14] | 0.00 | -0.08 | -0.442 | -0.20 (0.002) | -0.11 | -0.00 | +0.54 | +0.58 |
| over-split PI columns | +0.02 [-0.11, +0.18] | 0.00 | +0.08 | +0.028 | -0.16 (0.02) | -0.24 | +0.02 | +0.19 | +0.29 |
| p-distance bias, random pairs | -0.01 [-0.21, +0.17] | 0.00 | -0.01 | -0.015 | -0.06 (0.4) | -0.06 | -0.04 | +0.27 | +0.26 |

### Within-draw LODO cross-validated R² of small models

| predictors | CV R² |
|---|---|
| avgErr | +0.02 |
| SPFN | +0.02 |
| SPFP | -0.00 |
| SPFN + SPFP | +0.03 |
| split_res + misplaced_res | +0.02 |
| SPFN_pi + SPFP_pi | +0.03 |
| oversplit_cols + overmerge_cols | +0.00 |
| fitch_excess | +0.02 |
| fitch_excess_pi | +0.02 |
| SPFN + SPFP + fitch_excess | +0.03 |
| misplaced_res_pi | +0.02 |
| SPFP_pi | -0.00 |
| tax_err_p90 | +0.02 |
| TCerr | +0.02 |
| dist_mae_rand | +0.04 |
| dist_bias_rand | -0.02 |
| dist_rank_err | -0.00 |
| SPFN + SPFP + dist_mae_rand | +0.03 |
| sites_rand + dist_bias_rand | +0.01 |
| SPFN + SPFP + TCerr + len_ratio | +0.02 |

### Noise floor

- vote_hard vs its masked copy (3-8 of ~8,000 columns removed): n = 17, SD of ΔRF = 0.45, mean |ΔRF| = 0.22
- pairs of alignments of the same draw within 0.25 SP points (and 0.5 SPFN): n = 30, SD of ΔRF = 1.09
- SD of within-draw RF deviations: 0.75; SD of ΔRF vs MAGUS: 1.18
- implied ceiling on within-draw R² (if every alignment's RF carries independent noise of variance SD_mask²/2): 0.82

### Per-variant means, paired vs MAGUS (Δ in points)

| variant | n | ΔRF | ΔSP err | ΔSPFN | ΔSPFP | Δsplit res | Δmisplaced res | Δmisplaced PI | Δover-merged PI cols | ΔFitch excess (%) |
|---|---|---|---|---|---|---|---|---|---|---|
| es3 | 33 | -0.06 | -5.56 | -11.53 | +0.41 | -7.52 | -0.11 | +0.00 | -4.12 | -1.601 |
| es4 | 17 | -0.27 | -6.69 | -13.72 | +0.33 | -9.02 | -0.27 | -0.15 | -5.74 | -1.728 |
| hardf | 33 | -0.42 | -2.60 | -4.51 | -0.70 | -3.18 | -0.56 | -0.53 | -0.15 | -0.671 |
| recipe | 33 | -0.36 | -6.42 | -11.95 | -0.88 | -7.84 | -1.03 | -0.96 | -5.34 | -1.945 |
| split_magus | 5 | -1.10 | -8.50 | +0.00 | -17.01 | +0.00 | -10.56 | -10.74 | -77.98 | -13.879 |
| split_recipe | 5 | -1.89 | -15.39 | -13.77 | -17.01 | -9.25 | -10.56 | -10.74 | -77.98 | -13.916 |
| vote_hard | 17 | -0.18 | -6.55 | -12.11 | -0.99 | -8.25 | -1.01 | -0.98 | -3.31 | -1.334 |
| vote_hard-bb | 17 | -0.43 | -6.69 | -12.49 | -0.88 | -8.46 | -1.00 | -0.96 | -4.06 | -1.430 |
| vote_hard_mask | 17 | -0.29 | -6.55 | -12.11 | -0.99 | -8.23 | -1.02 | -0.98 | -3.37 | -1.399 |
| vote_soft | 17 | +0.32 | +0.18 | +1.57 | -1.20 | +0.63 | -0.62 | -0.73 | +3.15 | +0.590 |
| vote_soft-bb | 17 | -0.10 | +0.27 | +1.45 | -0.91 | +0.69 | -0.35 | -0.43 | +3.38 | +0.567 |
| vote_soft2 | 17 | +0.09 | -0.91 | -0.53 | -1.30 | -0.78 | -0.78 | -0.88 | +1.81 | +0.264 |
| vote_soft4 | 17 | -0.08 | -2.97 | -4.51 | -1.43 | -3.49 | -1.00 | -1.08 | -0.37 | -0.229 |

### Price of each error type (within-draw OLS, RF points per error point; 95% cluster-bootstrap CI)

| rows | model | coefficients | within R² | n (draws) |
|---|---|---|---|---|
| estimated only | avgErr | avgErr +0.053 [+0.007, +0.098] | 0.04 | 268 (33) |
| estimated only | SPFN + SPFP | SPFN +0.029 [+0.005, +0.050]; SPFP +0.164 [+0.020, +0.292] | 0.06 | 268 (33) |
| estimated only | split_res + misplaced_res | split_res +0.030 [-0.014, +0.070]; misplaced_res +0.228 [+0.004, +0.480] | 0.05 | 268 (33) |
| estimated + true alignment | avgErr | avgErr +0.132 [+0.102, +0.160] | 0.46 | 301 (33) |
| estimated + true alignment | SPFN + SPFP | SPFN +0.027 [+0.004, +0.050]; SPFP +0.144 [+0.089, +0.199] | 0.50 | 301 (33) |
| estimated + true alignment | split_res + misplaced_res | split_res +0.029 [-0.010, +0.068]; misplaced_res +0.261 [+0.159, +0.353] | 0.50 | 301 (33) |
| SIMHIGH R1-R5 (gcmtrees): estimated + true + oracle split(X) | avgErr | avgErr +0.163 [+0.126, +0.212] | 0.61 | 35 (5) |
| SIMHIGH R1-R5 (gcmtrees): estimated + true + oracle split(X) | SPFN + SPFP | SPFN +0.062 [+0.028, +0.110]; SPFP +0.106 [+0.093, +0.121] | 0.64 | 35 (5) |
| SIMHIGH R1-R5 (gcmtrees): estimated + true + oracle split(X) | split_res + misplaced_res | split_res +0.097 [+0.038, +0.177]; misplaced_res +0.180 [+0.157, +0.203] | 0.62 | 35 (5) |

### Observed ΔRF vs ΔRF predicted from error prices

| variant | n | observed ΔRF | predicted (prices from estimated alignments) | predicted (prices incl. true alignment) |
|---|---|---|---|---|
| es3 | 33 | -0.06 | -0.26 | -0.25 |
| es4 | 17 | -0.27 | -0.34 | -0.32 |
| hardf | 33 | -0.42 | -0.24 | -0.22 |
| recipe | 33 | -0.36 | -0.49 | -0.44 |
| split_magus | 5 | -1.10 | -2.78 | -2.45 |
| split_recipe | 5 | -1.89 | -3.18 | -2.82 |
| vote_hard | 17 | -0.18 | -0.51 | -0.47 |
| vote_hard-bb | 17 | -0.43 | -0.50 | -0.46 |
| vote_hard_mask | 17 | -0.29 | -0.51 | -0.47 |
| vote_soft | 17 | +0.32 | -0.15 | -0.13 |
| vote_soft-bb | 17 | -0.10 | -0.11 | -0.09 |
| vote_soft2 | 17 | +0.09 | -0.23 | -0.20 |
| vote_soft4 | 17 | -0.08 | -0.36 | -0.33 |
