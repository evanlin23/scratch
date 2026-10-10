### Pooled simulated data (FastMulRFS + DISCO, estimated gene trees)

#### Primary comparisons (Holm over 4)

| ASTRID-Pro vs | n | mean diff (FN rate) | 95% CI | W/T/L | p (Wilcoxon) | Holm p |
|---|---|---|---|---|---|---|
| astral-pro3 | 140 | -0.0135 | [-0.0179, -0.0092] | 84/30/26 | 1.4e-08 | 5.6e-08 |
| astrid-multi | 200 | -0.0034 | [-0.0060, -0.0009] | 75/72/53 | 0.023 | 0.068 |
| astrid-disco | 199 | +0.0009 | [-0.0013, +0.0031] | 49/86/64 | 0.32 | 0.64 |
| asteroid | 188 | +0.0002 | [-0.0024, +0.0027] | 58/68/62 | 0.98 | 0.98 |

#### Secondary comparisons (no correction)

| ASTRID-Pro vs | n | mean diff (FN rate) | 95% CI | W/T/L | p (Wilcoxon) |
|---|---|---|---|---|---|
| astrid-pro-r0 | 200 | +0.0003 | [-0.0003, +0.0010] | 9/179/12 | 0.75 |
| astrid-pro-s | 200 | -0.0008 | [-0.0019, +0.0002] | 24/162/14 | 0.23 |
| disco-astral | 121 | -0.0235 | [-0.0302, -0.0172] | 79/22/20 | 7.5e-11 |
| fastmulrfs | 121 | -0.0240 | [-0.0299, -0.0184] | 87/18/16 | 7.8e-13 |
| wqfm-gdl | 54 | -0.0042 | [-0.0088, +0.0002] | 22/14/18 | 0.092 |
| duploss2 | 54 | -0.1075 | [-0.1275, -0.0888] | 53/1/0 | 2.3e-10 |

### FastMulRFS data

#### FastMulRFS data: ASTRID-Pro vs each method

| ASTRID-Pro vs | n | mean diff (FN rate) | 95% CI | W/T/L | p (Wilcoxon) |
|---|---|---|---|---|---|
| astrid-pro-r0 | 120 | +0.0003 | [-0.0007, +0.0014] | 8/103/9 | 0.78 |
| astrid-pro-s | 120 | -0.0014 | [-0.0032, +0.0003] | 24/82/14 | 0.23 |
| astrid-multi | 120 | -0.0024 | [-0.0059, +0.0012] | 43/40/37 | 0.13 |
| astrid-disco | 120 | +0.0018 | [-0.0013, +0.0050] | 31/46/43 | 0.3 |
| asteroid | 120 | +0.0013 | [-0.0022, +0.0049] | 39/36/45 | 0.59 |
| astral-pro3 | 120 | -0.0150 | [-0.0201, -0.0101] | 76/21/23 | 5.2e-08 |
| disco-astral | 119 | -0.0237 | [-0.0307, -0.0172] | 77/22/20 | 7.9e-11 |
| fastmulrfs | 119 | -0.0246 | [-0.0306, -0.0189] | 87/17/15 | 4.9e-13 |
| wqfm-gdl | 54 | -0.0042 | [-0.0090, +0.0004] | 22/14/18 | 0.092 |
| duploss2 | 54 | -0.1075 | [-0.1279, -0.0878] | 53/1/0 | 2.3e-10 |

### DISCO data

#### DISCO data: ASTRID-Pro vs each method

| ASTRID-Pro vs | n | mean diff (FN rate) | 95% CI | W/T/L | p (Wilcoxon) |
|---|---|---|---|---|---|
| astrid-pro-r0 | 80 | +0.0003 | [-0.0003, +0.0008] | 1/76/3 | 0.75 |
| astrid-pro-s | 80 | +0.0000 | [+0.0000, +0.0000] | 0/80/0 | 1 |
| astrid-multi | 80 | -0.0048 | [-0.0085, -0.0015] | 32/32/16 | 0.016 |
| astrid-disco | 79 | -0.0005 | [-0.0036, +0.0022] | 18/40/21 | 0.83 |
| asteroid | 68 | -0.0017 | [-0.0045, +0.0011] | 19/32/17 | 0.27 |
| astral-pro3 | 20 | -0.0041 | [-0.0087, -0.0000] | 8/9/3 | 0.054 |
| disco-astral | 2 | -0.0102 | [-0.0102, -0.0102] | 2/0/0 | 0.5 |
| fastmulrfs | 2 | +0.0102 | [+0.0000, +0.0204] | 0/1/1 | 1 |

#### FastMulRFS data: mean FN rate by sequence length and #genes

| sqln | ngen | n | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl | duploss2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 100 | 60 | 0.042 | 0.042 | 0.042 | 0.046 | 0.041 | 0.040 | 0.048 | 0.050 (59) | 0.053 | **0.040 (27)** | 0.087 (27) |
| 25 | 100 | 60 | 0.106 | 0.105 | 0.109 | 0.107 | **0.104** | 0.105 | 0.130 | 0.145 | 0.144 (59) | 0.114 (27) | 0.273 (27) |

#### FastMulRFS data: mean FN rate by condition (all sqln, ngen)

| cond | n | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl | duploss2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dl1e-10_ps1e7 | 20 | 0.047 | 0.046 | 0.047 | 0.050 | 0.047 | 0.046 | 0.057 | 0.060 | 0.058 | **0.043 (10)** | 0.156 (10) |
| dl1e-10_ps5e7 | 20 | 0.069 | 0.069 | 0.067 | 0.068 | **0.066** | 0.070 | 0.070 | 0.074 | 0.081 | 0.068 (10) | 0.193 (10) |
| dl2e-10_ps1e7 | 20 | 0.065 | 0.064 | 0.069 | 0.065 | 0.064 | **0.063** | 0.076 | 0.082 | 0.089 | 0.069 (10) | 0.178 (10) |
| dl2e-10_ps5e7 | 20 | 0.070 | 0.070 | 0.071 | 0.069 | 0.069 | **0.066** | 0.084 | 0.087 (19) | 0.091 (19) | 0.079 (8) | 0.197 (8) |
| dl5e-10_ps1e7 | 20 | 0.090 | 0.092 | 0.094 | 0.099 | 0.093 | **0.086** | 0.121 | 0.136 | 0.128 | 0.110 (8) | 0.164 (8) |
| dl5e-10_ps5e7 | 20 | 0.103 | 0.102 | 0.105 | 0.108 | **0.094** | 0.105 | 0.126 | 0.148 | 0.142 | 0.104 (8) | 0.196 (8) |

#### DISCO data: mean FN rate by condition (100 bp, 1000 genes)

| cond | n | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs |
|---|---|---|---|---|---|---|---|---|---|---|
| default | 7 | 0.042 | 0.042 | 0.042 | 0.045 | 0.036 | 0.042 | **0.036 (2)** | 0.071 (1) | 0.061 (1) |
| gdl_1e-10_0 | 7 | 0.041 | 0.041 | 0.041 | 0.042 | 0.042 | **0.038** | 0.082 (2) | – | – |
| gdl_1e-10_05 | 7 | 0.044 | 0.042 | 0.044 | 0.044 | 0.044 | **0.041** | 0.046 (2) | – | – |
| gdl_1e-10_1 | 7 | 0.028 | 0.028 | 0.028 | **0.025** | 0.028 | 0.026 | 0.051 (2) | 0.071 (1) | 0.041 (1) |
| gdl_1e-9_0 | 6 | **0.058** | **0.058** | **0.058** | 0.092 | 0.077 | – | – | – | – |
| gdl_1e-9_05 | 7 | 0.031 | 0.029 | 0.031 | 0.039 | 0.029 (6) | 0.037 (3) | **0.000 (1)** | – | – |
| gdl_1e-9_1 | 7 | 0.034 | 0.034 | 0.034 | 0.045 | **0.031** | 0.047 | 0.051 (2) | – | – |
| gdl_5e-10_0 | 7 | 0.022 | 0.022 | 0.022 | 0.028 | 0.026 | 0.024 (5) | **0.020 (1)** | – | – |
| gdl_5e-10_05 | 7 | **0.054** | **0.054** | **0.054** | 0.060 | **0.054** | 0.063 | 0.056 (2) | – | – |
| ils_1e4 | 6 | 0.037 | 0.037 | 0.037 | **0.032** | 0.032 | 0.034 | 0.036 (2) | – | – |
| ils_2e8 | 6 | 0.053 | 0.053 | 0.053 | 0.054 | **0.048** | 0.048 | 0.051 (2) | – | – |
| missing_1000 | 6 | 0.039 | 0.039 | 0.039 | **0.034** | 0.043 | **0.034** | 0.041 (2) | – | – |

#### FastMulRFS data: mean runtime (s, 1 thread)

| ngen | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl | duploss2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 0.15 | 0.16 | 0.53 | 0.17 | 1.01 | 1.86 | 10.02 | 14.82 | 20.20 | 73.70 | 11.60 |

#### DISCO data: mean runtime (s, 1 thread)

| cond | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs |
|---|---|---|---|---|---|---|---|---|---|
| default | 0.82 | 0.83 | 3.07 | 0.89 | 8.07 | 11.01 | 384.71 | 891.78 | 404.98 |
| gdl_1e-10_0 | 0.83 | 0.85 | 3.11 | 1.03 | 6.91 | 1.03 | 556.47 | – | – |
| gdl_1e-10_05 | 0.86 | 0.88 | 2.89 | 0.88 | 6.18 | 8.93 | 458.58 | – | – |
| gdl_1e-10_1 | 0.69 | 0.76 | 2.75 | 0.90 | 5.67 | 11.91 | 308.91 | 636.33 | 726.43 |
| gdl_1e-9_0 | 41.12 | 39.98 | 168.97 | 162.71 | 294.54 | – | – | – | – |
| gdl_1e-9_05 | 8.72 | 8.91 | 40.42 | 66.07 | 78.14 | 25.02 | 586.02 | – | – |
| gdl_1e-9_1 | 1.36 | 1.39 | 5.05 | 3.43 | 11.93 | 13.84 | 565.49 | – | – |
| gdl_5e-10_0 | 3.61 | 3.60 | 13.10 | 9.53 | 31.54 | 12.39 | 783.47 | – | – |
| gdl_5e-10_05 | 1.49 | 1.54 | 5.61 | 2.84 | 14.54 | 16.64 | 1028.41 | – | – |
| ils_1e4 | 0.94 | 0.93 | 3.57 | 0.91 | 8.49 | 10.19 | 330.15 | – | – |
| ils_2e8 | 0.88 | 0.92 | 3.51 | 1.12 | 8.81 | 11.23 | 355.29 | – | – |
| missing_1000 | 0.50 | 0.53 | 1.86 | 0.67 | 4.37 | 3.75 | 129.99 | – | – |

### Scaling

#### Runtime (s) vs #species (species_1000 rep 01, ~400-435 genes)

| cond | astrid-pro | astrid-pro-s | astrid-multi | astrid-disco | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl |
|---|---|---|---|---|---|---|---|---|
| taxa_100 | 0.59 | 1.74 | 1.35 | 5.02 | 84.60 | 81.10 | 36.27 | 349.69 |
| taxa_1000 | 14.32 | 90.56 | 135.73 | 670.66 | – | – | – | – |
| taxa_200 | 1.49 | 6.19 | 2.61 | 9.17 | 541.56 | 485.31 | 175.77 | – |
| taxa_50 | 0.23 | 0.58 | 0.29 | 1.74 | 14.08 | 16.23 | 11.52 | 66.64 |
| taxa_500 | 3.63 | 26.62 | 42.11 | 48.26 | – | – | – | – |

#### FN rate vs #species

| cond | n | astrid-pro | astrid-pro-s | astrid-multi | astrid-disco | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl |
|---|---|---|---|---|---|---|---|---|---|
| taxa_100 | 1 | 0.082 | 0.082 | 0.093 | 0.082 | 0.082 | 0.093 | **0.072** | **0.072** |
| taxa_1000 | 1 | **0.066** | **0.066** | 0.085 | 0.077 | – | – | – | – |
| taxa_200 | 1 | **0.061** | **0.061** | 0.086 | 0.076 | 0.071 | 0.081 | 0.091 | – |
| taxa_50 | 1 | 0.085 | 0.085 | 0.106 | 0.064 | **0.043** | **0.043** | **0.043** | 0.064 |
| taxa_500 | 1 | 0.064 | 0.064 | 0.068 | **0.060** | – | – | – | – |

#### Runtime (s) vs #genes (gtrees_10000_l1, 100 species)

| ngen | astrid-pro | astrid-pro-s | astrid-multi | astrid-disco | astral-pro3 | fastmulrfs |
|---|---|---|---|---|---|---|
| 100 | 0.17 | 0.82 | 0.23 | 1.11 | 12.48 | 9.27 |
| 300 | 0.37 | 1.54 | 0.64 | 3.55 | 53.59 | 29.05 |
| 1000 | 1.18 | 4.47 | 1.27 | 9.51 | 340.44 | 136.85 |
| 3000 | 2.82 | 11.33 | 3.14 | 28.37 | 904.49 | – |
| 10000 | 8.00 | 36.16 | 9.18 | 95.29 | – | – |

#### FN rate vs #genes

| ngen | n | astrid-pro | astrid-pro-s | astrid-multi | astrid-disco | astral-pro3 | fastmulrfs |
|---|---|---|---|---|---|---|---|
| 100 | 2 | **0.051** | **0.051** | 0.056 | 0.071 | 0.066 | 0.168 |
| 1000 | 2 | 0.036 | 0.036 | 0.041 | **0.031** | 0.041 | 0.046 |
| 10000 | 2 | **0.031** | **0.031** | 0.036 | **0.031** | – | – |
| 300 | 2 | 0.041 | 0.041 | **0.036** | 0.051 | 0.061 | 0.087 |
| 3000 | 2 | 0.036 | 0.036 | 0.031 | 0.036 | **0.020 (1)** | – |

### Failed runs (timeouts and crashes; excluded from all paired tables)

| file | method | condition | error | n |
|---|---|---|---|---|
| disco_runs.jsonl | asteroid | gdl_1e-9_0 | CalledProcessError | 6 |
| disco_runs.jsonl | asteroid | gdl_1e-9_05 | CalledProcessError | 3 |
| disco_runs.jsonl | asteroid | gdl_5e-10_0 | CalledProcessError | 1 |
| disco_runs.jsonl | astral-pro3 | gdl_1e-9_0 | TimeoutExpired | 2 |
| disco_runs.jsonl | astral-pro3 | gdl_1e-9_05 | TimeoutExpired | 1 |
| disco_runs.jsonl | astral-pro3 | gdl_5e-10_0 | TimeoutExpired | 1 |
| disco_runs.jsonl | wqfm-gdl | gdl_1e-10_1 | CalledProcessError | 1 |
| empirical_runs.jsonl | asteroid | 1kp_c12 | CalledProcessError | 1 |
| empirical_runs.jsonl | asteroid | vertebrates188 | CalledProcessError | 1 |
| empirical_runs.jsonl | fastmulrfs | 1kp_c12 | CalledProcessError | 1 |
| empirical_runs.jsonl | wqfm-gdl | 1kp_c12 | TimeoutExpired | 1 |
| fmrfs_runs.jsonl | disco-astral | dl2e-10_ps5e7 | CalledProcessError | 1 |
| fmrfs_runs.jsonl | fastmulrfs | dl2e-10_ps5e7 | CalledProcessError | 1 |
| scaling_runs.jsonl | asteroid | taxa_200 | CalledProcessError | 1 |
| scaling_runs.jsonl | astral-pro3 | gtrees_10000_l1 | TimeoutExpired | 2 |
| scaling_runs.jsonl | astral-pro3 | taxa_1000 | TimeoutExpired | 1 |
| scaling_runs.jsonl | astral-pro3 | taxa_500 | TimeoutExpired | 1 |
| scaling_runs.jsonl | fastmulrfs | gtrees_10000_l1 | TimeoutExpired | 2 |
| scaling_runs.jsonl | fastmulrfs | taxa_1000 | TimeoutExpired | 1 |
| scaling_runs.jsonl | fastmulrfs | taxa_500 | TimeoutExpired | 1 |

### Empirical

#### FN rate vs reference tree

| cond | n | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1kp_c12 | 1 | **0.072** | **0.072** | 0.087 | 0.145 | 0.087 | – | **0.072** | **0.072** | – | – |
| fungi16 | 1 | **0.077** | **0.077** | **0.077** | **0.077** | **0.077** | **0.077** | **0.077** | **0.077** | **0.077** | **0.077** |
| vertebrates188 | 1 | **0.038** | **0.038** | **0.038** | **0.038** | **0.038** | – | – | – | – | – |

#### Runtime (s)

| cond | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl |
|---|---|---|---|---|---|---|---|---|---|---|
| 1kp_c12 | 6.89 | 6.69 | 26.28 | 6.83 | 46.92 | – | 3347.17 | 2962.00 | – | – |
| fungi16 | 0.53 | 0.53 | 1.22 | 0.69 | 2.88 | 6.02 | 11.55 | 23.05 | 8.94 | 17.07 |
| vertebrates188 | 21.51 | 22.54 | 79.32 | 22.98 | 168.53 | – | – | – | – | – |
