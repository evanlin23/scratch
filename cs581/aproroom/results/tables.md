## Q1. Error decomposition: FastMulRFS data (100 taxa, 100 genes, 6 conditions x 10 reps)

| method | L0 true tree+tags | L1 true tree | L2 est 100bp | L2 est 25bp | tagging (L1-L0) | GTEE 100bp (L2-L1) | GTEE 25bp |
|---|---|---|---|---|---|---|---|
| astrid-pro | 0.009 | 0.008 | 0.042 | 0.106 | -0.001 | 0.033 | 0.098 |
| astrid-disco | 0.009 | 0.009 | 0.041 | 0.104 | 0.000 | 0.031 | 0.095 |
| disco-astral | 0.009 | 0.009 | 0.051 | 0.145 | -0.000 | 0.042 | 0.136 |
| astrid-multi | – | 0.012 | 0.046 | 0.107 | – | 0.035 | 0.095 |
| asteroid | – | 0.010 | 0.040 | 0.105 | – | 0.030 | 0.095 |
| astral-pro3 | – | 0.008 | 0.048 | 0.130 | – | 0.040 | 0.122 |
| wqfm-gdl | – | 0.007 | 0.040 | 0.116 | – | 0.033 | 0.109 |

Mean FN rate by condition and level:

| level | cond | n | astrid-pro-tt | astrid-pro | astrid-disco-tt | astrid-disco | disco-astral-tt | disco-astral | astrid-multi | asteroid | astral-pro3 | wqfm-gdl |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| est100 | dl1e-10_ps1e7 | 10 | – | 0.026 | – | 0.029 | – | 0.034 | 0.027 | 0.028 | 0.033 | **0.023 (5)** |
| est100 | dl1e-10_ps5e7 | 10 | – | 0.042 | – | **0.039** | – | 0.047 | 0.044 | 0.040 | 0.045 | 0.047 (5) |
| est100 | dl2e-10_ps1e7 | 10 | – | 0.027 | – | 0.027 | – | 0.038 | 0.030 | 0.026 | 0.032 | **0.025 (5)** |
| est100 | dl2e-10_ps5e7 | 10 | – | 0.042 | – | 0.041 | – | 0.047 | 0.040 | **0.039** | 0.053 | 0.045 (5) |
| est100 | dl5e-10_ps1e7 | 10 | – | 0.051 | – | 0.054 | – | 0.063 | 0.063 | **0.041** | 0.054 | 0.045 (5) |
| est100 | dl5e-10_ps5e7 | 10 | – | 0.062 | – | **0.054** | – | 0.074 | 0.073 | 0.068 | 0.073 | 0.054 (5) |
| est25 | dl1e-10_ps1e7 | 10 | – | 0.068 | – | 0.065 | – | 0.086 | 0.073 | **0.064** | 0.080 | **0.064 (5)** |
| est25 | dl1e-10_ps5e7 | 10 | – | 0.095 | – | 0.094 | – | 0.100 | 0.091 | 0.099 | 0.095 | **0.089 (5)** |
| est25 | dl2e-10_ps1e7 | 10 | – | 0.103 | – | 0.101 | – | 0.126 | 0.100 | **0.100** | 0.120 | 0.113 (5) |
| est25 | dl2e-10_ps5e7 | 10 | – | 0.098 | – | 0.096 | – | 0.126 | 0.098 | **0.094** | 0.115 | 0.113 (5) |
| est25 | dl5e-10_ps1e7 | 10 | – | **0.130** | – | 0.133 | – | 0.208 | 0.135 | 0.131 | 0.189 | 0.165 (5) |
| est25 | dl5e-10_ps5e7 | 10 | – | 0.144 | – | **0.134** | – | 0.222 | 0.142 | 0.142 | 0.179 | 0.151 (5) |
| true | dl1e-10_ps1e7 | 10 | 0.001 | 0.001 | 0.001 | 0.001 | 0.001 | **0.000** | 0.001 | 0.001 | **0.000** | 0.002 (5) |
| true | dl1e-10_ps5e7 | 10 | 0.016 | 0.015 | 0.013 | 0.013 | 0.014 | 0.015 | 0.019 | 0.016 | 0.014 | **0.006 (5)** |
| true | dl2e-10_ps1e7 | 10 | **0.003** | **0.003** | **0.003** | **0.003** | 0.005 | 0.005 | 0.004 | **0.003** | 0.004 | 0.004 (5) |
| true | dl2e-10_ps5e7 | 10 | 0.013 | 0.013 | 0.015 | 0.018 | 0.013 | 0.014 | 0.013 | 0.016 | 0.013 | **0.012 (5)** |
| true | dl5e-10_ps1e7 | 10 | 0.003 | **0.001** | 0.003 | **0.001** | 0.003 | 0.003 | 0.004 | 0.003 | **0.001** | 0.002 (5) |
| true | dl5e-10_ps5e7 | 10 | 0.020 | 0.015 | 0.016 | 0.019 | 0.018 | 0.015 | 0.029 | 0.021 | 0.015 | **0.014 (5)** |

**FastMulRFS level true: ASTRID-Pro minus method**

| astrid-pro vs | n | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| astrid-disco | 60 | -0.0009 | 10/45/5 | 0.2 | 0.79 |
| asteroid | 60 | -0.0019 | 11/48/1 | 0.0059 | 0.035 |
| astral-pro3 | 60 | +0.0002 | 6/47/7 | 0.96 | 1 |
| wqfm-gdl | 30 | +0.0010 | 2/23/5 | 0.45 | 1 |
| astrid-multi | 60 | -0.0034 | 13/44/3 | 0.0083 | 0.042 |
| disco-astral | 60 | -0.0007 | 11/42/7 | 0.42 | 1 |

**FastMulRFS level est100: ASTRID-Pro minus method**

| astrid-pro vs | n | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| astrid-disco | 60 | +0.0010 | 13/32/15 | 0.51 | 1 |
| asteroid | 60 | +0.0012 | 17/23/20 | 0.44 | 1 |
| astral-pro3 | 60 | -0.0067 | 32/17/11 | 0.0032 | 0.016 |
| wqfm-gdl | 30 | -0.0003 | 8/11/11 | 0.85 | 1 |
| astrid-multi | 60 | -0.0046 | 21/27/12 | 0.01 | 0.042 |
| disco-astral | 60 | -0.0091 | 33/18/9 | 3.6e-05 | 0.00021 |

**FastMulRFS level est25: ASTRID-Pro minus method**

| astrid-pro vs | n | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| astrid-disco | 60 | +0.0026 | 18/14/28 | 0.37 | 1 |
| asteroid | 60 | +0.0014 | 22/13/25 | 0.87 | 1 |
| astral-pro3 | 60 | -0.0234 | 44/4/12 | 5.2e-06 | 2.6e-05 |
| wqfm-gdl | 30 | -0.0100 | 17/5/8 | 0.042 | 0.17 |
| astrid-multi | 60 | -0.0002 | 22/13/25 | 0.88 | 1 |
| disco-astral | 60 | -0.0381 | 45/4/11 | 1.3e-07 | 7.9e-07 |

**Level effects (W = first level has lower error)**

| method | contrast | n | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| astrid-pro | est100 - true | 60 | +0.0333 | 0/3/57 | 4.2e-11 |
| astrid-pro | est25 - est100 | 60 | +0.0648 | 0/0/60 | 1.6e-11 |
| astrid-pro | true tags: est - true | 60 | -0.0012 | 9/49/2 | 0.065 |
| astrid-disco | est100 - true | 60 | +0.0314 | 0/6/54 | 1.3e-10 |
| astrid-disco | est25 - est100 | 60 | +0.0632 | 0/0/60 | 1.5e-11 |
| astrid-disco | true tags: est - true | 60 | +0.0003 | 3/52/5 | 0.73 |
| asteroid | est100 - true | 60 | +0.0302 | 1/5/54 | 1.3e-10 |
| asteroid | est25 - est100 | 60 | +0.0646 | 0/1/59 | 2.3e-11 |
| astral-pro3 | est100 - true | 60 | +0.0402 | 0/0/60 | 1.4e-11 |
| astral-pro3 | est25 - est100 | 60 | +0.0814 | 1/1/58 | 2.6e-11 |
| wqfm-gdl | est100 - true | 30 | +0.0330 | 1/2/27 | 4.9e-06 |
| wqfm-gdl | est25 - est100 | 30 | +0.0759 | 0/0/30 | 1.7e-06 |
| astrid-multi | est100 - true | 60 | +0.0345 | 2/3/55 | 7.2e-11 |
| astrid-multi | est25 - est100 | 60 | +0.0603 | 0/1/59 | 2.3e-11 |
| disco-astral | est100 - true | 60 | +0.0418 | 0/0/60 | 1.4e-11 |
| disco-astral | est25 - est100 | 60 | +0.0938 | 1/1/58 | 2.5e-11 |
| disco-astral | true tags: est - true | 60 | -0.0002 | 5/50/5 | 0.88 |

**Headroom: best single method vs per-replicate oracle over the 5 main methods**

| level | n | best method | best mean FN | oracle mean FN | oracle gain | worst - best | oracle over 7 |
|---|---|---|---|---|---|---|---|
| true | 30 | wqfm-gdl | 0.007 | 0.006 | 0.001 | 0.003 | 0.005 |
| est100 | 30 | astrid-pro | 0.040 | 0.032 | 0.008 | 0.010 | 0.032 |
| est25 | 30 | astrid-disco | 0.099 | 0.087 | 0.012 | 0.029 | 0.085 |

**Q4 hybrid on FastMulRFS (diff = hybrid - ASTRAL-Pro3; seconds = mean wall time, 1 thread)**

| level | variant | n | FN variant | FN ASTRAL-Pro3 | diff | W/T/L | p | sec variant | sec ASTRAL-Pro3 |
|---|---|---|---|---|---|---|---|---|---|
| est100 | apro3-fast | 60 | 0.050 | 0.048 | +0.0017 | 0/53/7 | 0.016 | 1.093 | 6.239 |
| est100 | apro3-guide-fast | 60 | 0.048 | 0.048 | +0.0002 | 0/59/1 | 1 | 1.335 | 6.239 |
| est100 | apro3-guide | 60 | 0.048 | 0.048 | +0.0000 | 0/60/0 | 1 | 6.392 | 6.239 |
| est100 | astrid-pro | 60 | 0.042 | 0.048 | -0.0067 | 32/17/11 | 0.0032 | 0.167 | 6.239 |
| est25 | apro3-fast | 60 | 0.146 | 0.130 | +0.0165 | 5/23/32 | 3e-06 | 1.406 | 8.813 |
| est25 | apro3-guide-fast | 60 | 0.128 | 0.130 | -0.0021 | 14/37/9 | 0.31 | 1.708 | 8.813 |
| est25 | apro3-guide | 60 | 0.130 | 0.130 | +0.0002 | 0/59/1 | 1 | 8.971 | 8.813 |
| est25 | astrid-pro | 60 | 0.106 | 0.130 | -0.0234 | 44/4/12 | 5.2e-06 | 0.158 | 8.813 |

## Q1. Error decomposition: DISCO data (100 taxa, 1000 genes)

Estimated-tree (est100) results not rerun here are taken from the previous pilot (same software versions, same inputs; ASTRAL-Pro3 deterministic seed).

| cond | level | n | astrid-pro-tt | astrid-pro | astrid-disco-tt | astrid-disco | astrid-multi | asteroid | astral-pro3 | wqfm-gdl |
|---|---|---|---|---|---|---|---|---|---|---|
| default | est100 | 1 | – | 0.061 | – | **0.051** | **0.051** | **0.051** | 0.061 | – |
| default | true | 1 | **0.000** | **0.000** | **0.000** | **0.000** | **0.000** | **0.000** | **0.000** | 0.020 |
| gdl_1e-9_1 | est100 | 1 | – | **0.041** | – | 0.061 | 0.051 | 0.061 | 0.071 | – |
| gdl_1e-9_1 | true | 1 | **0.000** | **0.000** | **0.000** | **0.000** | **0.000** | **0.000** | **0.000** | – |
| gdl_5e-10_0 | est100 | 1 | – | **0.061** | – | **0.061** | **0.061** | – | – | – |
| gdl_5e-10_0 | true | 1 | **0.000** | **0.000** | **0.000** | **0.000** | **0.000** | 0.010 | – | – |

**DISCO level true**

| astrid-pro vs | n | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| astrid-disco | 3 | +0.0000 | 0/3/0 | 1 | 1 |
| asteroid | 3 | -0.0034 | 1/2/0 | 1 | 1 |
| astral-pro3 | 2 | +0.0000 | 0/2/0 | 1 | 1 |
| wqfm-gdl | 1 | -0.0204 | 1/0/0 | 1 | 1 |
| astrid-multi | 3 | +0.0000 | 0/3/0 | 1 | 1 |

Oracle over 5 methods (true, n=1): best single astrid-pro 0.0000, oracle 0.0000.

**DISCO level est100**

| astrid-pro vs | n | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| astrid-disco | 3 | -0.0034 | 1/1/1 | 1 | 1 |
| asteroid | 2 | -0.0051 | 1/0/1 | 1 | 1 |
| astral-pro3 | 2 | -0.0153 | 1/1/0 | 1 | 1 |
| astrid-multi | 3 | -0.0000 | 1/1/1 | 1 | 1 |

## Q1. DISCO: true-tree (L1) → estimated-tree (L2, 100 bp) FN rate per condition, fast methods, reps 01-07 (n in parentheses)

| cond | astrid-pro | astrid-disco | astrid-multi | asteroid |
|---|---|---|---|---|
| default | 0.003 → 0.042 (7) | 0.006 → 0.036 (7) | 0.007 → 0.045 (7) | 0.009 → 0.042 (7) |
| gdl_1e-10_0 | 0.004 → 0.041 (7) | 0.003 → 0.042 (7) | 0.004 → 0.042 (7) | 0.004 → 0.038 (7) |
| gdl_1e-10_05 | 0.006 → 0.044 (7) | 0.003 → 0.044 (7) | 0.004 → 0.044 (7) | 0.004 → 0.041 (7) |
| gdl_1e-10_1 | 0.006 → 0.028 (7) | 0.007 → 0.028 (7) | 0.006 → 0.025 (7) | 0.006 → 0.026 (7) |
| gdl_1e-9_0 | – | – | – | – |
| gdl_1e-9_05 | 0.001 → 0.031 (7) | 0.002 → 0.029 (6) | 0.009 → 0.039 (7) | 0.003 → 0.037 (3) |
| gdl_1e-9_1 | 0.007 → 0.034 (7) | 0.004 → 0.031 (7) | 0.013 → 0.045 (7) | 0.004 → 0.047 (7) |
| gdl_5e-10_0 | 0.001 → 0.022 (7) | 0.001 → 0.026 (7) | 0.007 → 0.028 (7) | 0.004 → 0.024 (5) |
| gdl_5e-10_05 | 0.006 → 0.054 (7) | 0.007 → 0.054 (7) | 0.009 → 0.060 (7) | 0.007 → 0.063 (7) |
| ils_1e4 | 0.000 → 0.037 (6) | 0.000 → 0.032 (6) | 0.000 → 0.032 (6) | 0.000 → 0.034 (6) |
| ils_2e8 | 0.026 → 0.053 (6) | 0.022 → 0.048 (6) | 0.036 → 0.054 (6) | 0.029 → 0.048 (6) |
| missing_1000 | – | – | – | – |

**DISCO true trees, all conditions**

| astrid-pro vs | n | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| astrid-disco | 70 | +0.0006 | 7/53/10 | 0.57 | 0.57 |
| asteroid | 66 | -0.0011 | 10/51/5 | 0.18 | 0.36 |
| astrid-multi | 70 | -0.0034 | 17/50/3 | 0.0041 | 0.012 |

## Tag accuracy on true gene trees (tagacc_fmrfs.jsonl)

| cond | reps | paralog pair frac | A-Pro3 false-orth | A-Pro3 false-para | DISCO/MinDup false-orth | DISCO/MinDup false-para |
|---|---|---|---|---|---|---|
| dl1e-10_ps1e7 | 10 | 0.056 | 0.011 | 0.014 | 0.011 | 0.014 |
| dl1e-10_ps5e7 | 10 | 0.044 | 0.043 | 0.043 | 0.041 | 0.042 |
| dl2e-10_ps1e7 | 10 | 0.105 | 0.020 | 0.027 | 0.012 | 0.025 |
| dl2e-10_ps5e7 | 10 | 0.093 | 0.042 | 0.077 | 0.040 | 0.075 |
| dl5e-10_ps1e7 | 10 | 0.218 | 0.043 | 0.083 | 0.034 | 0.069 |
| dl5e-10_ps5e7 | 10 | 0.199 | 0.043 | 0.156 | 0.042 | 0.157 |

## Tag accuracy on true gene trees (tagacc_disco.jsonl)

| cond | reps | paralog pair frac | A-Pro3 false-orth | A-Pro3 false-para | DISCO/MinDup false-orth | DISCO/MinDup false-para |
|---|---|---|---|---|---|---|
| default | 2 | 0.697 | 0.055 | 0.139 | 0.039 | 0.125 |
| gdl_1e-10_0 | 2 | 0.405 | 0.003 | 0.036 | 0.003 | 0.034 |
| gdl_1e-9_1 | 2 | 0.812 | 0.060 | 0.247 | 0.056 | 0.221 |
| gdl_5e-10_0 | 2 | 0.709 | 0.002 | 0.088 | 0.002 | 0.084 |
| ils_2e8 | 2 | 0.722 | 0.038 | 0.337 | 0.035 | 0.318 |

## Q2. DISCO regimes (estimated gene trees, 100 bp, 1000 genes)

| cond | n | astrid-pro | astrid-multi | astrid-disco | asteroid | astral-pro3 | wqfm-gdl | disco-astral |
|---|---|---|---|---|---|---|---|---|
| default | 7 | 0.042 | 0.045 | **0.036** | 0.042 | 0.041 (5) | – | 0.071 (1) |
| gdl_1e-10_0 | 7 | 0.041 | 0.042 | 0.042 | **0.038** | 0.053 (5) | – | – |
| gdl_1e-10_05 | 7 | 0.044 | 0.044 | 0.044 | **0.041** | 0.046 (2) | – | – |
| gdl_1e-10_1 | 7 | 0.028 | **0.025** | 0.028 | 0.026 | 0.051 (2) | – | 0.071 (1) |
| gdl_1e-9_0 | 6 | **0.058** | 0.092 | 0.077 | – | – | – | – |
| gdl_1e-9_05 | 7 | 0.031 | 0.039 | 0.029 (6) | 0.037 (3) | **0.000 (1)** | – | – |
| gdl_1e-9_1 | 7 | 0.034 | 0.045 | **0.031** | 0.047 | 0.041 (5) | – | – |
| gdl_5e-10_0 | 7 | 0.022 | 0.028 | 0.026 | 0.024 (5) | **0.020 (1)** | – | – |
| gdl_5e-10_05 | 7 | **0.054** | 0.060 | **0.054** | 0.063 | 0.056 (2) | – | – |
| ils_1e4 | 6 | 0.037 | 0.032 | 0.032 | 0.034 | **0.022 (5)** | – | – |
| ils_2e8 | 6 | 0.053 | 0.054 | 0.048 | 0.048 | **0.047 (5)** | – | – |
| missing_1000 | 6 | 0.039 | **0.034** | 0.043 | **0.034** | 0.051 (5) | – | – |

**DISCO est100, high GDL (dup >= 5e-10, excl. default)**

| astrid-pro vs | n | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| astrid-disco | 33 | -0.0034 | 11/15/7 | 0.69 | 0.88 |
| asteroid | 22 | -0.0111 | 15/5/2 | 0.0059 | 0.018 |
| astral-pro3 | 9 | -0.0057 | 4/4/1 | 0.44 | 0.88 |
| astrid-multi | 34 | -0.0126 | 23/7/4 | 0.0018 | 0.0073 |

**DISCO est100, low GDL (dup 1e-10)**

| astrid-pro vs | n | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| astrid-disco | 21 | -0.0005 | 4/13/4 | 0.88 | 1 |
| asteroid | 21 | +0.0024 | 0/16/5 | 0.062 | 0.19 |
| astral-pro3 | 9 | -0.0091 | 6/2/1 | 0.031 | 0.12 |
| astrid-multi | 21 | +0.0005 | 3/14/4 | 0.84 | 1 |

**DISCO est100, ILS (ils_1e4, ils_2e8)**

| astrid-pro vs | n | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| astrid-disco | 12 | +0.0051 | 0/7/5 | 0.062 | 0.19 |
| asteroid | 12 | +0.0043 | 2/5/5 | 0.28 | 0.56 |
| astral-pro3 | 10 | +0.0071 | 1/3/6 | 0.047 | 0.19 |
| astrid-multi | 12 | +0.0017 | 3/5/4 | 0.3 | 0.56 |

**DISCO est100, all 12 conditions**

| astrid-pro vs | n | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| astrid-disco | 79 | -0.0005 | 18/40/21 | 0.83 | 0.83 |
| asteroid | 68 | -0.0017 | 19/32/17 | 0.27 | 0.53 |
| astral-pro3 | 38 | -0.0038 | 15/15/8 | 0.17 | 0.5 |
| astrid-multi | 80 | -0.0048 | 32/32/16 | 0.016 | 0.065 |

**ASTRID-Pro minus method, per condition**

| cond | vs | n | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| default | astrid-disco | 7 | +0.0058 | 0/3/4 | 0.12 |
| default | asteroid | 7 | -0.0000 | 2/3/2 | 0.62 |
| default | astral-pro3 | 5 | -0.0061 | 2/3/0 | 0.5 |
| default | astrid-multi | 7 | -0.0029 | 3/2/2 | 0.38 |
| gdl_1e-10_0 | astrid-disco | 7 | -0.0015 | 1/6/0 | 1 |
| gdl_1e-10_0 | asteroid | 7 | +0.0029 | 0/5/2 | 0.5 |
| gdl_1e-10_0 | astral-pro3 | 5 | -0.0102 | 3/1/1 | 0.25 |
| gdl_1e-10_0 | astrid-multi | 7 | -0.0015 | 2/4/1 | 1 |
| gdl_1e-10_05 | astrid-disco | 7 | -0.0000 | 2/3/2 | 1 |
| gdl_1e-10_05 | asteroid | 7 | +0.0029 | 0/5/2 | 0.5 |
| gdl_1e-10_05 | astral-pro3 | 2 | -0.0102 | 2/0/0 | 0.5 |
| gdl_1e-10_05 | astrid-multi | 7 | -0.0000 | 1/5/1 | 1 |
| gdl_1e-10_1 | astrid-disco | 7 | +0.0000 | 1/4/2 | 1 |
| gdl_1e-10_1 | asteroid | 7 | +0.0015 | 0/6/1 | 1 |
| gdl_1e-10_1 | astral-pro3 | 2 | -0.0051 | 1/1/0 | 1 |
| gdl_1e-10_1 | astrid-multi | 7 | +0.0029 | 0/5/2 | 0.5 |
| gdl_1e-9_0 | astrid-disco | 6 | -0.0187 | 4/2/0 | 0.12 |
| gdl_1e-9_0 | astrid-multi | 6 | -0.0340 | 5/1/0 | 0.062 |
| gdl_1e-9_05 | astrid-disco | 6 | +0.0017 | 1/3/2 | 0.5 |
| gdl_1e-9_05 | asteroid | 3 | -0.0204 | 2/1/0 | 0.5 |
| gdl_1e-9_05 | astral-pro3 | 1 | +0.0000 | 0/1/0 | 1 |
| gdl_1e-9_05 | astrid-multi | 7 | -0.0087 | 6/0/1 | 0.17 |
| gdl_1e-9_1 | astrid-disco | 7 | +0.0029 | 2/1/4 | 0.75 |
| gdl_1e-9_1 | asteroid | 7 | -0.0131 | 5/0/2 | 0.22 |
| gdl_1e-9_1 | astral-pro3 | 5 | -0.0082 | 3/1/1 | 0.62 |
| gdl_1e-9_1 | astrid-multi | 7 | -0.0117 | 6/0/1 | 0.22 |
| gdl_5e-10_0 | astrid-disco | 7 | -0.0044 | 3/4/0 | 0.25 |
| gdl_5e-10_0 | asteroid | 5 | -0.0061 | 3/2/0 | 0.25 |
| gdl_5e-10_0 | astral-pro3 | 1 | +0.0000 | 0/1/0 | 1 |
| gdl_5e-10_0 | astrid-multi | 7 | -0.0058 | 3/3/1 | 0.62 |
| gdl_5e-10_05 | astrid-disco | 7 | +0.0000 | 1/5/1 | 1 |
| gdl_5e-10_05 | asteroid | 7 | -0.0087 | 5/2/0 | 0.062 |
| gdl_5e-10_05 | astral-pro3 | 2 | -0.0051 | 1/1/0 | 1 |
| gdl_5e-10_05 | astrid-multi | 7 | -0.0058 | 3/3/1 | 0.38 |
| ils_1e4 | astrid-disco | 6 | +0.0051 | 0/4/2 | 0.5 |
| ils_1e4 | asteroid | 6 | +0.0034 | 1/3/2 | 0.75 |
| ils_1e4 | astral-pro3 | 5 | +0.0082 | 0/2/3 | 0.25 |
| ils_1e4 | astrid-multi | 6 | +0.0051 | 0/3/3 | 0.25 |
| ils_2e8 | astrid-disco | 6 | +0.0051 | 0/3/3 | 0.25 |
| ils_2e8 | asteroid | 6 | +0.0051 | 1/2/3 | 0.5 |
| ils_2e8 | astral-pro3 | 5 | +0.0061 | 1/1/3 | 0.38 |
| ils_2e8 | astrid-multi | 6 | -0.0017 | 3/2/1 | 0.88 |
| missing_1000 | astrid-disco | 6 | -0.0034 | 3/2/1 | 0.88 |
| missing_1000 | asteroid | 6 | +0.0051 | 0/3/3 | 0.25 |
| missing_1000 | astral-pro3 | 5 | -0.0102 | 2/3/0 | 0.5 |
| missing_1000 | astrid-multi | 6 | +0.0051 | 0/4/2 | 0.5 |

## Q2. Few vs many genes (DISCO gtrees_10000_l1, 100 taxa, 100 bp)

| level | ngen | n | astrid-pro | astrid-multi | astrid-disco | asteroid | astral-pro3 | wqfm-gdl |
|---|---|---|---|---|---|---|---|---|
| est100 | 100 | 4 | 0.056 | 0.079 | 0.059 | 0.071 | 0.056 | **0.038** |
| est100 | 1000 | 4 | 0.026 | 0.028 | 0.026 | **0.020** | – | – |
| est100 | 10000 | 4 | **0.026** | 0.028 | **0.026** | – | – | – |
| est100 | 50 | 4 | 0.079 | 0.128 | **0.064** | 0.092 | 0.082 | 0.066 |
| true | 100 | 4 | 0.018 | – | **0.018** | – | – | – |
| true | 1000 | 4 | **0.000** | – | **0.000** | – | – | – |
| true | 10000 | 4 | **0.000** | – | 0.003 | – | – | – |

| ngen | ASTRID-Pro vs | n | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | astrid-disco | 4 | -0.0025 | 2/1/1 | 1 |
| 100 | asteroid | 4 | -0.0153 | 2/2/0 | 0.5 |
| 100 | astral-pro3 | 4 | +0.0000 | 1/1/2 | 1 |
| 100 | wqfm-gdl | 4 | +0.0179 | 0/0/4 | 0.12 |
| 1000 | astrid-disco | 4 | +0.0000 | 1/2/1 | 1 |
| 1000 | asteroid | 4 | +0.0051 | 0/2/2 | 0.5 |
| 10000 | astrid-disco | 4 | +0.0000 | 0/4/0 | 1 |
| 50 | astrid-disco | 4 | +0.0153 | 1/0/3 | 0.62 |
| 50 | asteroid | 4 | -0.0128 | 3/0/1 | 0.75 |
| 50 | astral-pro3 | 4 | -0.0025 | 2/0/2 | 0.88 |
| 50 | wqfm-gdl | 4 | +0.0128 | 0/2/2 | 0.5 |

## Q2. 1000 species (DISCO species_1000)

| level | ngen | n | astrid-pro-tt | astrid-pro | astrid-multi | astrid-disco | asteroid | astral-pro3 | wqfm-gdl |
|---|---|---|---|---|---|---|---|---|---|
| est100 | 1000 | 2 | – | 0.101 | 0.107 | **0.100** | – | – | – |
| true | 1000 | 1 | 0.007 | 0.006 | 0.013 | **0.003** | – | – | – |

**species_1000, level est100**

| astrid-pro vs | n | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| astrid-disco | 2 | +0.0010 | 1/0/1 | 1 | 1 |
| astrid-multi | 2 | -0.0065 | 2/0/0 | 0.5 | 1 |

**species_1000, level true**

| astrid-pro vs | n | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| astrid-disco | 1 | +0.0030 | 0/0/1 | 1 | 1 |
| astrid-multi | 1 | -0.0070 | 1/0/0 | 1 | 1 |

## Q4. Hybrid on DISCO (1000 genes, 100 bp; ASTRAL-Pro3 from previous pilot, mean 578 s)

| variant | n | FN | diff vs ASTRAL-Pro3 | W/T/L | p | sec |
|---|---|---|---|---|---|---|
| apro3-fast | 10 | 0.055 | +0.0000 | 0/10/0 | 1 | 53.764 |
| apro3-guide-fast | 10 | 0.055 | +0.0000 | 0/10/0 | 1 | 60.237 |
| astrid-pro | 10 | 0.049 | -0.0061 | 5/3/2 | 0.11 | 1.084 |

## Q3. Scaling (idle machine, one job at a time)

**wall time (s)**

|                       | asteroid@1t   | astral-pro3@4t   |   astrid-disco@1t |   astrid-pro@1t |   astrid-pro@4t | wqfm-gdl@4t   |
|:----------------------|:--------------|:-----------------|------------------:|----------------:|----------------:|:--------------|
| ('genes', '100')      | 0.8           | 4.2              |               1   |             0.2 |             0.1 | 62.7          |
| ('genes', '1000')     | 8.0           | 41.5             |               8   |             0.7 |             0.3 | –             |
| ('genes', '10000')    | –             | –                |              77.8 |             6.8 |             2.1 | –             |
| ('taxa_100', '1000')  | 8.7           | 34.1             |               7.6 |             0.7 |             0.3 | –             |
| ('taxa_1000', '1000') | –             | –                |             536.4 |            14.3 |             7.6 | –             |
| ('taxa_200', '1000')  | 46.5          | 164.7            |              17.6 |             2   |             1.1 | –             |
| ('taxa_500', '1000')  | –             | –                |              57.8 |             4.4 |             1.7 | –             |

**peak RSS (MB)**

|                       | asteroid@1t   | astral-pro3@4t   |   astrid-disco@1t |   astrid-pro@1t |   astrid-pro@4t | wqfm-gdl@4t   |
|:----------------------|:--------------|:-----------------|------------------:|----------------:|----------------:|:--------------|
| ('genes', '100')      | 140.6         | 13.9             |              45.5 |            13.6 |            13.6 | 2055.8        |
| ('genes', '1000')     | 1392.2        | 51.8             |              51.1 |            13.7 |            13.6 | –             |
| ('genes', '10000')    | –             | –                |             118.3 |            13.8 |            13.8 | –             |
| ('taxa_100', '1000')  | 1486.9        | 57.5             |              49.2 |            13.4 |            13.6 | –             |
| ('taxa_1000', '1000') | –             | –                |             604.5 |           113   |           211.8 | –             |
| ('taxa_200', '1000')  | 5930.2        | 107.7            |              90.7 |            13.7 |            21.5 | –             |
| ('taxa_500', '1000')  | –             | –                |             215.3 |            30.6 |            66.8 | –             |

**FN rate**

|                       | asteroid@1t   | astral-pro3@4t   |   astrid-disco@1t |   astrid-pro@1t |   astrid-pro@4t | wqfm-gdl@4t   |
|:----------------------|:--------------|:-----------------|------------------:|----------------:|----------------:|:--------------|
| ('genes', '100')      | 0.061         | 0.061            |             0.082 |           0.061 |           0.061 | 0.031         |
| ('genes', '1000')     | 0.061         | 0.061            |             0.061 |           0.061 |           0.061 | –             |
| ('genes', '10000')    | –             | –                |             0.061 |           0.061 |           0.061 | –             |
| ('taxa_100', '1000')  | 0.082         | 0.062            |             0.072 |           0.072 |           0.072 | –             |
| ('taxa_1000', '1000') | –             | –                |             0.072 |           0.07  |           0.07  | –             |
| ('taxa_200', '1000')  | 0.086         | 0.066            |             0.061 |           0.071 |           0.071 | –             |
| ('taxa_500', '1000')  | –             | –                |             0.058 |           0.062 |           0.062 | –             |

## Q4b. Search or objective? ASTRAL-Pro3 quartet score of the true tree vs ASTRAL-Pro3's tree (FastMulRFS, 100 genes, reps 01-05)

| level | n | true scores higher (search failure) | equal | true scores lower (objective failure) | score gap ASTRAL - true (%) | ASTRID-Pro tree scores lower than ASTRAL's | ... and has fewer FN | mean FN ASTRAL-Pro3 | mean FN ASTRID-Pro |
|---|---|---|---|---|---|---|---|---|---|
| est100 | 30 | 0 | 0 | 30 | 0.140 | 30 | 19 | 4.800 | 3.833 |
| est25 | 30 | 0 | 0 | 30 | 0.663 | 30 | 24 | 12.367 | 10.267 |

## Failed runs

| file | method | cond | error | n |
|---|---|---|---|---|
| disco_q1.jsonl | astral-pro3 | gdl_5e-10_0 | TimeoutExpired | 1 |
| disco_q2.jsonl | astral-pro3 | gdl_1e-9_0 | TimeoutExpired | 1 |
| scaling.jsonl | asteroid | taxa_100 | CalledProcessError | 1 |
| scaling.jsonl | asteroid | taxa_1000 | CalledProcessError | 1 |
| scaling.jsonl | asteroid | taxa_500 | CalledProcessError | 1 |
| scaling.jsonl | astral-pro3 | taxa_500 | TimeoutExpired | 1 |

