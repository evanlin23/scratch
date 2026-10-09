Replicates: 369 ok, 0 errors.


### short U[0.005,0.05] (Kim et al.), n=50

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 16 | 0.00 | 0.157 | 0.149 | 0.144 | 0.133 | 0.573 [0.046] | 0.168 | 0.166 | 0.153 | 0.156 | 0.160 | 0.161 | 5.4 | 1.06 |
| 200 | 16 | 0.00 | 0.068 | 0.064 | 0.056 | 0.045 | 0.379 [0.008] | 0.065 | 0.069 | 0.059 | 0.059 | 0.064 | 0.062 | 1.9 | 0.25 |
| 500 | 16 | 0.00 | 0.009 | 0.009 | 0.009 | 0.008 | 0.199 [0.003] | 0.011 | 0.011 | 0.011 | 0.011 | 0.009 | 0.009 | 1.5 | 0.12 |
| 1000 | 17 | 0.00 | 0.004 | 0.001 | 0.003 | 0.000 | 0.128 [0.007] | 0.010 | 0.010 | 0.009 | 0.009 | 0.004 | 0.001 | 1.0 | 0.29 |
| 2000 | 17 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.083 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.2 | 0.00 |
| 5000 | 17 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.028 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | +0.0106 | 1/9/6 | 0.042 |
| 100 | Forest+GTM(FastME) | FastME | +0.0093 | 2/9/5 | 0.11 |
| 100 | Forest+GTM(NJ) | FastME | +0.0239 | 1/5/10 | 0.0043 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0080 | 0/12/4 | 0.063 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0066 | 4/5/7 | 0.46 |
| 100 | Forest+GTM(NJ) | FastTree | +0.0346 | 2/3/11 | 0.0069 |
| 200 | Forest+GTM(NJ) | NJ | -0.0027 | 3/12/1 | 0.36 |
| 200 | Forest+GTM(FastME) | FastME | +0.0027 | 0/14/2 | 0.16 |
| 200 | Forest+GTM(NJ) | FastME | +0.0093 | 3/6/7 | 0.1 |
| 200 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0013 | 0/15/1 | 0.32 |
| 200 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0027 | 4/5/7 | 0.82 |
| 200 | Forest+GTM(NJ) | FastTree | +0.0199 | 3/4/9 | 0.04 |
| 500 | Forest+GTM(NJ) | NJ | +0.0013 | 0/15/1 | 0.32 |
| 500 | Forest+GTM(FastME) | FastME | +0.0013 | 0/15/1 | 0.32 |
| 500 | Forest+GTM(NJ) | FastME | +0.0013 | 0/15/1 | 0.32 |
| 500 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0013 | 0/15/1 | 0.32 |
| 500 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0013 | 1/13/2 | 0.56 |
| 500 | Forest+GTM(NJ) | FastTree | +0.0027 | 0/14/2 | 0.16 |
| 1000 | Forest+GTM(NJ) | NJ | +0.0063 | 0/14/3 | 0.1 |
| 1000 | Forest+GTM(FastME) | FastME | +0.0063 | 0/14/3 | 0.1 |
| 1000 | Forest+GTM(NJ) | FastME | +0.0075 | 0/13/4 | 0.063 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0063 | 0/14/3 | 0.1 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0088 | 0/12/5 | 0.038 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.0100 | 0/11/6 | 0.023 |
| 2000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/17/0 | 1 |
| 2000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/17/0 | 1 |
| 2000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/17/0 | 1 |
| 2000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/17/0 | 1 |
| 2000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/17/0 | 1 |
| 2000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/17/0 | 1 |
| 5000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/17/0 | 1 |
| 5000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/17/0 | 1 |
| 5000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/17/0 | 1 |
| 5000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/17/0 | 1 |
| 5000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/17/0 | 1 |
| 5000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/17/0 | 1 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.44 | 0.00 | 0.00 | 21.1 / 47 |
| 200 | 0.00 | 0.75 | 0.00 | 0.00 | 29.4 / 47 |
| 500 | 0.56 | 0.94 | 0.56 | 0.00 | 37.8 / 47 |
| 1000 | 0.82 | 0.82 | 0.82 | 0.00 | 41.3 / 47 |
| 2000 | 1.00 | 1.00 | 1.00 | 0.00 | 43.1 / 47 |
| 5000 | 1.00 | 1.00 | 1.00 | 0.47 | 45.7 / 47 |

### long U[0.05,0.1] (Kim et al.), n=100

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 16 | 0.00 | 0.092 | 0.085 | 0.035 | 0.020 | 0.654 [0.006] | 0.064 | 0.085 | 0.040 | 0.037 | 0.062 | 0.049 | 18.3 | 0.25 |
| 300 | 16 | 0.00 | 0.007 | 0.006 | 0.002 | 0.000 | 0.364 [0.002] | 0.005 | 0.008 | 0.003 | 0.003 | 0.004 | 0.004 | 8.1 | 0.12 |
| 1000 | 15 | 0.00 | 0.000 | 0.001 | 0.000 | 0.000 | 0.140 [0.001] | 0.001 | 0.001 | 0.001 | 0.001 | 0.000 | 0.000 | 3.5 | 0.07 |
| 3000 | 15 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.012 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |
| 10000 | 15 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.002 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |
| 100000 | 15 | 0.00 | 0.000 | 0.000 | 0.000 | — | 0.011 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.3 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | -0.0284 | 15/1/0 | 0.00064 |
| 100 | Forest+GTM(FastME) | FastME | +0.0052 | 1/8/7 | 0.02 |
| 100 | Forest+GTM(NJ) | FastME | +0.0290 | 0/4/12 | 0.0022 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0019 | 0/14/2 | 0.18 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0148 | 4/1/11 | 0.046 |
| 100 | Forest+GTM(NJ) | FastTree | +0.0438 | 0/1/15 | 0.00064 |
| 300 | Forest+GTM(NJ) | NJ | -0.0026 | 4/11/1 | 0.16 |
| 300 | Forest+GTM(FastME) | FastME | +0.0013 | 0/14/2 | 0.16 |
| 300 | Forest+GTM(NJ) | FastME | +0.0026 | 0/12/4 | 0.046 |
| 300 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0006 | 0/15/1 | 0.32 |
| 300 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0006 | 2/11/3 | 0.65 |
| 300 | Forest+GTM(NJ) | FastTree | +0.0045 | 0/10/6 | 0.02 |
| 1000 | Forest+GTM(NJ) | NJ | +0.0007 | 0/14/1 | 0.32 |
| 1000 | Forest+GTM(FastME) | FastME | +0.0007 | 0/14/1 | 0.32 |
| 1000 | Forest+GTM(NJ) | FastME | +0.0007 | 0/14/1 | 0.32 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0007 | 0/14/1 | 0.32 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0007 | 0/14/1 | 0.32 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.0007 | 0/14/1 | 0.32 |
| 3000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/15/0 | 1 |
| 3000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/15/0 | 1 |
| 3000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/15/0 | 1 |
| 3000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/15/0 | 1 |
| 3000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/15/0 | 1 |
| 3000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/15/0 | 1 |
| 10000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/15/0 | 1 |
| 10000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/15/0 | 1 |
| 10000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/15/0 | 1 |
| 10000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/15/0 | 1 |
| 10000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/15/0 | 1 |
| 10000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/15/0 | 1 |
| 100000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/15/0 | 1 |
| 100000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/15/0 | 1 |
| 100000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/15/0 | 1 |
| 100000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/15/0 | 1 |
| 100000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/15/0 | 1 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.81 | 0.00 | 0.00 | 33.8 / 97 |
| 300 | 0.38 | 0.88 | 0.38 | 0.00 | 61.8 / 97 |
| 1000 | 1.00 | 0.93 | 1.00 | 0.07 | 83.5 / 97 |
| 3000 | 1.00 | 1.00 | 1.00 | 0.67 | 95.9 / 97 |
| 10000 | 1.00 | 1.00 | 1.00 | 0.87 | 96.8 / 97 |
| 100000 | 1.00 | 1.00 | 1.00 | 0.87 | 95.9 / 97 |

### deep U[0.1,0.4], n=100

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 15 | 0.39 | 0.423 | 0.414 | 0.318 | 0.187 | 0.908 [0.006] | 0.402 | 0.412 | 0.311 | 0.315 | 0.404 | 0.399 | 42.1 | 0.07 |
| 300 | 15 | 0.31 | 0.298 | 0.280 | 0.122 | 0.024 | 0.766 [0.009] | 0.210 | 0.262 | 0.113 | 0.122 | 0.210 | 0.207 | 22.1 | 0.27 |
| 1000 | 15 | 0.23 | 0.254 | 0.186 | 0.047 | 0.001 | 0.674 [0.007] | 0.136 | 0.195 | 0.045 | 0.048 | 0.140 | 0.100 | 16.2 | 0.27 |
| 3000 | 15 | 0.16 | 0.193 | 0.142 | 0.018 | 0.000 | 0.615 [0.004] | 0.071 | 0.138 | 0.017 | 0.019 | 0.076 | 0.043 | 12.2 | 0.13 |
| 10000 | 15 | 0.09 | 0.115 | 0.090 | 0.003 | 0.000 | 0.429 [0.002] | 0.032 | 0.067 | 0.003 | 0.003 | 0.047 | 0.014 | 7.9 | 0.13 |
| 100000 | 15 | 0.03 | 0.049 | 0.034 | 0.000 | — | 0.244 [0.004] | 0.016 | 0.021 | 0.004 | 0.004 | 0.022 | 0.001 | 4.6 | 0.27 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | -0.0213 | 13/2/0 | 0.0014 |
| 100 | Forest+GTM(FastME) | FastME | -0.0069 | 7/7/1 | 0.025 |
| 100 | Forest+GTM(NJ) | FastME | +0.0838 | 2/1/12 | 0.0057 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0021 | 3/12/0 | 0.1 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0034 | 6/1/8 | 0.71 |
| 100 | Forest+GTM(NJ) | FastTree | +0.2151 | 0/0/15 | 0.00065 |
| 300 | Forest+GTM(NJ) | NJ | -0.0887 | 15/0/0 | 0.00065 |
| 300 | Forest+GTM(FastME) | FastME | -0.0082 | 9/4/2 | 0.067 |
| 300 | Forest+GTM(NJ) | FastME | +0.0880 | 2/1/12 | 0.0023 |
| 300 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0000 | 2/11/2 | 0.85 |
| 300 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0027 | 6/3/6 | 0.75 |
| 300 | Forest+GTM(NJ) | FastTree | +0.1856 | 0/0/15 | 0.00065 |
| 1000 | Forest+GTM(NJ) | NJ | -0.1175 | 15/0/0 | 0.00065 |
| 1000 | Forest+GTM(FastME) | FastME | -0.0021 | 5/8/2 | 0.31 |
| 1000 | Forest+GTM(NJ) | FastME | +0.0887 | 1/0/14 | 0.0048 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0034 | 2/11/2 | 0.47 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0364 | 3/0/12 | 0.011 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.1354 | 0/0/15 | 0.00063 |
| 3000 | Forest+GTM(NJ) | NJ | -0.1216 | 15/0/0 | 0.00065 |
| 3000 | Forest+GTM(FastME) | FastME | -0.0007 | 3/10/2 | 0.89 |
| 3000 | Forest+GTM(NJ) | FastME | +0.0536 | 0/2/13 | 0.0014 |
| 3000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0048 | 5/9/1 | 0.21 |
| 3000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0282 | 2/1/12 | 0.0089 |
| 3000 | Forest+GTM(NJ) | FastTree | +0.0715 | 0/0/15 | 0.00064 |
| 10000 | Forest+GTM(NJ) | NJ | -0.0832 | 14/1/0 | 0.00098 |
| 10000 | Forest+GTM(FastME) | FastME | +0.0007 | 0/14/1 | 0.32 |
| 10000 | Forest+GTM(NJ) | FastME | +0.0296 | 0/3/12 | 0.0021 |
| 10000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0144 | 5/10/0 | 0.042 |
| 10000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0186 | 2/2/11 | 0.0078 |
| 10000 | Forest+GTM(NJ) | FastTree | +0.0323 | 0/2/13 | 0.0014 |
| 100000 | Forest+GTM(NJ) | NJ | -0.0330 | 12/2/1 | 0.003 |
| 100000 | Forest+GTM(FastME) | FastME | +0.0041 | 0/11/4 | 0.063 |
| 100000 | Forest+GTM(NJ) | FastME | +0.0165 | 0/9/6 | 0.027 |
| 100000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0055 | 4/9/2 | 0.6 |
| 100000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0158 | 1/8/6 | 0.034 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.93 | 0.00 | 0.00 | 9.0 / 97 |
| 300 | 0.00 | 0.73 | 0.00 | 0.00 | 23.0 / 97 |
| 1000 | 0.00 | 0.87 | 0.00 | 0.00 | 31.9 / 97 |
| 3000 | 0.00 | 0.87 | 0.00 | 0.00 | 37.5 / 97 |
| 10000 | 0.00 | 0.87 | 0.00 | 0.00 | 55.5 / 97 |
| 100000 | 0.13 | 0.73 | 0.13 | 0.00 | 73.6 / 97 |

### ultrametric h=2, lognormal(1) rates, K2P, n=100

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 15 | 0.53 | 0.909 | 0.845 | 0.864 | 0.475 | 0.908 [0.166] | 0.762 | 0.821 | 0.745 | 0.794 | 0.815 | 0.879 | 16.8 | 2.13 |
| 300 | 15 | 0.47 | 0.890 | 0.805 | 0.847 | 0.322 | 0.883 [0.118] | 0.748 | 0.790 | 0.717 | 0.759 | 0.792 | 0.854 | 18.4 | 1.67 |
| 1000 | 15 | 0.44 | 0.845 | 0.757 | 0.784 | 0.216 | 0.847 [0.030] | 0.691 | 0.744 | 0.653 | 0.714 | 0.718 | 0.784 | 20.2 | 0.40 |
| 3000 | 15 | 0.35 | 0.757 | 0.682 | 0.669 | 0.132 | 0.788 [0.021] | 0.558 | 0.649 | 0.524 | 0.605 | 0.566 | 0.691 | 21.0 | 0.33 |
| 10000 | 14 | 0.26 | 0.715 | 0.620 | 0.580 | 0.102 | 0.744 [0.011] | 0.531 | 0.613 | 0.433 | 0.518 | 0.536 | 0.649 | 22.6 | 0.21 |
| 100000 | 14 | 0.16 | 0.522 | 0.432 | 0.371 | — | 0.628 [0.000] | 0.364 | 0.444 | 0.279 | 0.326 | 0.364 | 0.404 | 18.6 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | -0.1471 | 15/0/0 | 0.00065 |
| 100 | Forest+GTM(FastME) | FastME | -0.1189 | 15/0/0 | 0.00063 |
| 100 | Forest+GTM(NJ) | FastME | -0.1017 | 13/1/1 | 0.0014 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0529 | 11/3/1 | 0.0036 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1168 | 15/0/0 | 6.1e-05 |
| 100 | Forest+GTM(NJ) | FastTree | +0.2873 | 0/0/15 | 0.00065 |
| 300 | Forest+GTM(NJ) | NJ | -0.1423 | 15/0/0 | 0.00065 |
| 300 | Forest+GTM(FastME) | FastME | -0.1299 | 15/0/0 | 0.00065 |
| 300 | Forest+GTM(NJ) | FastME | -0.0990 | 13/0/2 | 0.0029 |
| 300 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0440 | 9/6/0 | 0.0076 |
| 300 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1058 | 15/0/0 | 0.00065 |
| 300 | Forest+GTM(NJ) | FastTree | +0.4254 | 0/0/15 | 0.00065 |
| 1000 | Forest+GTM(NJ) | NJ | -0.1540 | 15/0/0 | 0.00065 |
| 1000 | Forest+GTM(FastME) | FastME | -0.1306 | 15/0/0 | 0.00065 |
| 1000 | Forest+GTM(NJ) | FastME | -0.0921 | 11/2/2 | 0.0071 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0268 | 6/8/1 | 0.028 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0921 | 11/0/4 | 0.0064 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.4756 | 0/0/15 | 0.00065 |
| 3000 | Forest+GTM(NJ) | NJ | -0.1993 | 15/0/0 | 0.00065 |
| 3000 | Forest+GTM(FastME) | FastME | -0.1450 | 15/0/0 | 0.00064 |
| 3000 | Forest+GTM(NJ) | FastME | -0.1107 | 12/1/2 | 0.0019 |
| 3000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0082 | 2/13/0 | 0.18 |
| 3000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1326 | 13/1/1 | 0.0012 |
| 3000 | Forest+GTM(NJ) | FastTree | +0.4261 | 0/0/15 | 0.00065 |
| 10000 | Forest+GTM(NJ) | NJ | -0.1841 | 14/0/0 | 0.00097 |
| 10000 | Forest+GTM(FastME) | FastME | -0.1465 | 14/0/0 | 0.00012 |
| 10000 | Forest+GTM(NJ) | FastME | -0.0486 | 9/0/5 | 0.29 |
| 10000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0052 | 2/12/0 | 0.18 |
| 10000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1186 | 12/1/1 | 0.0019 |
| 10000 | Forest+GTM(NJ) | FastTree | +0.4286 | 0/0/14 | 0.00097 |
| 100000 | Forest+GTM(NJ) | NJ | -0.1583 | 14/0/0 | 0.00097 |
| 100000 | Forest+GTM(FastME) | FastME | -0.0920 | 14/0/0 | 0.00098 |
| 100000 | Forest+GTM(NJ) | FastME | -0.0074 | 7/1/6 | 0.62 |
| 100000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/14/0 | 1 |
| 100000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0398 | 12/0/2 | 0.024 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.33 | 0.00 | 0.00 | 11.1 / 97 |
| 300 | 0.00 | 0.47 | 0.00 | 0.00 | 13.0 / 97 |
| 1000 | 0.00 | 0.67 | 0.00 | 0.00 | 15.3 / 97 |
| 3000 | 0.00 | 0.73 | 0.00 | 0.00 | 20.9 / 97 |
| 10000 | 0.00 | 0.93 | 0.00 | 0.00 | 25.0 / 97 |
| 100000 | 0.00 | 1.00 | 0.00 | 0.00 | 36.1 / 97 |

### Runtime (mean seconds per replicate, single core)

| regime | n | k | NJ | FastME | FastTree | Forest (grid search) | Forest+GTM(NJ) total | Centroid dec.+GTM(NJ) |
|---|---|---|---|---|---|---|---|---|
| U:0.005:0.05 | 50 | 100 | 0.01 | 0.01 | 0.2 | 1.9 | 2.2 | 0.34 |
| U:0.005:0.05 | 50 | 200 | 0.01 | 0.01 | 0.4 | 1.5 | 1.8 | 0.31 |
| U:0.005:0.05 | 50 | 500 | 0.01 | 0.02 | 0.8 | 0.8 | 1.1 | 0.29 |
| U:0.005:0.05 | 50 | 1000 | 0.02 | 0.02 | 1.7 | 0.7 | 1.0 | 0.30 |
| U:0.005:0.05 | 50 | 2000 | 0.02 | 0.02 | 3.5 | 0.6 | 1.0 | 0.30 |
| U:0.005:0.05 | 50 | 5000 | 0.02 | 0.02 | 8.6 | 0.6 | 0.9 | 0.30 |
| U:0.05:0.1 | 100 | 100 | 0.08 | 0.06 | 0.4 | 11.6 | 12.0 | 0.40 |
| U:0.05:0.1 | 100 | 300 | 0.10 | 0.07 | 1.1 | 6.1 | 6.5 | 0.43 |
| U:0.05:0.1 | 100 | 1000 | 0.11 | 0.07 | 3.8 | 3.3 | 3.7 | 0.42 |
| U:0.05:0.1 | 100 | 3000 | 0.15 | 0.11 | 12.6 | 3.3 | 3.9 | 0.50 |
| U:0.05:0.1 | 100 | 10000 | 0.13 | 0.10 | 44.4 | 3.1 | 3.6 | 0.44 |
| U:0.05:0.1 | 100 | 100000 | 0.12 | 0.08 |  | 2.9 | 3.4 | 0.43 |
| U:0.1:0.4 | 100 | 100 | 0.06 | 0.12 | 0.6 | 8.7 | 9.0 | 0.36 |
| U:0.1:0.4 | 100 | 300 | 0.10 | 0.09 | 1.3 | 10.8 | 11.2 | 0.37 |
| U:0.1:0.4 | 100 | 1000 | 0.09 | 0.09 | 4.0 | 10.9 | 11.4 | 0.39 |
| U:0.1:0.4 | 100 | 3000 | 0.10 | 0.08 | 13.4 | 13.1 | 13.5 | 0.45 |
| U:0.1:0.4 | 100 | 10000 | 0.07 | 0.05 | 45.0 | 9.6 | 10.0 | 0.41 |
| U:0.1:0.4 | 100 | 100000 | 0.09 | 0.07 |  | 4.9 | 5.3 | 0.37 |
| UH:2.0:1.0 | 100 | 100 | 0.04 | 0.14 | 1.0 | 14.0 | 14.3 | 0.34 |
| UH:2.0:1.0 | 100 | 300 | 0.05 | 0.13 | 2.4 | 13.4 | 13.7 | 0.32 |
| UH:2.0:1.0 | 100 | 1000 | 0.05 | 0.13 | 6.6 | 19.9 | 20.3 | 0.36 |
| UH:2.0:1.0 | 100 | 3000 | 0.05 | 0.13 | 19.0 | 32.6 | 33.0 | 0.33 |
| UH:2.0:1.0 | 100 | 10000 | 0.04 | 0.11 | 66.7 | 56.5 | 56.9 | 0.36 |
| UH:2.0:1.0 | 100 | 100000 | 0.05 | 0.09 |  | 68.4 | 68.8 | 0.34 |

### Follow-up: saturation handling and FastME-guided controls (6 replicates, 0 errors)

| regime | k | reps | NJ_cap2 | NJ_cap1.2 | NJ_cap5 | NJ_pcap | FastME_cap2 | FastME_cap1.2 | FastME_cap5 | FastME_pcap | FGTM_FastME | CompGTM_FastME | DecGTM_FastME_25 | DecGTM_FastME_50 | FGTM_FastME_pcap |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| U:0.1:0.4 | 300 | 2 | 0.247 | 0.082 | 0.825 | 0.222 | 0.046 | 0.021 | 0.263 | 0.082 | 0.057 | 0.046 | 0.046 | 0.046 | 0.077 |
| U:0.1:0.4 | 3000 | 2 | 0.196 | 0.062 | 0.562 | 0.046 | 0.010 | 0.000 | 0.046 | 0.000 | 0.010 | 0.010 | 0.010 | 0.015 | 0.000 |
| U:0.1:0.4 | 100000 | 1 | 0.072 | 0.000 | 0.216 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| UH:2.0:1.0 | 300 | 1 | 0.876 | 0.866 | 0.897 | 0.340 | 0.742 | 0.753 | 0.773 | 0.309 | 0.608 | 0.608 | 0.742 | 0.753 | 0.289 |

| regime | k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|
| U:0.1:0.4 | 300 | FGTM_FastME_pcap | FastME_pcap | -0.0052 | 1/1/0 | 1 |
| U:0.1:0.4 | 300 | FGTM_FastME | FastME_cap2 | +0.0103 | 0/1/1 | 1 |
| U:0.1:0.4 | 300 | FGTM_FastME | CompGTM_FastME | +0.0103 | 0/0/2 | 0.5 |
| U:0.1:0.4 | 300 | FGTM_FastME | DecGTM_FastME_25 | +0.0103 | 1/0/1 | 1 |
| U:0.1:0.4 | 300 | FastME_pcap | FastME_cap2 | +0.0361 | 0/1/1 | 1 |
| U:0.1:0.4 | 3000 | FGTM_FastME_pcap | FastME_pcap | +0.0000 | 0/2/0 | 1 |
| U:0.1:0.4 | 3000 | FGTM_FastME | FastME_cap2 | +0.0000 | 0/2/0 | 1 |
| U:0.1:0.4 | 3000 | FGTM_FastME | CompGTM_FastME | +0.0000 | 0/2/0 | 1 |
| U:0.1:0.4 | 3000 | FGTM_FastME | DecGTM_FastME_25 | +0.0000 | 0/2/0 | 1 |
| U:0.1:0.4 | 3000 | FastME_pcap | FastME_cap2 | -0.0103 | 2/0/0 | 0.5 |
| U:0.1:0.4 | 100000 | FGTM_FastME_pcap | FastME_pcap | +0.0000 | 0/1/0 | 1 |
| U:0.1:0.4 | 100000 | FGTM_FastME | FastME_cap2 | +0.0000 | 0/1/0 | 1 |
| U:0.1:0.4 | 100000 | FGTM_FastME | CompGTM_FastME | +0.0000 | 0/1/0 | 1 |
| U:0.1:0.4 | 100000 | FGTM_FastME | DecGTM_FastME_25 | +0.0000 | 0/1/0 | 1 |
| U:0.1:0.4 | 100000 | FastME_pcap | FastME_cap2 | +0.0000 | 0/1/0 | 1 |
| UH:2.0:1.0 | 300 | FGTM_FastME_pcap | FastME_pcap | -0.0206 | 1/0/0 | 1 |
| UH:2.0:1.0 | 300 | FGTM_FastME | FastME_cap2 | -0.1340 | 1/0/0 | 1 |
| UH:2.0:1.0 | 300 | FGTM_FastME | CompGTM_FastME | +0.0000 | 0/1/0 | 1 |
| UH:2.0:1.0 | 300 | FGTM_FastME | DecGTM_FastME_25 | -0.1340 | 1/0/0 | 1 |
| UH:2.0:1.0 | 300 | FastME_pcap | FastME_cap2 | -0.4330 | 1/0/0 | 1 |
