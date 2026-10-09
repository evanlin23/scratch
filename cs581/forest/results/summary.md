Replicates: 242 ok, 0 errors.


### short U[0.005,0.05] (Kim et al.), n=50

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 11 | 0.00 | 0.164 | 0.159 | 0.157 | 0.137 | 0.613 [0.051] | 0.176 | 0.172 | 0.162 | 0.166 | 0.170 | 0.174 | 6.7 | 1.18 |
| 200 | 11 | 0.00 | 0.072 | 0.066 | 0.060 | 0.048 | 0.350 [0.012] | 0.068 | 0.074 | 0.064 | 0.064 | 0.066 | 0.070 | 1.5 | 0.36 |
| 500 | 11 | 0.00 | 0.008 | 0.010 | 0.008 | 0.006 | 0.209 [0.005] | 0.010 | 0.010 | 0.010 | 0.010 | 0.008 | 0.008 | 1.7 | 0.18 |
| 1000 | 11 | 0.00 | 0.004 | 0.002 | 0.004 | 0.000 | 0.137 [0.007] | 0.010 | 0.010 | 0.010 | 0.010 | 0.004 | 0.002 | 1.0 | 0.27 |
| 2000 | 11 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.081 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |
| 5000 | 10 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.017 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | +0.0116 | 0/7/4 | 0.12 |
| 100 | Forest+GTM(FastME) | FastME | +0.0058 | 2/6/3 | 0.38 |
| 100 | Forest+GTM(NJ) | FastME | +0.0193 | 1/4/6 | 0.031 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0058 | 0/9/2 | 0.5 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0019 | 3/4/4 | 0.98 |
| 100 | Forest+GTM(NJ) | FastTree | +0.0387 | 1/1/9 | 0.0078 |
| 200 | Forest+GTM(NJ) | NJ | -0.0039 | 3/7/1 | 0.5 |
| 200 | Forest+GTM(FastME) | FastME | +0.0039 | 0/9/2 | 0.5 |
| 200 | Forest+GTM(NJ) | FastME | +0.0077 | 2/4/5 | 0.12 |
| 200 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0019 | 0/10/1 | 1 |
| 200 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0019 | 3/4/4 | 0.98 |
| 200 | Forest+GTM(NJ) | FastTree | +0.0193 | 2/3/6 | 0.1 |
| 500 | Forest+GTM(NJ) | NJ | +0.0019 | 0/10/1 | 1 |
| 500 | Forest+GTM(FastME) | FastME | +0.0019 | 0/10/1 | 1 |
| 500 | Forest+GTM(NJ) | FastME | +0.0019 | 0/10/1 | 1 |
| 500 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0019 | 0/10/1 | 1 |
| 500 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0019 | 1/8/2 | 1 |
| 500 | Forest+GTM(NJ) | FastTree | +0.0039 | 0/9/2 | 0.5 |
| 1000 | Forest+GTM(NJ) | NJ | +0.0058 | 0/9/2 | 0.5 |
| 1000 | Forest+GTM(FastME) | FastME | +0.0058 | 0/9/2 | 0.5 |
| 1000 | Forest+GTM(NJ) | FastME | +0.0058 | 0/9/2 | 0.5 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0058 | 0/9/2 | 0.5 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0077 | 0/8/3 | 0.25 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.0097 | 0/7/4 | 0.12 |
| 2000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/11/0 | 1 |
| 2000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/11/0 | 1 |
| 2000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/11/0 | 1 |
| 2000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/11/0 | 1 |
| 2000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/11/0 | 1 |
| 2000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/11/0 | 1 |
| 5000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/10/0 | 1 |
| 5000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/10/0 | 1 |
| 5000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/10/0 | 1 |
| 5000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/10/0 | 1 |
| 5000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/10/0 | 1 |
| 5000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/10/0 | 1 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.45 | 0.00 | 0.00 | 19.4 / 47 |
| 200 | 0.00 | 0.64 | 0.00 | 0.00 | 30.9 / 47 |
| 500 | 0.64 | 0.91 | 0.64 | 0.00 | 37.4 / 47 |
| 1000 | 0.82 | 0.82 | 0.82 | 0.00 | 40.8 / 47 |
| 2000 | 1.00 | 1.00 | 1.00 | 0.00 | 43.2 / 47 |
| 5000 | 1.00 | 1.00 | 1.00 | 0.60 | 46.2 / 47 |

### long U[0.05,0.1] (Kim et al.), n=100

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 10 | 0.00 | 0.094 | 0.085 | 0.041 | 0.027 | 0.669 [0.005] | 0.070 | 0.088 | 0.047 | 0.043 | 0.069 | 0.058 | 19.6 | 0.20 |
| 300 | 10 | 0.00 | 0.008 | 0.008 | 0.003 | 0.000 | 0.398 [0.003] | 0.005 | 0.010 | 0.004 | 0.004 | 0.004 | 0.005 | 9.0 | 0.20 |
| 1000 | 10 | 0.00 | 0.000 | 0.001 | 0.000 | 0.000 | 0.137 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 3.0 | 0.00 |
| 3000 | 10 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.012 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |
| 10000 | 10 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.001 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |
| 100000 | 10 | 0.00 | 0.000 | 0.000 | 0.000 | — | 0.001 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | -0.0237 | 9/1/0 | 0.0039 |
| 100 | Forest+GTM(FastME) | FastME | +0.0062 | 1/3/6 | 0.047 |
| 100 | Forest+GTM(NJ) | FastME | +0.0289 | 0/2/8 | 0.0078 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0010 | 0/9/1 | 1 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0124 | 3/0/7 | 0.27 |
| 100 | Forest+GTM(NJ) | FastTree | +0.0433 | 0/1/9 | 0.0039 |
| 300 | Forest+GTM(NJ) | NJ | -0.0031 | 3/6/1 | 0.5 |
| 300 | Forest+GTM(FastME) | FastME | +0.0010 | 0/9/1 | 1 |
| 300 | Forest+GTM(NJ) | FastME | +0.0021 | 0/8/2 | 0.5 |
| 300 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0010 | 0/9/1 | 1 |
| 300 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 1/8/1 | 1 |
| 300 | Forest+GTM(NJ) | FastTree | +0.0052 | 0/6/4 | 0.12 |
| 1000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/10/0 | 1 |
| 1000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/10/0 | 1 |
| 1000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/10/0 | 1 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/10/0 | 1 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/10/0 | 1 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/10/0 | 1 |
| 3000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/10/0 | 1 |
| 3000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/10/0 | 1 |
| 3000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/10/0 | 1 |
| 3000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/10/0 | 1 |
| 3000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/10/0 | 1 |
| 3000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/10/0 | 1 |
| 10000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/10/0 | 1 |
| 10000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/10/0 | 1 |
| 10000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/10/0 | 1 |
| 10000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/10/0 | 1 |
| 10000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/10/0 | 1 |
| 10000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/10/0 | 1 |
| 100000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/10/0 | 1 |
| 100000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/10/0 | 1 |
| 100000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/10/0 | 1 |
| 100000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/10/0 | 1 |
| 100000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/10/0 | 1 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.80 | 0.00 | 0.00 | 32.3 / 97 |
| 300 | 0.30 | 0.80 | 0.30 | 0.00 | 58.6 / 97 |
| 1000 | 1.00 | 1.00 | 1.00 | 0.10 | 83.7 / 97 |
| 3000 | 1.00 | 1.00 | 1.00 | 0.70 | 95.8 / 97 |
| 10000 | 1.00 | 1.00 | 1.00 | 0.90 | 96.9 / 97 |
| 100000 | 1.00 | 1.00 | 1.00 | 0.90 | 96.9 / 97 |

### deep U[0.1,0.4], n=100

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 10 | 0.37 | 0.414 | 0.419 | 0.331 | 0.199 | 0.898 [0.008] | 0.394 | 0.403 | 0.322 | 0.327 | 0.396 | 0.397 | 41.0 | 0.10 |
| 300 | 10 | 0.32 | 0.296 | 0.259 | 0.101 | 0.019 | 0.763 [0.010] | 0.208 | 0.260 | 0.092 | 0.102 | 0.206 | 0.205 | 23.2 | 0.30 |
| 1000 | 10 | 0.24 | 0.264 | 0.197 | 0.057 | 0.000 | 0.671 [0.006] | 0.145 | 0.206 | 0.054 | 0.057 | 0.152 | 0.098 | 15.3 | 0.30 |
| 3000 | 10 | 0.16 | 0.190 | 0.140 | 0.015 | 0.000 | 0.644 [0.006] | 0.067 | 0.136 | 0.015 | 0.018 | 0.071 | 0.039 | 13.2 | 0.20 |
| 10000 | 10 | 0.10 | 0.126 | 0.102 | 0.003 | 0.000 | 0.426 [0.002] | 0.030 | 0.063 | 0.003 | 0.003 | 0.051 | 0.011 | 6.9 | 0.10 |
| 100000 | 10 | 0.03 | 0.052 | 0.037 | 0.000 | — | 0.243 [0.005] | 0.014 | 0.021 | 0.004 | 0.004 | 0.021 | 0.000 | 4.9 | 0.30 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | -0.0206 | 9/1/0 | 0.0039 |
| 100 | Forest+GTM(FastME) | FastME | -0.0093 | 6/4/0 | 0.031 |
| 100 | Forest+GTM(NJ) | FastME | +0.0629 | 1/1/8 | 0.012 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0021 | 2/8/0 | 0.5 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0031 | 4/1/5 | 0.98 |
| 100 | Forest+GTM(NJ) | FastTree | +0.1948 | 0/0/10 | 0.002 |
| 300 | Forest+GTM(NJ) | NJ | -0.0876 | 10/0/0 | 0.002 |
| 300 | Forest+GTM(FastME) | FastME | -0.0093 | 6/3/1 | 0.11 |
| 300 | Forest+GTM(NJ) | FastME | +0.1072 | 1/0/9 | 0.0059 |
| 300 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0021 | 0/8/2 | 0.5 |
| 300 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0031 | 4/3/3 | 0.84 |
| 300 | Forest+GTM(NJ) | FastTree | +0.1897 | 0/0/10 | 0.002 |
| 1000 | Forest+GTM(NJ) | NJ | -0.1186 | 10/0/0 | 0.002 |
| 1000 | Forest+GTM(FastME) | FastME | -0.0031 | 4/5/1 | 0.5 |
| 1000 | Forest+GTM(NJ) | FastME | +0.0887 | 1/0/9 | 0.02 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0062 | 2/7/1 | 0.5 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0474 | 2/0/8 | 0.014 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.1454 | 0/0/10 | 0.002 |
| 3000 | Forest+GTM(NJ) | NJ | -0.1227 | 10/0/0 | 0.002 |
| 3000 | Forest+GTM(FastME) | FastME | +0.0000 | 2/6/2 | 0.75 |
| 3000 | Forest+GTM(NJ) | FastME | +0.0515 | 0/2/8 | 0.0078 |
| 3000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0041 | 2/7/1 | 0.75 |
| 3000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0278 | 2/0/8 | 0.064 |
| 3000 | Forest+GTM(NJ) | FastTree | +0.0670 | 0/0/10 | 0.002 |
| 10000 | Forest+GTM(NJ) | NJ | -0.0959 | 9/1/0 | 0.0039 |
| 10000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/10/0 | 1 |
| 10000 | Forest+GTM(NJ) | FastME | +0.0268 | 0/3/7 | 0.016 |
| 10000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0206 | 4/6/0 | 0.12 |
| 10000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0186 | 1/2/7 | 0.031 |
| 10000 | Forest+GTM(NJ) | FastTree | +0.0299 | 0/2/8 | 0.0078 |
| 100000 | Forest+GTM(NJ) | NJ | -0.0371 | 8/1/1 | 0.0078 |
| 100000 | Forest+GTM(FastME) | FastME | +0.0041 | 0/7/3 | 0.25 |
| 100000 | Forest+GTM(NJ) | FastME | +0.0144 | 0/6/4 | 0.12 |
| 100000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0062 | 2/6/2 | 1 |
| 100000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0144 | 0/6/4 | 0.12 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.90 | 0.00 | 0.00 | 10.0 / 97 |
| 300 | 0.00 | 0.70 | 0.00 | 0.00 | 23.3 / 97 |
| 1000 | 0.00 | 0.90 | 0.00 | 0.00 | 32.2 / 97 |
| 3000 | 0.00 | 0.80 | 0.00 | 0.00 | 34.7 / 97 |
| 10000 | 0.00 | 0.90 | 0.00 | 0.00 | 55.8 / 97 |
| 100000 | 0.10 | 0.70 | 0.10 | 0.00 | 73.7 / 97 |

### ultrametric h=2, lognormal(1) rates, K2P, n=100

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 10 | 0.56 | 0.918 | 0.851 | 0.868 | 0.486 | 0.910 [0.210] | 0.770 | 0.843 | 0.749 | 0.800 | 0.803 | 0.885 | 21.1 | 2.70 |
| 300 | 10 | 0.46 | 0.882 | 0.776 | 0.821 | 0.303 | 0.871 [0.093] | 0.733 | 0.781 | 0.687 | 0.745 | 0.763 | 0.835 | 19.9 | 1.50 |
| 1000 | 10 | 0.46 | 0.840 | 0.764 | 0.801 | 0.221 | 0.837 [0.028] | 0.691 | 0.732 | 0.661 | 0.714 | 0.711 | 0.789 | 22.4 | 0.40 |
| 3000 | 9 | 0.35 | 0.758 | 0.663 | 0.690 | 0.129 | 0.786 [0.013] | 0.541 | 0.643 | 0.541 | 0.623 | 0.549 | 0.710 | 20.2 | 0.22 |
| 10000 | 9 | 0.26 | 0.712 | 0.617 | 0.550 | 0.087 | 0.751 [0.000] | 0.528 | 0.614 | 0.408 | 0.495 | 0.530 | 0.643 | 23.1 | 0.00 |
| 100000 | 9 | 0.17 | 0.517 | 0.418 | 0.387 | — | 0.656 [0.000] | 0.362 | 0.451 | 0.286 | 0.339 | 0.362 | 0.397 | 17.3 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | -0.1474 | 10/0/0 | 0.002 |
| 100 | Forest+GTM(FastME) | FastME | -0.1186 | 10/0/0 | 0.002 |
| 100 | Forest+GTM(NJ) | FastME | -0.0979 | 8/1/1 | 0.012 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0330 | 7/2/1 | 0.031 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1144 | 10/0/0 | 0.002 |
| 100 | Forest+GTM(NJ) | FastTree | +0.2845 | 0/0/10 | 0.002 |
| 300 | Forest+GTM(NJ) | NJ | -0.1495 | 10/0/0 | 0.002 |
| 300 | Forest+GTM(FastME) | FastME | -0.1340 | 10/0/0 | 0.002 |
| 300 | Forest+GTM(NJ) | FastME | -0.0876 | 8/0/2 | 0.018 |
| 300 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0299 | 5/5/0 | 0.062 |
| 300 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1021 | 10/0/0 | 0.002 |
| 300 | Forest+GTM(NJ) | FastTree | +0.4299 | 0/0/10 | 0.002 |
| 1000 | Forest+GTM(NJ) | NJ | -0.1495 | 10/0/0 | 0.002 |
| 1000 | Forest+GTM(FastME) | FastME | -0.1402 | 10/0/0 | 0.002 |
| 1000 | Forest+GTM(NJ) | FastME | -0.1103 | 8/2/0 | 0.0078 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0206 | 3/6/1 | 0.25 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0979 | 8/0/2 | 0.0098 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.4701 | 0/0/10 | 0.002 |
| 3000 | Forest+GTM(NJ) | NJ | -0.2176 | 9/0/0 | 0.0039 |
| 3000 | Forest+GTM(FastME) | FastME | -0.1489 | 9/0/0 | 0.0039 |
| 3000 | Forest+GTM(NJ) | FastME | -0.1489 | 8/0/1 | 0.0078 |
| 3000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0080 | 1/8/0 | 1 |
| 3000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1695 | 8/1/0 | 0.0078 |
| 3000 | Forest+GTM(NJ) | FastTree | +0.4112 | 0/0/9 | 0.0039 |
| 10000 | Forest+GTM(NJ) | NJ | -0.1844 | 9/0/0 | 0.0039 |
| 10000 | Forest+GTM(FastME) | FastME | -0.1420 | 9/0/0 | 0.0039 |
| 10000 | Forest+GTM(NJ) | FastME | -0.0218 | 5/0/4 | 0.84 |
| 10000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0023 | 1/8/0 | 1 |
| 10000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1145 | 8/1/0 | 0.0078 |
| 10000 | Forest+GTM(NJ) | FastTree | +0.4410 | 0/0/9 | 0.0039 |
| 100000 | Forest+GTM(NJ) | NJ | -0.1546 | 9/0/0 | 0.0039 |
| 100000 | Forest+GTM(FastME) | FastME | -0.1008 | 9/0/0 | 0.0039 |
| 100000 | Forest+GTM(NJ) | FastME | -0.0252 | 5/1/3 | 0.46 |
| 100000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/9/0 | 1 |
| 100000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0355 | 7/0/2 | 0.16 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.20 | 0.00 | 0.00 | 11.4 / 97 |
| 300 | 0.00 | 0.60 | 0.00 | 0.00 | 14.0 / 97 |
| 1000 | 0.00 | 0.70 | 0.00 | 0.00 | 16.2 / 97 |
| 3000 | 0.00 | 0.78 | 0.00 | 0.00 | 21.0 / 97 |
| 10000 | 0.00 | 1.00 | 0.00 | 0.00 | 24.1 / 97 |
| 100000 | 0.00 | 1.00 | 0.00 | 0.00 | 33.3 / 97 |

### Runtime (mean seconds per replicate, single core)

| regime | n | k | NJ | FastME | FastTree | Forest (grid search) | Forest+GTM(NJ) total | Centroid dec.+GTM(NJ) |
|---|---|---|---|---|---|---|---|---|
| U:0.005:0.05 | 50 | 100 | 0.01 | 0.01 | 0.2 | 1.8 | 2.1 | 0.30 |
| U:0.005:0.05 | 50 | 200 | 0.01 | 0.01 | 0.3 | 1.5 | 1.8 | 0.30 |
| U:0.005:0.05 | 50 | 500 | 0.01 | 0.02 | 0.8 | 0.9 | 1.2 | 0.31 |
| U:0.005:0.05 | 50 | 1000 | 0.02 | 0.02 | 1.6 | 0.7 | 1.0 | 0.30 |
| U:0.005:0.05 | 50 | 2000 | 0.02 | 0.02 | 3.4 | 0.6 | 0.9 | 0.27 |
| U:0.005:0.05 | 50 | 5000 | 0.02 | 0.02 | 8.2 | 0.6 | 0.9 | 0.29 |
| U:0.05:0.1 | 100 | 100 | 0.08 | 0.07 | 0.4 | 11.0 | 11.5 | 0.41 |
| U:0.05:0.1 | 100 | 300 | 0.09 | 0.06 | 1.1 | 6.2 | 6.7 | 0.44 |
| U:0.05:0.1 | 100 | 1000 | 0.12 | 0.07 | 3.7 | 3.4 | 3.8 | 0.45 |
| U:0.05:0.1 | 100 | 3000 | 0.17 | 0.12 | 12.8 | 3.4 | 4.1 | 0.55 |
| U:0.05:0.1 | 100 | 10000 | 0.14 | 0.10 | 44.6 | 3.3 | 3.8 | 0.45 |
| U:0.05:0.1 | 100 | 100000 | 0.13 | 0.09 |  | 2.8 | 3.3 | 0.43 |
| U:0.1:0.4 | 100 | 100 | 0.07 | 0.12 | 0.6 | 8.0 | 8.4 | 0.36 |
| U:0.1:0.4 | 100 | 300 | 0.10 | 0.10 | 1.3 | 9.7 | 10.1 | 0.36 |
| U:0.1:0.4 | 100 | 1000 | 0.11 | 0.10 | 4.0 | 11.1 | 11.6 | 0.39 |
| U:0.1:0.4 | 100 | 3000 | 0.10 | 0.08 | 13.5 | 14.2 | 14.6 | 0.45 |
| U:0.1:0.4 | 100 | 10000 | 0.07 | 0.06 | 46.2 | 9.4 | 9.8 | 0.43 |
| U:0.1:0.4 | 100 | 100000 | 0.10 | 0.08 |  | 5.4 | 5.8 | 0.37 |
| UH:2.0:1.0 | 100 | 100 | 0.04 | 0.13 | 1.0 | 8.4 | 8.8 | 0.32 |
| UH:2.0:1.0 | 100 | 300 | 0.05 | 0.13 | 2.4 | 15.2 | 15.5 | 0.31 |
| UH:2.0:1.0 | 100 | 1000 | 0.05 | 0.13 | 6.8 | 17.5 | 17.9 | 0.37 |
| UH:2.0:1.0 | 100 | 3000 | 0.06 | 0.14 | 19.9 | 29.0 | 29.4 | 0.36 |
| UH:2.0:1.0 | 100 | 10000 | 0.04 | 0.11 | 66.1 | 60.2 | 60.6 | 0.35 |
| UH:2.0:1.0 | 100 | 100000 | 0.05 | 0.09 |  | 65.1 | 65.4 | 0.35 |
