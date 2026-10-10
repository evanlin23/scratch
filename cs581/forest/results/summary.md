Replicates: 480 ok, 0 errors.


### short U[0.005,0.05] (Kim et al.), n=50

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 20 | 0.00 | 0.153 | 0.147 | 0.141 | 0.129 | 0.570 [0.053] | 0.165 | 0.163 | 0.153 | 0.154 | 0.155 | 0.159 | 5.3 | 1.20 |
| 200 | 20 | 0.00 | 0.077 | 0.073 | 0.061 | 0.051 | 0.372 [0.011] | 0.076 | 0.079 | 0.066 | 0.066 | 0.073 | 0.070 | 1.8 | 0.35 |
| 500 | 20 | 0.00 | 0.015 | 0.013 | 0.013 | 0.009 | 0.210 [0.003] | 0.016 | 0.016 | 0.014 | 0.014 | 0.015 | 0.013 | 1.4 | 0.10 |
| 1000 | 20 | 0.00 | 0.006 | 0.004 | 0.004 | 0.000 | 0.141 [0.007] | 0.013 | 0.013 | 0.011 | 0.011 | 0.006 | 0.003 | 1.0 | 0.30 |
| 2000 | 20 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.079 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.1 | 0.00 |
| 5000 | 20 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.023 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | +0.0117 | 1/10/9 | 0.0084 |
| 100 | Forest+GTM(FastME) | FastME | +0.0117 | 2/11/7 | 0.044 |
| 100 | Forest+GTM(NJ) | FastME | +0.0234 | 1/7/12 | 0.0018 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0096 | 0/14/6 | 0.026 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0064 | 5/7/8 | 0.3 |
| 100 | Forest+GTM(NJ) | FastTree | +0.0362 | 2/3/15 | 0.0013 |
| 200 | Forest+GTM(NJ) | NJ | -0.0011 | 3/15/2 | 0.5 |
| 200 | Forest+GTM(FastME) | FastME | +0.0053 | 0/15/5 | 0.038 |
| 200 | Forest+GTM(NJ) | FastME | +0.0149 | 3/6/11 | 0.011 |
| 200 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0021 | 0/18/2 | 0.18 |
| 200 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0053 | 4/7/9 | 0.6 |
| 200 | Forest+GTM(NJ) | FastTree | +0.0245 | 4/4/12 | 0.015 |
| 500 | Forest+GTM(NJ) | NJ | +0.0011 | 0/19/1 | 0.32 |
| 500 | Forest+GTM(FastME) | FastME | +0.0011 | 0/19/1 | 0.32 |
| 500 | Forest+GTM(NJ) | FastME | +0.0032 | 0/17/3 | 0.1 |
| 500 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0011 | 0/19/1 | 0.32 |
| 500 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0032 | 1/15/4 | 0.16 |
| 500 | Forest+GTM(NJ) | FastTree | +0.0074 | 0/16/4 | 0.059 |
| 1000 | Forest+GTM(NJ) | NJ | +0.0064 | 0/16/4 | 0.063 |
| 1000 | Forest+GTM(FastME) | FastME | +0.0064 | 0/16/4 | 0.063 |
| 1000 | Forest+GTM(NJ) | FastME | +0.0085 | 0/14/6 | 0.023 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0064 | 0/16/4 | 0.063 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0096 | 0/13/7 | 0.014 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.0128 | 0/12/8 | 0.0097 |
| 2000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/20/0 | 1 |
| 2000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/20/0 | 1 |
| 2000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/20/0 | 1 |
| 2000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/20/0 | 1 |
| 2000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/20/0 | 1 |
| 2000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/20/0 | 1 |
| 5000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/20/0 | 1 |
| 5000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/20/0 | 1 |
| 5000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/20/0 | 1 |
| 5000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/20/0 | 1 |
| 5000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/20/0 | 1 |
| 5000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/20/0 | 1 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.40 | 0.00 | 0.00 | 21.4 / 47 |
| 200 | 0.00 | 0.65 | 0.00 | 0.00 | 29.9 / 47 |
| 500 | 0.50 | 0.95 | 0.50 | 0.00 | 37.2 / 47 |
| 1000 | 0.75 | 0.80 | 0.75 | 0.00 | 40.6 / 47 |
| 2000 | 1.00 | 1.00 | 1.00 | 0.00 | 43.3 / 47 |
| 5000 | 1.00 | 1.00 | 1.00 | 0.45 | 45.9 / 47 |

### long U[0.05,0.1] (Kim et al.), n=100

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 20 | 0.00 | 0.099 | 0.093 | 0.038 | 0.021 | 0.641 [0.006] | 0.068 | 0.087 | 0.042 | 0.039 | 0.066 | 0.053 | 17.6 | 0.25 |
| 300 | 20 | 0.00 | 0.007 | 0.006 | 0.002 | 0.000 | 0.366 [0.002] | 0.004 | 0.007 | 0.003 | 0.002 | 0.003 | 0.004 | 8.2 | 0.10 |
| 1000 | 20 | 0.00 | 0.000 | 0.001 | 0.000 | 0.000 | 0.126 [0.001] | 0.001 | 0.001 | 0.001 | 0.001 | 0.000 | 0.000 | 3.1 | 0.05 |
| 3000 | 20 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.015 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |
| 10000 | 20 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.003 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |
| 100000 | 20 | 0.00 | 0.000 | 0.000 | 0.000 | — | 0.009 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.2 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | -0.0314 | 18/2/0 | 0.00019 |
| 100 | Forest+GTM(FastME) | FastME | +0.0041 | 3/9/8 | 0.09 |
| 100 | Forest+GTM(NJ) | FastME | +0.0299 | 0/4/16 | 0.00043 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0015 | 0/18/2 | 0.18 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0144 | 4/2/14 | 0.02 |
| 100 | Forest+GTM(NJ) | FastTree | +0.0464 | 0/1/19 | 0.00013 |
| 300 | Forest+GTM(NJ) | NJ | -0.0031 | 6/13/1 | 0.058 |
| 300 | Forest+GTM(FastME) | FastME | +0.0010 | 0/18/2 | 0.16 |
| 300 | Forest+GTM(NJ) | FastME | +0.0021 | 0/16/4 | 0.046 |
| 300 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0005 | 0/19/1 | 0.32 |
| 300 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 3/14/3 | 1 |
| 300 | Forest+GTM(NJ) | FastTree | +0.0036 | 0/14/6 | 0.02 |
| 1000 | Forest+GTM(NJ) | NJ | +0.0005 | 0/19/1 | 0.32 |
| 1000 | Forest+GTM(FastME) | FastME | +0.0005 | 0/19/1 | 0.32 |
| 1000 | Forest+GTM(NJ) | FastME | +0.0005 | 0/19/1 | 0.32 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0005 | 0/19/1 | 0.32 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0005 | 0/19/1 | 0.32 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.0005 | 0/19/1 | 0.32 |
| 3000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/20/0 | 1 |
| 3000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/20/0 | 1 |
| 3000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/20/0 | 1 |
| 3000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/20/0 | 1 |
| 3000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/20/0 | 1 |
| 3000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/20/0 | 1 |
| 10000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/20/0 | 1 |
| 10000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/20/0 | 1 |
| 10000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/20/0 | 1 |
| 10000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/20/0 | 1 |
| 10000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/20/0 | 1 |
| 10000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/20/0 | 1 |
| 100000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/20/0 | 1 |
| 100000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/20/0 | 1 |
| 100000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/20/0 | 1 |
| 100000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/20/0 | 1 |
| 100000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/20/0 | 1 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.80 | 0.00 | 0.00 | 35.0 / 97 |
| 300 | 0.40 | 0.90 | 0.40 | 0.00 | 61.5 / 97 |
| 1000 | 1.00 | 0.95 | 1.00 | 0.05 | 84.8 / 97 |
| 3000 | 1.00 | 1.00 | 1.00 | 0.60 | 95.5 / 97 |
| 10000 | 1.00 | 1.00 | 1.00 | 0.80 | 96.7 / 97 |
| 100000 | 1.00 | 1.00 | 1.00 | 0.85 | 96.2 / 97 |

### deep U[0.1,0.4], n=100

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 20 | 0.39 | 0.430 | 0.430 | 0.318 | 0.193 | 0.908 [0.006] | 0.409 | 0.419 | 0.312 | 0.315 | 0.410 | 0.403 | 42.3 | 0.10 |
| 300 | 20 | 0.32 | 0.297 | 0.286 | 0.119 | 0.023 | 0.786 [0.007] | 0.215 | 0.265 | 0.112 | 0.119 | 0.216 | 0.211 | 22.9 | 0.20 |
| 1000 | 20 | 0.23 | 0.241 | 0.185 | 0.043 | 0.001 | 0.665 [0.005] | 0.135 | 0.190 | 0.041 | 0.044 | 0.138 | 0.098 | 16.4 | 0.20 |
| 3000 | 20 | 0.17 | 0.204 | 0.151 | 0.014 | 0.000 | 0.608 [0.004] | 0.084 | 0.143 | 0.014 | 0.016 | 0.088 | 0.053 | 13.4 | 0.15 |
| 10000 | 20 | 0.09 | 0.110 | 0.083 | 0.003 | 0.000 | 0.426 [0.003] | 0.028 | 0.059 | 0.003 | 0.003 | 0.044 | 0.013 | 8.1 | 0.15 |
| 100000 | 20 | 0.03 | 0.046 | 0.029 | 0.000 | — | 0.237 [0.003] | 0.015 | 0.018 | 0.003 | 0.003 | 0.022 | 0.002 | 4.5 | 0.20 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | -0.0211 | 16/4/0 | 0.0004 |
| 100 | Forest+GTM(FastME) | FastME | -0.0057 | 10/8/2 | 0.059 |
| 100 | Forest+GTM(NJ) | FastME | +0.0907 | 2/1/17 | 0.00067 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0010 | 3/16/1 | 0.45 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0057 | 7/2/11 | 0.57 |
| 100 | Forest+GTM(NJ) | FastTree | +0.2155 | 0/0/20 | 8.8e-05 |
| 300 | Forest+GTM(NJ) | NJ | -0.0814 | 20/0/0 | 8.8e-05 |
| 300 | Forest+GTM(FastME) | FastME | -0.0072 | 11/6/3 | 0.032 |
| 300 | Forest+GTM(NJ) | FastME | +0.0964 | 3/1/16 | 0.00083 |
| 300 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0005 | 3/15/2 | 0.49 |
| 300 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0046 | 8/3/9 | 0.54 |
| 300 | Forest+GTM(NJ) | FastTree | +0.1928 | 0/0/20 | 8.7e-05 |
| 1000 | Forest+GTM(NJ) | NJ | -0.1062 | 20/0/0 | 8.8e-05 |
| 1000 | Forest+GTM(FastME) | FastME | -0.0026 | 7/11/2 | 0.17 |
| 1000 | Forest+GTM(NJ) | FastME | +0.0918 | 1/0/19 | 0.00054 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0031 | 3/15/2 | 0.28 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0366 | 4/0/16 | 0.0032 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.1345 | 0/0/20 | 8.5e-05 |
| 3000 | Forest+GTM(NJ) | NJ | -0.1201 | 20/0/0 | 8.8e-05 |
| 3000 | Forest+GTM(FastME) | FastME | +0.0000 | 3/14/3 | 0.59 |
| 3000 | Forest+GTM(NJ) | FastME | +0.0691 | 0/2/18 | 0.00019 |
| 3000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0041 | 6/13/1 | 0.15 |
| 3000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0304 | 3/1/16 | 0.0058 |
| 3000 | Forest+GTM(NJ) | FastTree | +0.0835 | 0/0/20 | 8.7e-05 |
| 10000 | Forest+GTM(NJ) | NJ | -0.0814 | 19/1/0 | 0.00013 |
| 10000 | Forest+GTM(FastME) | FastME | +0.0005 | 0/19/1 | 0.32 |
| 10000 | Forest+GTM(NJ) | FastME | +0.0258 | 0/4/16 | 0.00042 |
| 10000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0160 | 6/13/1 | 0.063 |
| 10000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0149 | 3/4/13 | 0.0065 |
| 10000 | Forest+GTM(NJ) | FastTree | +0.0284 | 0/3/17 | 0.00028 |
| 100000 | Forest+GTM(NJ) | NJ | -0.0314 | 17/2/1 | 0.00038 |
| 100000 | Forest+GTM(FastME) | FastME | +0.0031 | 0/16/4 | 0.063 |
| 100000 | Forest+GTM(NJ) | FastME | +0.0149 | 0/12/8 | 0.011 |
| 100000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0072 | 6/12/2 | 0.16 |
| 100000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0134 | 1/11/8 | 0.015 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.90 | 0.00 | 0.00 | 9.0 / 97 |
| 300 | 0.00 | 0.80 | 0.00 | 0.00 | 21.0 / 97 |
| 1000 | 0.00 | 0.90 | 0.00 | 0.00 | 32.7 / 97 |
| 3000 | 0.00 | 0.85 | 0.00 | 0.00 | 38.2 / 97 |
| 10000 | 0.00 | 0.85 | 0.00 | 0.00 | 55.8 / 97 |
| 100000 | 0.10 | 0.80 | 0.10 | 0.00 | 74.2 / 97 |

### ultrametric h=2, lognormal(1) rates, K2P, n=100

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 20 | 0.53 | 0.913 | 0.854 | 0.880 | 0.484 | 0.910 [0.184] | 0.775 | 0.827 | 0.761 | 0.804 | 0.829 | 0.889 | 15.6 | 2.45 |
| 300 | 20 | 0.48 | 0.891 | 0.811 | 0.849 | 0.319 | 0.893 [0.163] | 0.755 | 0.796 | 0.726 | 0.769 | 0.798 | 0.862 | 15.2 | 2.05 |
| 1000 | 20 | 0.44 | 0.847 | 0.756 | 0.777 | 0.216 | 0.849 [0.030] | 0.689 | 0.745 | 0.642 | 0.710 | 0.716 | 0.795 | 21.8 | 0.40 |
| 3000 | 20 | 0.35 | 0.764 | 0.675 | 0.666 | 0.135 | 0.786 [0.016] | 0.562 | 0.654 | 0.523 | 0.601 | 0.569 | 0.698 | 22.1 | 0.25 |
| 10000 | 20 | 0.23 | 0.693 | 0.614 | 0.537 | 0.097 | 0.740 [0.010] | 0.504 | 0.590 | 0.405 | 0.488 | 0.507 | 0.603 | 21.9 | 0.20 |
| 100000 | 20 | 0.19 | 0.571 | 0.470 | 0.425 | — | 0.649 [0.000] | 0.408 | 0.482 | 0.319 | 0.371 | 0.411 | 0.454 | 19.4 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | -0.1381 | 20/0/0 | 8.8e-05 |
| 100 | Forest+GTM(FastME) | FastME | -0.1186 | 20/0/0 | 8.7e-05 |
| 100 | Forest+GTM(NJ) | FastME | -0.1052 | 17/1/2 | 0.00023 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0541 | 15/4/1 | 0.00062 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1144 | 19/0/1 | 0.0001 |
| 100 | Forest+GTM(NJ) | FastTree | +0.2907 | 0/0/20 | 8.8e-05 |
| 300 | Forest+GTM(NJ) | NJ | -0.1366 | 20/0/0 | 8.8e-05 |
| 300 | Forest+GTM(FastME) | FastME | -0.1232 | 20/0/0 | 8.8e-05 |
| 300 | Forest+GTM(NJ) | FastME | -0.0948 | 17/0/3 | 0.00077 |
| 300 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0433 | 13/7/0 | 0.0014 |
| 300 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1072 | 20/0/0 | 8.7e-05 |
| 300 | Forest+GTM(NJ) | FastTree | +0.4356 | 0/0/20 | 8.8e-05 |
| 1000 | Forest+GTM(NJ) | NJ | -0.1588 | 20/0/0 | 8.7e-05 |
| 1000 | Forest+GTM(FastME) | FastME | -0.1351 | 20/0/0 | 8.8e-05 |
| 1000 | Forest+GTM(NJ) | FastME | -0.0887 | 15/3/2 | 0.0019 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0273 | 7/12/1 | 0.017 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1062 | 16/0/4 | 0.00051 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.4722 | 0/0/20 | 8.8e-05 |
| 3000 | Forest+GTM(NJ) | NJ | -0.2021 | 20/0/0 | 8.8e-05 |
| 3000 | Forest+GTM(FastME) | FastME | -0.1433 | 20/0/0 | 8.8e-05 |
| 3000 | Forest+GTM(NJ) | FastME | -0.1041 | 16/1/3 | 0.00054 |
| 3000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0062 | 2/18/0 | 0.18 |
| 3000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1361 | 18/1/1 | 0.00015 |
| 3000 | Forest+GTM(NJ) | FastTree | +0.4273 | 0/0/20 | 8.8e-05 |
| 10000 | Forest+GTM(NJ) | NJ | -0.1892 | 20/0/0 | 8.8e-05 |
| 10000 | Forest+GTM(FastME) | FastME | -0.1325 | 20/0/0 | 1.9e-06 |
| 10000 | Forest+GTM(NJ) | FastME | -0.0335 | 12/0/8 | 0.34 |
| 10000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0031 | 2/17/1 | 0.29 |
| 10000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0990 | 17/1/2 | 0.00031 |
| 10000 | Forest+GTM(NJ) | FastTree | +0.4067 | 0/0/20 | 8.8e-05 |
| 100000 | Forest+GTM(NJ) | NJ | -0.1629 | 20/0/0 | 8.8e-05 |
| 100000 | Forest+GTM(FastME) | FastME | -0.1067 | 20/0/0 | 8.8e-05 |
| 100000 | Forest+GTM(NJ) | FastME | -0.0170 | 10/1/9 | 0.34 |
| 100000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0031 | 1/19/0 | 0.32 |
| 100000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0454 | 17/0/3 | 0.0076 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.30 | 0.00 | 0.00 | 11.2 / 97 |
| 300 | 0.00 | 0.35 | 0.00 | 0.00 | 12.4 / 97 |
| 1000 | 0.00 | 0.65 | 0.00 | 0.00 | 15.1 / 97 |
| 3000 | 0.00 | 0.80 | 0.00 | 0.00 | 21.0 / 97 |
| 10000 | 0.00 | 0.90 | 0.00 | 0.00 | 25.4 / 97 |
| 100000 | 0.00 | 1.00 | 0.00 | 0.00 | 34.0 / 97 |

### Runtime (mean seconds per replicate, single core)

| regime | n | k | NJ | FastME | FastTree | Forest (grid search) | Forest+GTM(NJ) total | Centroid dec.+GTM(NJ) |
|---|---|---|---|---|---|---|---|---|
| U:0.005:0.05 | 50 | 100 | 0.01 | 0.01 | 0.2 | 1.8 | 2.2 | 0.33 |
| U:0.005:0.05 | 50 | 200 | 0.01 | 0.01 | 0.4 | 1.4 | 1.7 | 0.30 |
| U:0.005:0.05 | 50 | 500 | 0.01 | 0.02 | 0.8 | 0.8 | 1.1 | 0.28 |
| U:0.005:0.05 | 50 | 1000 | 0.01 | 0.02 | 1.6 | 0.6 | 1.0 | 0.29 |
| U:0.005:0.05 | 50 | 2000 | 0.02 | 0.01 | 3.3 | 0.6 | 0.9 | 0.28 |
| U:0.005:0.05 | 50 | 5000 | 0.02 | 0.02 | 8.3 | 0.6 | 0.9 | 0.28 |
| U:0.05:0.1 | 100 | 100 | 0.07 | 0.06 | 0.4 | 11.4 | 11.8 | 0.38 |
| U:0.05:0.1 | 100 | 300 | 0.10 | 0.07 | 1.1 | 5.8 | 6.3 | 0.42 |
| U:0.05:0.1 | 100 | 1000 | 0.12 | 0.08 | 3.8 | 3.3 | 3.8 | 0.44 |
| U:0.05:0.1 | 100 | 3000 | 0.14 | 0.10 | 12.8 | 3.3 | 3.9 | 0.49 |
| U:0.05:0.1 | 100 | 10000 | 0.13 | 0.09 | 44.2 | 3.1 | 3.5 | 0.44 |
| U:0.05:0.1 | 100 | 100000 | 0.12 | 0.08 |  | 2.8 | 3.3 | 0.40 |
| U:0.1:0.4 | 100 | 100 | 0.06 | 0.11 | 0.6 | 8.6 | 9.0 | 0.36 |
| U:0.1:0.4 | 100 | 300 | 0.09 | 0.10 | 1.3 | 9.9 | 10.3 | 0.37 |
| U:0.1:0.4 | 100 | 1000 | 0.08 | 0.08 | 4.0 | 10.7 | 11.2 | 0.39 |
| U:0.1:0.4 | 100 | 3000 | 0.10 | 0.08 | 13.3 | 12.9 | 13.3 | 0.45 |
| U:0.1:0.4 | 100 | 10000 | 0.06 | 0.05 | 45.3 | 10.1 | 10.5 | 0.38 |
| U:0.1:0.4 | 100 | 100000 | 0.09 | 0.07 |  | 5.1 | 5.5 | 0.37 |
| UH:2.0:1.0 | 100 | 100 | 0.04 | 0.14 | 1.1 | 13.7 | 14.0 | 0.34 |
| UH:2.0:1.0 | 100 | 300 | 0.05 | 0.14 | 2.3 | 12.8 | 13.2 | 0.33 |
| UH:2.0:1.0 | 100 | 1000 | 0.05 | 0.14 | 6.5 | 18.9 | 19.3 | 0.35 |
| UH:2.0:1.0 | 100 | 3000 | 0.05 | 0.12 | 19.6 | 32.6 | 32.9 | 0.32 |
| UH:2.0:1.0 | 100 | 10000 | 0.04 | 0.11 | 66.2 | 66.3 | 66.6 | 0.35 |
| UH:2.0:1.0 | 100 | 100000 | 0.05 | 0.10 |  | 63.6 | 63.9 | 0.33 |

### Large trees (n=500; reduced grid; 20 replicates, 0 errors)

| regime | k | reps | NJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(FastME) | Centroid dec.+GTM(NJ) | Forest #comp | t_Forest (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| U:0.1:0.4 | 300 | 5 | 0.301 | 0.219 | 0.038 | 1.000 | 0.301 | 0.219 | 0.263 | 421 | 130 |
| U:0.1:0.4 | 3000 | 5 | 0.246 | 0.048 | 0.010 | 0.976 | 0.241 | 0.048 | 0.106 | 294 | 161 |
| UH:2.0:1.0 | 300 | 5 | 0.887 | 0.745 | 0.277 | 0.940 | 0.672 | 0.617 | 0.810 | 130 | 128 |
| UH:2.0:1.0 | 3000 | 5 | 0.817 | 0.626 | 0.133 | 0.923 | 0.630 | 0.507 | 0.727 | 200 | 321 |

Paired tests (n=500, all k pooled):

| A | B | mean diff | W/T/L | p | N |
|---|---|---|---|---|---|
| Forest+GTM(NJ) | NJ | -0.1017 | 14/6/0 | 0.00098 | 20 |
| Forest+GTM(FastME) | FastME | -0.0617 | 9/11/0 | 0.0076 | 20 |
| Forest+GTM(NJ) | FastME | +0.0515 | 8/0/12 | 0.19 | 20 |
| Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0016 | 2/18/0 | 0.18 | 20 |
| Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0156 | 10/0/10 | 0.62 | 20 |
| Forest+GTM(NJ) | FastTree | +0.3465 | 0/0/20 | 8.8e-05 | 20 |

### Follow-up: saturation handling and FastME-guided controls (60 replicates, 0 errors)

| regime | k | reps | NJ_cap2 | NJ_cap1.2 | NJ_cap5 | NJ_pcap | FastME_cap2 | FastME_cap1.2 | FastME_cap5 | FastME_pcap | FGTM_FastME | CompGTM_FastME | DecGTM_FastME_25 | DecGTM_FastME_50 | FGTM_FastME_pcap |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| U:0.1:0.4 | 300 | 10 | 0.296 | 0.121 | 0.748 | 0.212 | 0.101 | 0.037 | 0.335 | 0.098 | 0.092 | 0.090 | 0.109 | 0.112 | 0.097 |
| U:0.1:0.4 | 3000 | 10 | 0.190 | 0.060 | 0.531 | 0.067 | 0.015 | 0.001 | 0.053 | 0.002 | 0.015 | 0.013 | 0.011 | 0.009 | 0.004 |
| U:0.1:0.4 | 100000 | 10 | 0.052 | 0.009 | 0.144 | 0.002 | 0.000 | 0.000 | 0.001 | 0.000 | 0.004 | 0.000 | 0.000 | 0.000 | 0.004 |
| UH:2.0:1.0 | 300 | 10 | 0.882 | 0.834 | 0.921 | 0.431 | 0.821 | 0.765 | 0.869 | 0.373 | 0.687 | 0.711 | 0.797 | 0.813 | 0.390 |
| UH:2.0:1.0 | 3000 | 10 | 0.753 | 0.706 | 0.808 | 0.315 | 0.666 | 0.589 | 0.701 | 0.251 | 0.519 | 0.522 | 0.631 | 0.651 | 0.238 |
| UH:2.0:1.0 | 100000 | 10 | 0.513 | 0.423 | 0.598 | 0.194 | 0.373 | 0.330 | 0.448 | 0.145 | 0.275 | 0.275 | 0.334 | 0.364 | 0.118 |

| regime | k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|
| U:0.1:0.4 | 300 | FGTM_FastME_pcap | FastME_pcap | -0.0010 | 2/7/1 | 0.75 |
| U:0.1:0.4 | 300 | FGTM_FastME | FastME_cap2 | -0.0093 | 6/3/1 | 0.11 |
| U:0.1:0.4 | 300 | FGTM_FastME | CompGTM_FastME | +0.0021 | 0/8/2 | 0.5 |
| U:0.1:0.4 | 300 | FGTM_FastME | DecGTM_FastME_25 | -0.0175 | 7/2/1 | 0.07 |
| U:0.1:0.4 | 300 | FastME_pcap | FastME_cap2 | -0.0031 | 4/2/4 | 0.91 |
| U:0.1:0.4 | 3000 | FGTM_FastME_pcap | FastME_pcap | +0.0021 | 0/8/2 | 0.5 |
| U:0.1:0.4 | 3000 | FGTM_FastME | FastME_cap2 | +0.0000 | 2/6/2 | 0.75 |
| U:0.1:0.4 | 3000 | FGTM_FastME | CompGTM_FastME | +0.0021 | 0/8/2 | 0.5 |
| U:0.1:0.4 | 3000 | FGTM_FastME | DecGTM_FastME_25 | +0.0041 | 0/6/4 | 0.12 |
| U:0.1:0.4 | 3000 | FastME_pcap | FastME_cap2 | -0.0134 | 8/1/1 | 0.023 |
| U:0.1:0.4 | 100000 | FGTM_FastME_pcap | FastME_pcap | +0.0041 | 0/7/3 | 0.25 |
| U:0.1:0.4 | 100000 | FGTM_FastME | FastME_cap2 | +0.0041 | 0/7/3 | 0.25 |
| U:0.1:0.4 | 100000 | FGTM_FastME | CompGTM_FastME | +0.0041 | 0/7/3 | 0.25 |
| U:0.1:0.4 | 100000 | FGTM_FastME | DecGTM_FastME_25 | +0.0041 | 0/7/3 | 0.25 |
| U:0.1:0.4 | 100000 | FastME_pcap | FastME_cap2 | +0.0000 | 0/10/0 | 1 |
| UH:2.0:1.0 | 300 | FGTM_FastME_pcap | FastME_pcap | +0.0165 | 4/2/4 | 0.31 |
| UH:2.0:1.0 | 300 | FGTM_FastME | FastME_cap2 | -0.1340 | 10/0/0 | 0.002 |
| UH:2.0:1.0 | 300 | FGTM_FastME | CompGTM_FastME | -0.0247 | 5/4/1 | 0.094 |
| UH:2.0:1.0 | 300 | FGTM_FastME | DecGTM_FastME_25 | -0.1103 | 10/0/0 | 0.002 |
| UH:2.0:1.0 | 300 | FastME_pcap | FastME_cap2 | -0.4474 | 10/0/0 | 0.002 |
| UH:2.0:1.0 | 3000 | FGTM_FastME_pcap | FastME_pcap | -0.0124 | 6/2/2 | 0.12 |
| UH:2.0:1.0 | 3000 | FGTM_FastME | FastME_cap2 | -0.1474 | 10/0/0 | 0.002 |
| UH:2.0:1.0 | 3000 | FGTM_FastME | CompGTM_FastME | -0.0031 | 1/8/1 | 1 |
| UH:2.0:1.0 | 3000 | FGTM_FastME | DecGTM_FastME_25 | -0.1124 | 10/0/0 | 0.002 |
| UH:2.0:1.0 | 3000 | FastME_pcap | FastME_cap2 | -0.4155 | 10/0/0 | 0.002 |
| UH:2.0:1.0 | 100000 | FGTM_FastME_pcap | FastME_pcap | -0.0278 | 9/0/1 | 0.0039 |
| UH:2.0:1.0 | 100000 | FGTM_FastME | FastME_cap2 | -0.0979 | 10/0/0 | 0.002 |
| UH:2.0:1.0 | 100000 | FGTM_FastME | CompGTM_FastME | +0.0000 | 0/10/0 | 1 |
| UH:2.0:1.0 | 100000 | FGTM_FastME | DecGTM_FastME_25 | -0.0588 | 10/0/0 | 0.002 |
| UH:2.0:1.0 | 100000 | FastME_pcap | FastME_cap2 | -0.2278 | 10/0/0 | 0.002 |
