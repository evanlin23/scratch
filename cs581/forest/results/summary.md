Replicates: 47 ok, 0 errors.


### short U[0.005,0.05] (Kim et al.), n=50

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 3 | 0.00 | 0.135 | 0.135 | 0.106 | 0.106 | 0.574 [0.040] | 0.142 | 0.142 | 0.121 | 0.113 | 0.135 | 0.135 | 6.3 | 0.67 |
| 200 | 2 | 0.00 | 0.021 | 0.021 | 0.000 | 0.011 | 0.362 [0.015] | 0.021 | 0.032 | 0.011 | 0.011 | 0.011 | 0.000 | 1.5 | 0.50 |
| 500 | 2 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.287 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.011 | 1.0 | 0.00 |
| 1000 | 2 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.074 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |
| 2000 | 2 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.128 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |
| 5000 | 2 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.043 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | +0.0071 | 0/2/1 | 1 |
| 100 | Forest+GTM(FastME) | FastME | +0.0142 | 0/2/1 | 1 |
| 100 | Forest+GTM(NJ) | FastME | +0.0355 | 0/1/2 | 0.5 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0071 | 0/2/1 | 1 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0071 | 1/1/1 | 1 |
| 100 | Forest+GTM(NJ) | FastTree | +0.0355 | 1/0/2 | 0.5 |
| 200 | Forest+GTM(NJ) | NJ | +0.0000 | 1/0/1 | 1 |
| 200 | Forest+GTM(FastME) | FastME | +0.0106 | 0/1/1 | 1 |
| 200 | Forest+GTM(NJ) | FastME | +0.0213 | 0/1/1 | 1 |
| 200 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0106 | 0/1/1 | 1 |
| 200 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0213 | 0/1/1 | 1 |
| 200 | Forest+GTM(NJ) | FastTree | +0.0106 | 1/0/1 | 1 |
| 500 | Forest+GTM(NJ) | NJ | +0.0000 | 0/2/0 | 1 |
| 500 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 500 | Forest+GTM(NJ) | FastME | +0.0000 | 0/2/0 | 1 |
| 500 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 500 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0106 | 1/1/0 | 1 |
| 500 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/2/0 | 1 |
| 1000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/2/0 | 1 |
| 1000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 1000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/2/0 | 1 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/2/0 | 1 |
| 2000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/2/0 | 1 |
| 2000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 2000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/2/0 | 1 |
| 2000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 2000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 2000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/2/0 | 1 |
| 5000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/2/0 | 1 |
| 5000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 5000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/2/0 | 1 |
| 5000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 5000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 5000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/2/0 | 1 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.33 | 0.00 | 0.00 | 20.7 / 47 |
| 200 | 0.00 | 0.50 | 0.00 | 0.00 | 30.5 / 47 |
| 500 | 1.00 | 1.00 | 1.00 | 0.00 | 33.5 / 47 |
| 1000 | 1.00 | 1.00 | 1.00 | 0.00 | 43.5 / 47 |
| 2000 | 1.00 | 1.00 | 1.00 | 0.00 | 41.0 / 47 |
| 5000 | 1.00 | 1.00 | 1.00 | 0.50 | 45.0 / 47 |

### long U[0.05,0.1] (Kim et al.), n=100

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 2 | 0.00 | 0.098 | 0.103 | 0.021 | 0.031 | 0.660 [0.000] | 0.057 | 0.082 | 0.026 | 0.021 | 0.057 | 0.052 | 17.5 | 0.00 |
| 300 | 2 | 0.00 | 0.005 | 0.005 | 0.000 | 0.000 | 0.593 [0.000] | 0.005 | 0.005 | 0.000 | 0.000 | 0.005 | 0.000 | 17.5 | 0.00 |
| 1000 | 2 | 0.00 | 0.000 | 0.005 | 0.000 | 0.000 | 0.309 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 10.0 | 0.00 |
| 3000 | 2 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |
| 10000 | 2 | 0.00 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |
| 100000 | 2 | 0.00 | 0.000 | 0.000 | 0.000 | — | 0.000 [0.000] | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.0 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | -0.0412 | 2/0/0 | 0.5 |
| 100 | Forest+GTM(FastME) | FastME | +0.0052 | 0/1/1 | 1 |
| 100 | Forest+GTM(NJ) | FastME | +0.0361 | 0/0/2 | 0.5 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0052 | 1/0/1 | 1 |
| 100 | Forest+GTM(NJ) | FastTree | +0.0258 | 0/1/1 | 1 |
| 300 | Forest+GTM(NJ) | NJ | +0.0000 | 0/2/0 | 1 |
| 300 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 300 | Forest+GTM(NJ) | FastME | +0.0052 | 0/1/1 | 1 |
| 300 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 300 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0052 | 0/1/1 | 1 |
| 300 | Forest+GTM(NJ) | FastTree | +0.0052 | 0/1/1 | 1 |
| 1000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/2/0 | 1 |
| 1000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 1000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/2/0 | 1 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/2/0 | 1 |
| 3000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/2/0 | 1 |
| 3000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 3000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/2/0 | 1 |
| 3000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 3000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 3000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/2/0 | 1 |
| 10000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/2/0 | 1 |
| 10000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 10000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/2/0 | 1 |
| 10000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 10000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 10000 | Forest+GTM(NJ) | FastTree | +0.0000 | 0/2/0 | 1 |
| 100000 | Forest+GTM(NJ) | NJ | +0.0000 | 0/2/0 | 1 |
| 100000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 100000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/2/0 | 1 |
| 100000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 100000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/2/0 | 1 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 1.00 | 0.00 | 0.00 | 33.0 / 97 |
| 300 | 0.50 | 1.00 | 0.50 | 0.00 | 39.5 / 97 |
| 1000 | 1.00 | 1.00 | 1.00 | 0.00 | 67.0 / 97 |
| 3000 | 1.00 | 1.00 | 1.00 | 1.00 | 97.0 / 97 |
| 10000 | 1.00 | 1.00 | 1.00 | 1.00 | 97.0 / 97 |
| 100000 | 1.00 | 1.00 | 1.00 | 1.00 | 97.0 / 97 |

### deep U[0.1,0.4], n=100

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 2 | 0.36 | 0.448 | 0.474 | 0.320 | 0.211 | 0.938 [0.000] | 0.433 | 0.443 | 0.314 | 0.314 | 0.438 | 0.361 | 47.0 | 0.00 |
| 300 | 2 | 0.29 | 0.247 | 0.196 | 0.046 | 0.015 | 0.773 [0.018] | 0.186 | 0.227 | 0.046 | 0.052 | 0.180 | 0.160 | 25.5 | 0.50 |
| 1000 | 2 | 0.26 | 0.314 | 0.263 | 0.134 | 0.000 | 0.716 [0.000] | 0.170 | 0.242 | 0.124 | 0.124 | 0.191 | 0.077 | 16.0 | 0.00 |
| 3000 | 2 | 0.16 | 0.196 | 0.155 | 0.010 | 0.000 | 0.675 [0.000] | 0.119 | 0.160 | 0.010 | 0.010 | 0.119 | 0.052 | 21.5 | 0.00 |
| 10000 | 2 | 0.10 | 0.088 | 0.088 | 0.005 | 0.000 | 0.515 [0.000] | 0.057 | 0.057 | 0.005 | 0.005 | 0.057 | 0.015 | 14.0 | 0.00 |
| 100000 | 2 | 0.02 | 0.067 | 0.052 | 0.000 | — | 0.294 [0.000] | 0.000 | 0.021 | 0.000 | 0.000 | 0.000 | 0.000 | 7.5 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | -0.0155 | 2/0/0 | 0.5 |
| 100 | Forest+GTM(FastME) | FastME | -0.0052 | 1/1/0 | 1 |
| 100 | Forest+GTM(NJ) | FastME | +0.1134 | 0/0/2 | 0.5 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0052 | 1/1/0 | 1 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0722 | 0/0/2 | 0.5 |
| 100 | Forest+GTM(NJ) | FastTree | +0.2216 | 0/0/2 | 0.5 |
| 300 | Forest+GTM(NJ) | NJ | -0.0619 | 2/0/0 | 0.5 |
| 300 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 300 | Forest+GTM(NJ) | FastME | +0.1392 | 0/0/2 | 0.5 |
| 300 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0052 | 0/1/1 | 1 |
| 300 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0258 | 0/0/2 | 0.5 |
| 300 | Forest+GTM(NJ) | FastTree | +0.1701 | 0/0/2 | 0.5 |
| 1000 | Forest+GTM(NJ) | NJ | -0.1443 | 2/0/0 | 0.5 |
| 1000 | Forest+GTM(FastME) | FastME | -0.0103 | 2/0/0 | 0.5 |
| 1000 | Forest+GTM(NJ) | FastME | +0.0361 | 1/0/1 | 1 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0206 | 1/1/0 | 1 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0928 | 0/0/2 | 0.5 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.1701 | 0/0/2 | 0.5 |
| 3000 | Forest+GTM(NJ) | NJ | -0.0773 | 2/0/0 | 0.5 |
| 3000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 3000 | Forest+GTM(NJ) | FastME | +0.1082 | 0/0/2 | 0.5 |
| 3000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 3000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0670 | 0/0/2 | 0.5 |
| 3000 | Forest+GTM(NJ) | FastTree | +0.1186 | 0/0/2 | 0.5 |
| 10000 | Forest+GTM(NJ) | NJ | -0.0309 | 1/1/0 | 1 |
| 10000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 10000 | Forest+GTM(NJ) | FastME | +0.0515 | 0/0/2 | 0.5 |
| 10000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 10000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0412 | 0/0/2 | 0.5 |
| 10000 | Forest+GTM(NJ) | FastTree | +0.0567 | 0/0/2 | 0.5 |
| 100000 | Forest+GTM(NJ) | NJ | -0.0670 | 2/0/0 | 0.5 |
| 100000 | Forest+GTM(FastME) | FastME | +0.0000 | 0/2/0 | 1 |
| 100000 | Forest+GTM(NJ) | FastME | +0.0000 | 0/2/0 | 1 |
| 100000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/2/0 | 1 |
| 100000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | +0.0000 | 0/2/0 | 1 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 1.00 | 0.00 | 0.00 | 6.0 / 97 |
| 300 | 0.00 | 0.50 | 0.00 | 0.00 | 22.5 / 97 |
| 1000 | 0.00 | 1.00 | 0.00 | 0.00 | 27.5 / 97 |
| 3000 | 0.00 | 1.00 | 0.00 | 0.00 | 31.5 / 97 |
| 10000 | 0.00 | 1.00 | 0.00 | 0.00 | 47.0 / 97 |
| 100000 | 0.00 | 1.00 | 0.00 | 0.00 | 68.5 / 97 |

### ultrametric h=2, lognormal(1) rates, K2P, n=100

Mean FN rate (FP rate for Forest in brackets; Forest FN = 1 - correct splits/(n-3)).

| k | reps | sat | NJ | BIONJ | FastME | FastTree | Forest | Forest+GTM(NJ) | Forest+GTM(NJ, induced) | Forest+GTM(FastME) | Forest+GTM(FastME, induced) | Forest comps only+GTM(NJ) | Centroid dec.+GTM(NJ) | Forest #comp | Forest false splits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 2 | 0.62 | 0.923 | 0.840 | 0.918 | 0.500 | 0.918 [0.056] | 0.809 | 0.876 | 0.825 | 0.840 | 0.845 | 0.902 | 27.5 | 0.50 |
| 300 | 2 | 0.46 | 0.820 | 0.747 | 0.763 | 0.304 | 0.876 [0.136] | 0.670 | 0.722 | 0.629 | 0.665 | 0.691 | 0.722 | 14.0 | 1.50 |
| 1000 | 2 | 0.47 | 0.861 | 0.758 | 0.861 | 0.237 | 0.892 [0.056] | 0.737 | 0.784 | 0.696 | 0.768 | 0.763 | 0.851 | 21.0 | 0.50 |
| 3000 | 2 | 0.40 | 0.784 | 0.722 | 0.763 | 0.139 | 0.866 [0.038] | 0.613 | 0.665 | 0.557 | 0.660 | 0.680 | 0.763 | 15.0 | 0.50 |
| 10000 | 1 | 0.18 | 0.722 | 0.629 | 0.526 | 0.113 | 0.825 [0.000] | 0.567 | 0.649 | 0.402 | 0.505 | 0.567 | 0.649 | 26.0 | 0.00 |
| 100000 | 1 | 0.22 | 0.639 | 0.515 | 0.619 | — | 0.639 [0.000] | 0.515 | 0.567 | 0.464 | 0.526 | 0.515 | 0.536 | 21.0 | 0.00 |

Paired two-sided Wilcoxon on FN rate (A − B; negative = A better). W/T/L = A better / tie / A worse, tie band |diff| < 0.5/(n−3) (less than one split).

| k | A | B | mean diff | W/T/L | p |
|---|---|---|---|---|---|
| 100 | Forest+GTM(NJ) | NJ | -0.1134 | 2/0/0 | 0.5 |
| 100 | Forest+GTM(FastME) | FastME | -0.0928 | 2/0/0 | 0.5 |
| 100 | Forest+GTM(NJ) | FastME | -0.1082 | 2/0/0 | 0.5 |
| 100 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0361 | 1/1/0 | 1 |
| 100 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0928 | 2/0/0 | 0.5 |
| 100 | Forest+GTM(NJ) | FastTree | +0.3093 | 0/0/2 | 0.5 |
| 300 | Forest+GTM(NJ) | NJ | -0.1495 | 2/0/0 | 0.5 |
| 300 | Forest+GTM(FastME) | FastME | -0.1340 | 2/0/0 | 0.5 |
| 300 | Forest+GTM(NJ) | FastME | -0.0928 | 2/0/0 | 0.5 |
| 300 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0206 | 1/1/0 | 1 |
| 300 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0515 | 2/0/0 | 0.5 |
| 300 | Forest+GTM(NJ) | FastTree | +0.3660 | 0/0/2 | 0.5 |
| 1000 | Forest+GTM(NJ) | NJ | -0.1237 | 2/0/0 | 0.5 |
| 1000 | Forest+GTM(FastME) | FastME | -0.1649 | 2/0/0 | 0.5 |
| 1000 | Forest+GTM(NJ) | FastME | -0.1237 | 2/0/0 | 0.5 |
| 1000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0258 | 1/1/0 | 1 |
| 1000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1134 | 2/0/0 | 0.5 |
| 1000 | Forest+GTM(NJ) | FastTree | +0.5000 | 0/0/2 | 0.5 |
| 3000 | Forest+GTM(NJ) | NJ | -0.1701 | 2/0/0 | 0.5 |
| 3000 | Forest+GTM(FastME) | FastME | -0.2062 | 2/0/0 | 0.5 |
| 3000 | Forest+GTM(NJ) | FastME | -0.1495 | 2/0/0 | 0.5 |
| 3000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | -0.0670 | 1/1/0 | 1 |
| 3000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.1495 | 2/0/0 | 0.5 |
| 3000 | Forest+GTM(NJ) | FastTree | +0.4742 | 0/0/2 | 0.5 |
| 10000 | Forest+GTM(NJ) | NJ | -0.1546 | 1/0/0 | 1 |
| 10000 | Forest+GTM(FastME) | FastME | -0.1237 | 1/0/0 | 1 |
| 10000 | Forest+GTM(NJ) | FastME | +0.0412 | 0/0/1 | 1 |
| 10000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/1/0 | 1 |
| 10000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0825 | 1/0/0 | 1 |
| 10000 | Forest+GTM(NJ) | FastTree | +0.4536 | 0/0/1 | 1 |
| 100000 | Forest+GTM(NJ) | NJ | -0.1237 | 1/0/0 | 1 |
| 100000 | Forest+GTM(FastME) | FastME | -0.1546 | 1/0/0 | 1 |
| 100000 | Forest+GTM(NJ) | FastME | -0.1031 | 1/0/0 | 1 |
| 100000 | Forest+GTM(NJ) | Forest comps only+GTM(NJ) | +0.0000 | 0/1/0 | 1 |
| 100000 | Forest+GTM(NJ) | Centroid dec.+GTM(NJ) | -0.0206 | 1/0/0 | 1 |

Kim et al. style success rates: P(no incorrect split) and P(fully resolved and correct).

| k | NJ no-false | Forest no-false | NJ fully correct | Forest fully correct | Forest mean #splits |
|---|---|---|---|---|---|
| 100 | 0.00 | 0.50 | 0.00 | 0.00 | 8.5 / 97 |
| 300 | 0.00 | 0.50 | 0.00 | 0.00 | 13.5 / 97 |
| 1000 | 0.00 | 0.50 | 0.00 | 0.00 | 11.0 / 97 |
| 3000 | 0.00 | 0.50 | 0.00 | 0.00 | 13.5 / 97 |
| 10000 | 0.00 | 1.00 | 0.00 | 0.00 | 17.0 / 97 |
| 100000 | 0.00 | 1.00 | 0.00 | 0.00 | 35.0 / 97 |

### Runtime (mean seconds per replicate, single core)

| regime | n | k | NJ | FastME | FastTree | Forest (grid search) | Forest+GTM(NJ) total | Centroid dec.+GTM(NJ) |
|---|---|---|---|---|---|---|---|---|
| U:0.005:0.05 | 50 | 100 | 0.01 | 0.01 | 0.1 | 0.6 | 0.8 | 0.21 |
| U:0.005:0.05 | 50 | 200 | 0.02 | 0.01 | 0.2 | 0.4 | 0.5 | 0.15 |
| U:0.005:0.05 | 50 | 500 | 0.01 | 0.01 | 0.4 | 0.3 | 0.5 | 0.20 |
| U:0.005:0.05 | 50 | 1000 | 0.01 | 0.01 | 0.8 | 0.2 | 0.4 | 0.14 |
| U:0.005:0.05 | 50 | 2000 | 0.00 | 0.00 | 1.6 | 0.2 | 0.3 | 0.13 |
| U:0.005:0.05 | 50 | 5000 | 0.01 | 0.01 | 4.2 | 0.2 | 0.4 | 0.14 |
| U:0.05:0.1 | 100 | 100 | 0.06 | 0.04 | 0.3 | 4.0 | 4.2 | 0.26 |
| U:0.05:0.1 | 100 | 300 | 0.07 | 0.05 | 0.5 | 2.4 | 2.7 | 0.26 |
| U:0.05:0.1 | 100 | 1000 | 0.13 | 0.06 | 1.9 | 1.7 | 2.0 | 0.27 |
| U:0.05:0.1 | 100 | 3000 | 0.08 | 0.05 | 6.6 | 1.3 | 1.6 | 0.23 |
| U:0.05:0.1 | 100 | 10000 | 0.05 | 0.03 | 23.1 | 1.5 | 1.7 | 0.19 |
| U:0.05:0.1 | 100 | 100000 | 0.05 | 0.05 |  | 1.1 | 1.3 | 0.19 |
| U:0.1:0.4 | 100 | 100 | 0.03 | 0.05 | 0.2 | 4.3 | 4.5 | 0.17 |
| U:0.1:0.4 | 100 | 300 | 0.05 | 0.07 | 0.7 | 4.5 | 4.7 | 0.20 |
| U:0.1:0.4 | 100 | 1000 | 0.08 | 0.07 | 2.0 | 3.8 | 4.0 | 0.23 |
| U:0.1:0.4 | 100 | 3000 | 0.05 | 0.04 | 6.4 | 5.7 | 5.9 | 0.20 |
| U:0.1:0.4 | 100 | 10000 | 0.02 | 0.02 | 23.1 | 3.6 | 3.8 | 0.18 |
| U:0.1:0.4 | 100 | 100000 | 0.02 | 0.02 |  | 2.3 | 2.5 | 0.16 |
| UH:2.0:1.0 | 100 | 100 | 0.02 | 0.07 | 0.5 | 2.3 | 2.5 | 0.23 |
| UH:2.0:1.0 | 100 | 300 | 0.02 | 0.10 | 1.1 | 4.0 | 4.2 | 0.17 |
| UH:2.0:1.0 | 100 | 1000 | 0.02 | 0.07 | 3.6 | 4.2 | 4.4 | 0.17 |
| UH:2.0:1.0 | 100 | 3000 | 0.05 | 0.06 | 10.6 | 8.9 | 9.1 | 0.29 |
| UH:2.0:1.0 | 100 | 10000 | 0.01 | 0.03 | 32.4 | 31.1 | 31.3 | 0.17 |
| UH:2.0:1.0 | 100 | 100000 | 0.02 | 0.08 |  | 23.1 | 23.4 | 0.16 |
