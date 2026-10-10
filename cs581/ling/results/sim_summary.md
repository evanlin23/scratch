# Simulation results

180 replicate runs (60 train, 120 test); n = 24 languages.

## Cap selection on training replicates (mean FN pooled over conditions)

| MC | Cap2 | Cap3 | Cap4 | MP | WMC | WCap2 | WCap3 | WCap4 |
|---|---|---|---|---|---|---|---|---|
| 0.0936 | 0.1040 | 0.1024 | 0.1063 | 0.1079 | 0.1079 | 0.1103 | 0.1183 | 0.1167 |

Selected: unweighted **MC**, weighted **WMC** (pilot "new" method = selected cap).

## Mean FN rate on TEST replicates (lower is better)

| method | poly-high | moderate | clean | homoplasy | borrow3 | lexonly | all |
|---|---|---|---|---|---|---|---|
| NJ | 0.131 | 0.083 | 0.098 | 0.143 | 0.133 | 0.129 | 0.119 |
| MP | 0.140 | 0.100 | 0.057 | 0.143 | 0.152 | 0.133 | 0.121 |
| MP-poly | 0.190 | 0.129 | 0.055 | 0.212 | 0.171 | 0.212 | 0.162 |
| MC | 0.136 | 0.079 | 0.057 | 0.129 | 0.131 | 0.124 | 0.109 |
| WMC | 0.143 | 0.079 | 0.057 | 0.140 | 0.140 | 0.124 | 0.114 |
| Cap2 | 0.148 | 0.095 | 0.057 | 0.133 | 0.136 | 0.133 | 0.117 |
| WCap2 | 0.148 | 0.098 | 0.057 | 0.167 | 0.145 | 0.131 | 0.124 |
| Cap3 | 0.143 | 0.095 | 0.055 | 0.152 | 0.143 | 0.140 | 0.121 |
| WCap3 | 0.155 | 0.105 | 0.057 | 0.174 | 0.164 | 0.143 | 0.133 |
| Cap4 | 0.143 | 0.090 | 0.057 | 0.152 | 0.138 | 0.140 | 0.120 |
| WCap4 | 0.169 | 0.110 | 0.055 | 0.181 | 0.162 | 0.143 | 0.137 |
| ML-Mk | 0.157 | 0.117 | 0.086 | 0.157 | 0.138 | 0.148 | 0.134 |
| ML-bin | 0.086 | 0.076 | 0.083 | 0.086 | 0.095 | 0.093 | 0.087 |
| (reps) | 20 | 20 | 20 | 20 | 20 | 20 | 120 |

## Simulated data summaries (all reps; compare with IE: lexical ~14.3 states/char, morph ~9-10, phon 2.3)

| cond | states_L | states_M | states_P | poly_char_frac | poly_cell_frac | borrow events |
|---|---|---|---|---|---|---|
| poly-high | 12.55 | 9.38 | 2.55 | 0.92 | 0.35 | 26.6 |
| moderate | 12.34 | 9.75 | 2.57 | 0.81 | 0.18 | 27.0 |
| clean | 11.47 | 9.82 | 2.75 | 0.00 | 0.00 | 0.0 |
| homoplasy | 11.84 | 8.66 | 2.19 | 0.81 | 0.18 | 26.5 |
| borrow3 | 12.41 | 9.46 | 2.66 | 0.81 | 0.18 | 117.6 |
| lexonly | 11.94 | nan | nan | 0.80 | 0.17 | 49.0 |

## Paired Wilcoxon signed-rank tests on TEST reps (A vs B; mean dFN = A - B; W/T/L = A better/tie/worse)

| A | B | scope | mean dFN | W/T/L | p |
|---|---|---|---|---|---|
| MC | MP | all | -0.0119 | 39/64/17 | 0.00692 |
| MC | MP | poly-high | -0.0048 | 7/8/5 | 0.594 |
| MC | MP | moderate | -0.0214 | 9/9/2 | 0.0869 |
| MC | MP | clean | +0.0000 | 0/20/0 | nan |
| MC | MP | homoplasy | -0.0143 | 9/7/4 | 0.281 |
| MC | MP | borrow3 | -0.0214 | 9/10/1 | 0.0176 |
| MC | MP | lexonly | -0.0095 | 5/10/5 | 0.607 |
| MC | ML-Mk | all | -0.0246 | 58/36/26 | 0.000963 |
| MC | ML-Mk | poly-high | -0.0214 | 10/4/6 | 0.182 |
| MC | ML-Mk | moderate | -0.0381 | 12/4/4 | 0.0351 |
| MC | ML-Mk | clean | -0.0286 | 10/7/3 | 0.0344 |
| MC | ML-Mk | homoplasy | -0.0286 | 9/8/3 | 0.116 |
| MC | ML-Mk | borrow3 | -0.0072 | 7/7/6 | 0.564 |
| MC | ML-Mk | lexonly | -0.0238 | 10/6/4 | 0.282 |
| MC | ML-bin | all | +0.0226 | 36/29/55 | 0.00504 |
| MC | ML-bin | poly-high | +0.0500 | 3/6/11 | 0.0174 |
| MC | ML-bin | moderate | +0.0024 | 9/4/7 | 0.875 |
| MC | ML-bin | clean | -0.0262 | 8/8/4 | 0.238 |
| MC | ML-bin | homoplasy | +0.0429 | 5/4/11 | 0.0455 |
| MC | ML-bin | borrow3 | +0.0357 | 5/3/12 | 0.128 |
| MC | ML-bin | lexonly | +0.0310 | 6/4/10 | 0.253 |
| MC | NJ | all | -0.0103 | 49/31/40 | 0.256 |
| MC | NJ | poly-high | +0.0048 | 8/2/10 | 0.965 |
| MC | NJ | moderate | -0.0048 | 8/5/7 | 0.977 |
| MC | NJ | clean | -0.0405 | 9/11/0 | 0.00391 |
| MC | NJ | homoplasy | -0.0143 | 9/1/10 | 0.49 |
| MC | NJ | borrow3 | -0.0024 | 7/8/5 | 1 |
| MC | NJ | lexonly | -0.0048 | 8/4/8 | 0.697 |
| WMC | MC | all | +0.0048 | 9/100/11 | 0.12 |
| WMC | MC | poly-high | +0.0071 | 3/14/3 | 0.656 |
| WMC | MC | moderate | +0.0000 | 3/16/1 | nan |
| WMC | MC | clean | +0.0000 | 0/20/0 | nan |
| WMC | MC | homoplasy | +0.0119 | 1/16/3 | nan |
| WMC | MC | borrow3 | +0.0095 | 2/14/4 | 0.469 |
| WMC | MC | lexonly | +0.0000 | 0/20/0 | nan |
| MP-poly | MP | all | +0.0405 | 20/39/61 | 3.29e-06 |
| MP-poly | MP | poly-high | +0.0500 | 5/2/13 | 0.0569 |
| MP-poly | MP | moderate | +0.0286 | 5/4/11 | 0.207 |
| MP-poly | MP | clean | -0.0024 | 1/19/0 | nan |
| MP-poly | MP | homoplasy | +0.0690 | 1/4/15 | 0.00127 |
| MP-poly | MP | borrow3 | +0.0191 | 6/5/9 | 0.169 |
| MP-poly | MP | lexonly | +0.0786 | 2/5/13 | 0.00687 |
| ML-Mk | ML-bin | all | +0.0472 | 22/26/72 | 6.68e-08 |
| ML-Mk | ML-bin | poly-high | +0.0714 | 3/4/13 | 0.00461 |
| ML-Mk | ML-bin | moderate | +0.0405 | 3/4/13 | 0.0141 |
| ML-Mk | ML-bin | clean | +0.0024 | 8/5/7 | 0.931 |
| ML-Mk | ML-bin | homoplasy | +0.0714 | 1/4/15 | 0.000784 |
| ML-Mk | ML-bin | borrow3 | +0.0429 | 2/5/13 | 0.0294 |
| ML-Mk | ML-bin | lexonly | +0.0548 | 5/4/11 | 0.0643 |
| MP | ML-Mk | all | -0.0127 | 53/36/31 | 0.0658 |
| MP | ML-Mk | poly-high | -0.0167 | 10/5/5 | 0.372 |
| MP | ML-Mk | moderate | -0.0167 | 10/5/5 | 0.276 |
| MP | ML-Mk | clean | -0.0286 | 10/7/3 | 0.0344 |
| MP | ML-Mk | homoplasy | -0.0143 | 10/4/6 | 0.621 |
| MP | ML-Mk | borrow3 | +0.0143 | 6/7/7 | 0.692 |
| MP | ML-Mk | lexonly | -0.0143 | 7/8/5 | 0.537 |
| MP | NJ | all | +0.0016 | 41/37/42 | 0.921 |
| MP | NJ | poly-high | +0.0095 | 7/3/10 | 0.867 |
| MP | NJ | moderate | +0.0167 | 6/6/8 | 0.203 |
| MP | NJ | clean | -0.0405 | 9/11/0 | 0.00391 |
| MP | NJ | homoplasy | +0.0000 | 8/4/8 | 0.917 |
| MP | NJ | borrow3 | +0.0191 | 4/7/9 | 0.359 |
| MP | NJ | lexonly | +0.0048 | 7/6/7 | 0.727 |
| Cap2 | MC | all | +0.0079 | 12/77/31 | 0.0394 |
| Cap2 | MC | poly-high | +0.0119 | 4/9/7 | 0.186 |
| Cap2 | MC | moderate | +0.0167 | 1/12/7 | 0.141 |
| Cap2 | MC | clean | +0.0000 | 0/20/0 | nan |
| Cap2 | MC | homoplasy | +0.0048 | 2/13/5 | 0.797 |
| Cap2 | MC | borrow3 | +0.0048 | 3/12/5 | 1 |
| Cap2 | MC | lexonly | +0.0095 | 2/11/7 | 0.312 |
| Cap2 | MP | all | -0.0040 | 30/69/21 | 0.303 |
| Cap2 | MP | poly-high | +0.0071 | 6/7/7 | 0.555 |
| Cap2 | MP | moderate | -0.0048 | 4/13/3 | 0.828 |
| Cap2 | MP | clean | +0.0000 | 0/20/0 | nan |
| Cap2 | MP | homoplasy | -0.0095 | 9/6/5 | 0.369 |
| Cap2 | MP | borrow3 | -0.0167 | 6/13/1 | 0.0625 |
| Cap2 | MP | lexonly | -0.0000 | 5/10/5 | 0.9 |
| Cap3 | MP | all | +0.0004 | 14/91/15 | 0.638 |
| Cap3 | MP | poly-high | +0.0024 | 4/12/4 | 0.906 |
| Cap3 | MP | moderate | -0.0048 | 1/19/0 | nan |
| Cap3 | MP | clean | -0.0024 | 1/19/0 | nan |
| Cap3 | MP | homoplasy | +0.0095 | 1/14/5 | 0.406 |
| Cap3 | MP | borrow3 | -0.0095 | 4/15/1 | 0.125 |
| Cap3 | MP | lexonly | +0.0071 | 3/12/5 | 0.727 |

## Mean runtime per replicate (seconds, one core)

| NJ | MP | MP-poly | MC | WMC | Cap2 | WCap2 | Cap3 | WCap3 | Cap4 | WCap4 | ML-Mk | ML-bin |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.02 | 0.60 | 0.54 | 0.62 | 0.63 | 0.61 | 0.62 | 0.60 | 0.60 | 0.59 | 0.60 | 63.76 | 5.01 |
