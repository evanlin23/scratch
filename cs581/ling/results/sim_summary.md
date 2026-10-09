# Simulation results

60 replicate runs (60 train, 0 test); n = 24 languages.

## Cap selection on training replicates (mean FN pooled over conditions)

| MC | Cap2 | Cap3 | Cap4 | MP | WMC | WCap2 | WCap3 | WCap4 |
|---|---|---|---|---|---|---|---|---|
| 0.0936 | 0.1040 | 0.1024 | 0.1063 | 0.1079 | 0.1079 | 0.1103 | 0.1183 | 0.1167 |

Selected: unweighted **MC**, weighted **WMC** (pilot "new" method = selected cap).

## Mean FN rate on TEST replicates (lower is better)

| method | poly-high | moderate | clean | homoplasy | borrow3 | lexonly | all |
|---|---|---|---|---|---|---|---|
| NJ | - | - | - | - | - | - | nan |
| MP | - | - | - | - | - | - | nan |
| MP-poly | - | - | - | - | - | - | nan |
| MC | - | - | - | - | - | - | nan |
| WMC | - | - | - | - | - | - | nan |
| Cap2 | - | - | - | - | - | - | nan |
| WCap2 | - | - | - | - | - | - | nan |
| Cap3 | - | - | - | - | - | - | nan |
| WCap3 | - | - | - | - | - | - | nan |
| Cap4 | - | - | - | - | - | - | nan |
| WCap4 | - | - | - | - | - | - | nan |
| ML-Mk | - | - | - | - | - | - | nan |
| ML-bin | - | - | - | - | - | - | nan |
| (reps) | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

## Simulated data summaries (all reps; compare with IE: lexical ~14.3 states/char, morph ~9-10, phon 2.3)

| cond | states_L | states_M | states_P | poly_char_frac | poly_cell_frac | borrow events |
|---|---|---|---|---|---|---|
| poly-high | 12.26 | 9.66 | 2.46 | 0.92 | 0.34 | 26.2 |
| moderate | 12.42 | 10.15 | 2.39 | 0.81 | 0.18 | 27.9 |
| clean | 11.44 | 10.24 | 2.90 | 0.00 | 0.00 | 0.0 |
| homoplasy | 12.53 | 8.87 | 2.18 | 0.83 | 0.18 | 25.7 |
| borrow3 | 12.40 | 9.45 | 2.62 | 0.81 | 0.18 | 114.5 |
| lexonly | 11.95 | nan | nan | 0.80 | 0.17 | 50.1 |

## Paired Wilcoxon signed-rank tests on TEST reps (A vs B; mean dFN = A - B; W/T/L = A better/tie/worse)

| A | B | scope | mean dFN | W/T/L | p |
|---|---|---|---|---|---|

## Mean runtime per replicate (seconds, one core)

| NJ | MP | MP-poly | MC | WMC | Cap2 | WCap2 | Cap3 | WCap3 | Cap4 | WCap4 | ML-Mk | ML-bin |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.02 | 0.63 | 0.55 | 0.64 | 0.64 | 0.62 | 0.62 | 0.63 | 0.61 | 0.60 | 0.60 | 64.03 | 4.91 |
