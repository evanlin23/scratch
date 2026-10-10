# Supplement: binary encoding vs polymorphism information (TEST reps)

Replicates regenerated from the same seeds; NJ FN reproduced exactly in 120/120 replicates.

| method | clean | moderate | poly-high | homoplasy | borrow3 | lexonly | all |
|---|---|---|---|---|---|---|---|
| MC | 0.057 | 0.079 | 0.136 | 0.129 | 0.131 | 0.124 | 0.109 |
| MP | 0.057 | 0.100 | 0.140 | 0.143 | 0.152 | 0.133 | 0.121 |
| ML-Mk | 0.086 | 0.117 | 0.157 | 0.157 | 0.138 | 0.148 | 0.134 |
| ML-bin | 0.083 | 0.076 | 0.086 | 0.086 | 0.095 | 0.093 | 0.087 |
| ML-bin-res | 0.083 | 0.093 | 0.143 | 0.117 | 0.112 | 0.129 | 0.113 |
| MP-bin | 0.129 | 0.098 | 0.114 | 0.138 | 0.131 | 0.133 | 0.124 |
| MP-bin-res | 0.126 | 0.090 | 0.152 | 0.155 | 0.150 | 0.143 | 0.136 |

| A | B | mean dFN | W/T/L | p (Wilcoxon, exact ties dropped) |
|---|---|---|---|---|
| ML-bin | ML-bin-res | -0.0262 | 57/50/13 | 3.83e-06 |
| ML-bin-res | MC | +0.0036 | 36/38/46 | 0.798 |
| ML-bin-res | ML-Mk | -0.0210 | 56/29/35 | 0.0051 |
| MP-bin | MP | +0.0028 | 46/25/49 | 0.78 |
| MP-bin | ML-bin | +0.0373 | 20/37/63 | 3.57e-07 |
| MP-bin-res | MP | +0.0151 | 36/27/57 | 0.052 |
| MP-bin | MC | +0.0147 | 39/24/57 | 0.105 |
