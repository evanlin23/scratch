
### 1000M2-HF (simulated, our fragmentation) — true alignment

| method | n | mean FN | Δ FN vs RAxML-NG (W/T/L) | Wilcoxon p | mean CPU min | mean ΔlnL vs RAxML-NG |
|---|---|---|---|---|---|---|
| FastTree 2 | 3 | 48.1% | +23.9 (0/0/3) | 0.25 | 2 | – |
| IQ-TREE 3 --fast | 3 | 34.6% | +10.5 (0/0/3) | 0.25 | 7 | – |
| RAxML-NG --fast | 3 | 26.9% | +2.7 (1/0/2) | 0.50 | 16 | -141.4 |
| RAxML-NG (1 start) | 3 | 24.2% | – | - | 40 | +0.0 |

### 16S.M, real 16S rRNA, random fragments; reference = CRW tree (47% resolved) — true alignment

| method | n | mean FN | Δ FN vs RAxML-NG (W/T/L) | Wilcoxon p | mean CPU min | mean ΔlnL vs RAxML-NG |
|---|---|---|---|---|---|---|
| FastTree 2 | 5 | 52.6% | +19.9 (0/0/5) | 0.06 | 2 | – |
| IQ-TREE 3 --fast | 5 | 39.3% | +6.7 (0/0/5) | 0.06 | 6 | – |
| RAxML-NG --fast | 5 | 34.2% | +1.5 (0/1/4) | 0.12 | 8 | -130.5 |
| RAxML-NG (1 start) | 5 | 32.7% | – | - | 28 | +0.0 |
| **ours-fast** (constrained) | 5 | 39.2% | +6.6 (0/0/5) | 0.06 | 8 | -994.1 |
| **ours-accurate** (+ fast polish) | 5 | 33.9% | +1.2 (1/0/4) | 0.19 | 16 | -96.5 |

### 16S.M, real 16S rRNA, amplicon-like fragments; reference = CRW tree — true alignment

| method | n | mean FN | Δ FN vs RAxML-NG (W/T/L) | Wilcoxon p | mean CPU min | mean ΔlnL vs RAxML-NG |
|---|---|---|---|---|---|---|
| FastTree 2 | 5 | 35.4% | +7.7 (0/0/5) | 0.06 | 2 | – |
| IQ-TREE 3 --fast | 5 | 33.9% | +6.2 (0/0/5) | 0.06 | 5 | – |
| RAxML-NG --fast | 5 | 27.0% | -0.7 (3/1/1) | 0.62 | 8 | -45.5 |
| RAxML-NG (1 start) | 5 | 27.6% | – | - | 30 | +0.0 |
| **ours-fast** (constrained) | 5 | 29.2% | +1.5 (0/2/3) | 0.25 | 8 | -456.1 |
| **ours-accurate** (+ fast polish) | 5 | 26.7% | -1.0 (4/0/1) | 0.44 | 16 | -22.7 |

### 1000M1-HF, Park et al. 2021 published inputs (simulated) — true alignment

| method | n | mean FN | Δ FN vs RAxML-NG (W/T/L) | Wilcoxon p | mean CPU min | mean ΔlnL vs RAxML-NG |
|---|---|---|---|---|---|---|
| FastTree 2 | 5 | 48.9% | +24.2 (0/0/5) | 0.06 | 2 | – |
| IQ-TREE 3 --fast | 5 | 37.5% | +12.8 (0/0/5) | 0.06 | 7 | – |
| RAxML-NG --fast | 5 | 29.8% | +5.1 (0/0/5) | 0.06 | 17 | -246.2 |
| RAxML-NG (1 start) | 5 | 24.7% | – | - | 35 | +0.0 |
| **ours-fast** (constrained) | 5 | 25.7% | +0.9 (0/0/5) | 0.06 | 14 | -23.5 |
| **ours-accurate** (+ fast polish) | 5 | 25.1% | +0.4 (1/1/3) | 0.62 | 24 | +1.4 |
| published RAxML-NG (20 starts, 24 h cap) | 5 | 24.9% | +0.2 (1/1/3) | 0.88 | – | – |
| published IQ-TREE 2 (default) | 5 | 30.2% | +5.4 (0/0/5) | 0.06 | – | – |
| published GTM (IQ-TREE start) | 5 | 28.4% | +3.6 (0/0/5) | 0.06 | – | – |

### 1000M1-HF, Park et al. 2021 published inputs (simulated) — UPP alignment

| method | n | mean FN | Δ FN vs RAxML-NG (W/T/L) | Wilcoxon p | mean CPU min | mean ΔlnL vs RAxML-NG |
|---|---|---|---|---|---|---|
| FastTree 2 | 1 | 91.6% | -0.8 (1/0/0) | - | 2 | – |
| RAxML-NG --fast | 1 | 91.8% | -0.6 (1/0/0) | - | 19 | -1166.5 |
| RAxML-NG (1 start) | 1 | 92.4% | – | - | 43 | +0.0 |
| **ours-fast** (constrained) | 1 | 90.8% | -1.6 (1/0/0) | - | 16 | -767.0 |
| **ours-accurate** (+ fast polish) | 1 | 91.1% | -1.3 (1/0/0) | - | 29 | -62.2 |

UPP alignment cost (not included above): 103 CPU min per replicate.

### RNASim1K-HF (simulated, our fragmentation) — true alignment

| method | n | mean FN | Δ FN vs RAxML-NG (W/T/L) | Wilcoxon p | mean CPU min | mean ΔlnL vs RAxML-NG |
|---|---|---|---|---|---|---|
| FastTree 2 | 3 | 59.3% | +24.0 (0/0/3) | 0.25 | 3 | – |
| IQ-TREE 3 --fast | 3 | 46.1% | +10.8 (0/0/3) | 0.25 | 6 | – |
| RAxML-NG --fast | 3 | 35.4% | +0.1 (2/0/1) | 1.00 | 19 | +0.7 |
| RAxML-NG (1 start) | 3 | 35.3% | – | - | 42 | +0.0 |

### Pooled over all held-out replicates (paired vs RAxML-NG, 1 start)

| method | n pairs | mean Δ FN | W/T/L | Wilcoxon p | mean CPU ratio vs RAxML-NG |
|---|---|---|---|---|---|
| FastTree 2 | 22 | +18.27 | 1/0/21 | 9.54e-07 | 0.07 |
| IQ-TREE 3 --fast | 21 | +9.15 | 0/0/21 | 5.94e-05 | 0.18 |
| RAxML-NG --fast | 22 | +1.70 | 7/2/13 | 0.0111 | 0.37 |
| **ours-fast** (constrained) | 16 | +2.71 | 1/2/13 | 0.00428 | 0.33 |
| **ours-accurate** (+ fast polish) | 16 | +0.13 | 7/1/8 | 0.733 | 0.60 |
