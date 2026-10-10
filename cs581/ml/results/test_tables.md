
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

### 16S.M, real 16S rRNA, amplicon-like fragments; reference = CRW tree — true alignment

| method | n | mean FN | Δ FN vs RAxML-NG (W/T/L) | Wilcoxon p | mean CPU min | mean ΔlnL vs RAxML-NG |
|---|---|---|---|---|---|---|
| FastTree 2 | 5 | 35.4% | +7.7 (0/0/5) | 0.06 | 2 | – |
| IQ-TREE 3 --fast | 5 | 33.9% | +6.2 (0/0/5) | 0.06 | 5 | – |
| RAxML-NG --fast | 5 | 27.0% | -0.7 (3/1/1) | 0.62 | 8 | -45.5 |
| RAxML-NG (1 start) | 5 | 27.6% | – | - | 30 | +0.0 |

### 1000M1-HF, Park et al. 2021 published inputs (simulated) — true alignment

| method | n | mean FN | Δ FN vs RAxML-NG (W/T/L) | Wilcoxon p | mean CPU min | mean ΔlnL vs RAxML-NG |
|---|---|---|---|---|---|---|
| FastTree 2 | 5 | 48.9% | +24.2 (0/0/5) | 0.06 | 2 | – |
| IQ-TREE 3 --fast | 5 | 37.5% | +12.8 (0/0/5) | 0.06 | 7 | – |
| RAxML-NG --fast | 5 | 29.8% | +5.1 (0/0/5) | 0.06 | 17 | -246.2 |
| RAxML-NG (1 start) | 5 | 24.7% | – | - | 35 | +0.0 |
| published RAxML-NG (20 starts, 24 h cap) | 5 | 24.9% | +0.2 (1/1/3) | 0.88 | – | – |
| published IQ-TREE 2 (default) | 5 | 30.2% | +5.4 (0/0/5) | 0.06 | – | – |
| published GTM (IQ-TREE start) | 5 | 28.4% | +3.6 (0/0/5) | 0.06 | – | – |

### RNASim1K-HF (simulated, our fragmentation) — true alignment

| method | n | mean FN | Δ FN vs RAxML-NG (W/T/L) | Wilcoxon p | mean CPU min | mean ΔlnL vs RAxML-NG |
|---|---|---|---|---|---|---|
| FastTree 2 | 3 | 59.3% | +24.9 (0/0/1) | - | 3 | – |
| IQ-TREE 3 --fast | 3 | 46.1% | +13.9 (0/0/1) | - | 6 | – |
| RAxML-NG --fast | 2 | 36.0% | -0.3 (1/0/0) | - | 20 | +3.8 |
| RAxML-NG (1 start) | 1 | 34.0% | – | - | 44 | +0.0 |

### Pooled over all held-out replicates (paired vs RAxML-NG, 1 start)

| method | n pairs | mean Δ FN | W/T/L | Wilcoxon p | mean CPU ratio vs RAxML-NG |
|---|---|---|---|---|---|
| FastTree 2 | 19 | +18.72 | 0/0/19 | 3.81e-06 | 0.06 |
| IQ-TREE 3 --fast | 19 | +9.14 | 0/0/19 | 0.000131 | 0.19 |
| RAxML-NG --fast | 19 | +1.97 | 5/2/12 | 0.00987 | 0.36 |
