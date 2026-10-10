## Tree error (FN) on the MAGUS-paper data, mean over R0-R2

| alignment | method | 1000M2 | 1000M3 | 1000L1 | RNASim | ParkRNASim1000 | Park1000M1HF |
|---|---|---|---|---|---|---|---|
| gcm | fasttree | 10.8% (n=3) | 7.7% (n=3) | 11.2% (n=3) | 16.1% (n=3) | - | - |
| gcm | fasttree_pasta | 10.8% (n=3) | 7.5% (n=3) | 11.2% (n=3) | 16.1% (n=3) | - | - |
| gcm | iqtree_fast | 10.5% (n=3) | 8.1% (n=3) | 10.9% (n=3) | 16.5% (n=3) | - | - |
| gcm_mask50 | fasttree | 11.1% (n=3) | 7.8% (n=3) | 10.9% (n=3) | 16.6% (n=3) | - | - |
| gcm_mask70 | fasttree | 11.6% (n=3) | 7.8% (n=3) | 11.3% (n=3) | 16.5% (n=3) | - | - |
| gcm_oracle50 | fasttree | 10.7% (n=3) | 8.1% (n=3) | 11.3% (n=3) | 16.1% (n=3) | - | - |
| gcm_oracle70 | fasttree | 10.2% (n=3) | 8.0% (n=3) | 11.1% (n=3) | 17.5% (n=3) | - | - |
| gcm_pasta | fasttree | 9.8% (n=3) | 8.0% (n=3) | 11.4% (n=3) | 15.6% (n=3) | - | - |
| pasta_align | fasttree | 10.2% (n=3) | 7.2% (n=3) | 12.0% (n=3) | 15.0% (n=3) | - | - |
| pasta_align | fasttree_pasta | 10.1% (n=3) | 7.1% (n=3) | 11.6% (n=3) | 15.3% (n=3) | - | - |
| pasta_align | iqtree_fast | 10.9% (n=3) | 7.7% (n=3) | 11.8% (n=3) | 16.4% (n=3) | - | - |
| true_align | fasttree | 9.4% (n=3) | 8.0% (n=3) | 9.3% (n=3) | 15.1% (n=3) | 14.9% (n=5) | 48.9% (n=5) |
| true_align | fasttree_pasta | 9.3% (n=3) | 8.0% (n=3) | 9.8% (n=3) | 14.8% (n=3) | - | - |
| true_align | iqtree_fast | 9.5% (n=3) | 8.1% (n=3) | 10.2% (n=3) | 15.7% (n=3) | 15.6% (n=5) | 37.5% (n=5) |


## CPU seconds

| alignment | method | 1000M2 | 1000M3 | 1000L1 | RNASim | ParkRNASim1000 | Park1000M1HF |
|---|---|---|---|---|---|---|---|
| gcm | fasttree | - | - | - | - | - | - |
| gcm | fasttree_pasta | - | - | - | - | - | - |
| gcm | iqtree_fast | 150 (n=3) | 98 (n=3) | 166 (n=3) | 166 (n=3) | - | - |
| pasta_align | fasttree | - | - | - | - | - | - |
| pasta_align | fasttree_pasta | - | - | - | - | - | - |
| pasta_align | iqtree_fast | 137 (n=3) | 85 (n=3) | 166 (n=3) | 142 (n=3) | - | - |
| true_align | fasttree | - | - | - | - | 159 (n=5) | 112 (n=5) |
| true_align | fasttree_pasta | - | - | - | - | - | - |
| true_align | iqtree_fast | 149 (n=3) | 94 (n=3) | 170 (n=3) | 140 (n=3) | 137 (n=5) | 433 (n=5) |


## Validation: published Park et al. 2021 trees rescored

| method (published tree) | FN per replicate 1-5 | mean FN |
|---|---|---|
| FastTree | 17.2 14.8 14.9 14.1 13.6 | 14.9% |
| IQ-TREE 2 | 16.4 15.0 15.3 14.7 13.8 | 15.1% |
| RAxML-NG (24h cap, lastTree) | 16.9 15.2 14.7 14.6 14.2 | 15.1% |
| GTM (FastTree start, 500) | 17.1 13.5 14.9 14.4 13.5 | 14.7% |
| GTM (IQ-TREE start, 500) | 15.9 13.5 15.0 13.8 13.4 | 14.4% |
| GTM (IQ-TREE start, 120) | 16.5 13.9 14.0 14.4 12.7 | 14.3% |

1000M1-HF (Park et al. 2021 Table 5)

| method (published tree) | FN per replicate 0-4 | mean FN |
|---|---|---|
| FastTree | 53.9 47.0 50.4 49.9 55.9 | 51.4% |
| IQ-TREE 2 | 27.4 25.8 31.7 29.2 36.7 | 30.2% |
| RAxML-NG (24h cap, lastTree) | 22.9 24.0 23.5 24.4 29.8 | 24.9% |
| GTM (FastTree start, 500) | 42.5 38.5 43.0 39.0 49.2 | 42.4% |
| GTM (IQ-TREE start, 500) | 26.2 25.3 28.3 28.2 33.8 | 28.4% |


## Validation: PASTA's published tree vs our FastTree on the PASTA alignment

| dataset | rep | published pasta.tre FN | our FastTree -fastest FN | our FastTree (default) FN | RF(pasta.tre, ours -fastest) |
|---|---|---|---|---|---|
| 1000M2 | 0 | 11.9% | 11.3% | 11.3% | 1.5% |
| 1000M2 | 1 | 12.6% | 10.1% | 10.4% | 5.7% |
| 1000M2 | 2 | 9.1% | 9.1% | 9.0% | 0.7% |
| 1000M3 | 0 | 7.0% | 7.5% | 7.5% | 0.5% |
| 1000M3 | 1 | 5.8% | 5.5% | 5.7% | 1.0% |
| 1000M3 | 2 | 8.9% | 8.5% | 8.5% | 1.8% |
| 1000L1 | 0 | 10.9% | 11.0% | 10.9% | 2.4% |
| 1000L1 | 1 | 15.4% | 13.5% | 13.6% | 6.5% |
| 1000L1 | 2 | 12.3% | 10.5% | 11.5% | 7.7% |
| RNASim | 0 | 13.9% | 13.5% | 13.0% | 2.6% |
| RNASim | 1 | 14.6% | 15.0% | 14.7% | 3.8% |
| RNASim | 2 | 17.1% | 17.3% | 17.3% | 2.9% |


### Pilot: 1000M1-HF (true alignment, half the sequences fragmentary), paired vs RAxML-NG (1 parsimony start)

| method | n | mean FN | FN vs raxmlng (mean paired diff, wins/ties/losses) | mean CPU min | mean lnL - raxmlng lnL |
|---|---|---|---|---|---|
| fasttree | 5 | 48.9% | +24.5 pts (0/0/3) | 1.9 | - |
| iqtree_fast | 5 | 37.5% | +13.7 pts (0/0/3) | 7.2 | - |
| raxmlng | 3 | 23.6% | - | 35.9 | - |
| raxmlng_from_iqtree_fast | 2 | 22.0% | -0.5 pts (1/0/1) | 48.0 | +0.3 |
| frag_constr_fasttree | 2 | 22.4% | -0.1 pts (1/0/1) | 12.2 | -33.5 |
| frag_polish_fasttree | 2 | 22.2% | -0.4 pts (1/0/1) | 49.9 | +1.4 |

### Pilot: start trees, 1000M2 R0, MAGUS alignment

| method | n | mean FN | FN vs raxmlng (mean paired diff, wins/ties/losses) | mean CPU min | mean lnL - raxmlng lnL |
|---|---|---|---|---|---|
| fasttree | 1 | 10.7% | +0.8 pts (0/0/1) | 1.8 | - |
| iqtree_fast | 1 | 9.7% | -0.2 pts (1/0/0) | 2.6 | - |
| raxmlng | 1 | 9.9% | - | 45.1 | - |
| raxmlng_ft | 1 | 9.5% | -0.4 pts (1/0/0) | 49.1 | -0.3 |

### Pilot: start trees, 1000M3 R0, MAGUS alignment

| method | n | mean FN | FN vs raxmlng (mean paired diff, wins/ties/losses) | mean CPU min | mean lnL - raxmlng lnL |
|---|---|---|---|---|---|
| fasttree | 1 | 7.7% | -0.3 pts (1/0/0) | 1.2 | - |
| iqtree_fast | 1 | 7.7% | -0.3 pts (1/0/0) | 1.4 | - |
| raxmlng | 1 | 8.0% | - | 35.0 | - |
| raxmlng_ft | 1 | 7.2% | -0.8 pts (1/0/0) | 25.9 | -0.7 |

### Pilot: start trees, 1000L1 R0, MAGUS alignment

| method | n | mean FN | FN vs raxmlng (mean paired diff, wins/ties/losses) | mean CPU min | mean lnL - raxmlng lnL |
|---|---|---|---|---|---|
| fasttree | 1 | 11.7% | -0.5 pts (1/0/0) | 1.5 | - |
| iqtree_fast | 1 | 11.7% | -0.5 pts (1/0/0) | 2.5 | - |
| raxmlng | 1 | 12.2% | - | 38.9 | - |
| raxmlng_ft | 1 | 11.3% | -0.9 pts (1/0/0) | 39.8 | +2.2 |

### Pilot: start trees, RNASim R0, MAGUS alignment

| method | n | mean FN | FN vs raxmlng (mean paired diff, wins/ties/losses) | mean CPU min | mean lnL - raxmlng lnL |
|---|---|---|---|---|---|
| fasttree | 1 | 14.4% | - | 6.9 | - |
| iqtree_fast | 1 | 15.3% | - | 3.1 | - |

### Pilot: column masking / alignment ensembles (FastTree), FN change vs the unmasked MAGUS alignment

| alignment given to FastTree | 1000M2 | 1000M3 | 1000L1 | RNASim | all 12 (wins/ties/losses) |
|---|---|---|---|---|---|
| gcm_mask50 | +0.3 | +0.1 | -0.3 | +0.5 | +0.2 pts (4/0/8) |
| gcm_mask70 | +0.8 | +0.1 | +0.2 | +0.4 | +0.4 pts (3/1/8) |
| gcm_oracle50 | -0.0 | +0.4 | +0.1 | -0.0 | +0.1 pts (4/1/7) |
| gcm_oracle70 | -0.5 | +0.3 | -0.1 | +1.4 | +0.3 pts (5/0/7) |
| gcm_pasta | -0.9 | +0.3 | +0.3 | -0.5 | -0.2 pts (7/1/4) |
| pasta_align | -0.6 | -0.5 | +0.8 | -1.1 | -0.3 pts (9/0/3) |
| true_align | -1.4 | +0.3 | -1.8 | -1.0 | -1.0 pts (9/0/3) |
