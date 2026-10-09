# CAMUS pilot: aggregated results

Error = (FN+FP)/2 of softwired clusters (equal to PhyloNet `CmpNets -m cluster`).
Paired by replicate; diff = variant - default (negative = better).
W/T/L = variant better / tie (|diff| < 1e-9) / worse. p = two-sided Wilcoxon signed-rank (zero differences dropped). Training condition (used to pick z settings): n25 FastTree (reps 20-39).

## Reproduction against the published CAMUS networks (V2 data, t = 0.5)

| condition | reps | our default = published (k=1 network identical) | published FN/FP (k=1) | ours FN/FP (k=1) |
|---|---|---|---|---|
| n25 fasttree | 4 | 4/4 | 0.101 / 0.058 | 0.101 / 0.058 |

## Variants vs default: 1-reticulation network (paper's protocol)

| variant | condition | n | default err | variant err | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|---|
| t0 | n25 fasttree (train) | 4 | 0.079 | 0.161 | +0.0816 | 0/0/4 | 0.125 |
| t0.2 | n25 fasttree (train) | 4 | 0.079 | 0.091 | +0.0121 | 2/1/1 | 1 |
| t0.3 | n25 fasttree (train) | 4 | 0.079 | 0.103 | +0.0235 | 0/3/1 | 1 |

## Variants vs default: network with the true number of reticulations

| variant | condition | n | default err | variant err | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|---|
| t0 | n25 fasttree (train) | 4 | 0.079 | 0.161 | +0.0816 | 0/0/4 | 0.125 |
| t0.2 | n25 fasttree (train) | 4 | 0.079 | 0.091 | +0.0121 | 2/1/1 | 1 |
| t0.3 | n25 fasttree (train) | 4 | 0.079 | 0.103 | +0.0235 | 0/3/1 | 1 |

## Base-tree error (softwired = tree clusters vs true network clusters) and runtime

| condition | variant | base FN | base FP | k=1 FN | k=1 FP | mean CAMUS time (s, 1 thread) |
|---|---|---|---|---|---|---|
| n25 fasttree | default | 0.154 | 0.031 | 0.101 | 0.058 | 8.5 |

## Selecting the base tree by the CAMUS objective (no truth used)

For each replicate, among the candidate 1-reticulation networks (default, tqmc, wastral, swap), pick the one maximizing the number of filtered (t=0.5) gene-tree quartets displayed. Ties -> default.

| condition | n | default err | selected err | mean diff | W/T/L | p | oracle-best err |
|---|---|---|---|---|---|---|---|

## Threshold t: per-condition best fixed t and per-replicate oracle (upper bound for adaptive t)

| condition | t=0 | t=0.2 | t=0.3 | t=0.5 | t=0.7 | t=0.8 | per-rep oracle t |
|---|---|---|---|---|---|---|---|

