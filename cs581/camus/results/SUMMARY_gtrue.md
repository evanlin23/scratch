# CAMUS pilot: aggregated results

Error = (FN+FP)/2 of softwired clusters (equal to PhyloNet `CmpNets -m cluster`).
Paired by replicate; diff = variant - default (negative = better).
W/T/L = variant better / tie (|diff| < 1e-9) / worse. p = two-sided Wilcoxon signed-rank (zero differences dropped). Training condition (used to pick z settings): n25 FastTree (reps 20-39).

## Reproduction against the published CAMUS networks (V2 data, t = 0.5)

| condition | reps | our default = published (k=1 network identical) | published FN/FP (k=1) | ours FN/FP (k=1) |
|---|---|---|---|---|

## Variants vs default: 1-reticulation network (paper's protocol)

| variant | condition | n | default err | variant err | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|---|
| swap | n25 true-gt | 20 | 0.103 | 0.102 | -0.0012 | 3/15/2 | 0.812 |
| swap | **held-out pooled** | 20 | 0.103 | 0.102 | -0.0012 | 3/15/2 | 0.812 |
| t0.2 | n25 true-gt | 20 | 0.103 | 0.088 | -0.0146 | 8/10/2 | 0.105 |
| t0.2 | **held-out pooled** | 20 | 0.103 | 0.088 | -0.0146 | 8/10/2 | 0.105 |
| tqmc | n25 true-gt | 20 | 0.103 | 0.105 | +0.0017 | 4/12/4 | 0.844 |
| tqmc | **held-out pooled** | 20 | 0.103 | 0.105 | +0.0017 | 4/12/4 | 0.844 |
| true_major | n25 true-gt | 20 | 0.103 | 0.056 | -0.0469 | 15/3/2 | 0.000381 |
| true_major | **held-out pooled** | 20 | 0.103 | 0.056 | -0.0469 | 15/3/2 | 0.000381 |
| z3only | n25 true-gt | 20 | 0.103 | 0.083 | -0.0204 | 9/10/1 | 0.0137 |
| z3only | **held-out pooled** | 20 | 0.103 | 0.083 | -0.0204 | 9/10/1 | 0.0137 |
| z5only | n25 true-gt | 20 | 0.103 | 0.087 | -0.0164 | 7/13/0 | 0.0156 |
| z5only | **held-out pooled** | 20 | 0.103 | 0.087 | -0.0164 | 7/13/0 | 0.0156 |

## Variants vs default: network with the true number of reticulations

| variant | condition | n | default err | variant err | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|---|
| swap | n25 true-gt | 20 | 0.100 | 0.102 | +0.0020 | 3/14/3 | 0.844 |
| swap | **held-out pooled** | 20 | 0.100 | 0.102 | +0.0020 | 3/14/3 | 0.844 |
| t0.2 | n25 true-gt | 20 | 0.100 | 0.086 | -0.0146 | 8/10/2 | 0.105 |
| t0.2 | **held-out pooled** | 20 | 0.100 | 0.086 | -0.0146 | 8/10/2 | 0.105 |
| tqmc | n25 true-gt | 20 | 0.100 | 0.100 | -0.0000 | 4/12/4 | 0.945 |
| tqmc | **held-out pooled** | 20 | 0.100 | 0.100 | -0.0000 | 4/12/4 | 0.945 |
| true_major | n25 true-gt | 20 | 0.100 | 0.053 | -0.0471 | 15/3/2 | 0.000381 |
| true_major | **held-out pooled** | 20 | 0.100 | 0.053 | -0.0471 | 15/3/2 | 0.000381 |
| z3only | n25 true-gt | 20 | 0.100 | 0.078 | -0.0222 | 10/9/1 | 0.00684 |
| z3only | **held-out pooled** | 20 | 0.100 | 0.078 | -0.0222 | 10/9/1 | 0.00684 |
| z5only | n25 true-gt | 20 | 0.100 | 0.084 | -0.0164 | 7/13/0 | 0.0156 |
| z5only | **held-out pooled** | 20 | 0.100 | 0.084 | -0.0164 | 7/13/0 | 0.0156 |

## Base-tree error (softwired = tree clusters vs true network clusters) and runtime

| condition | variant | base FN | base FP | k=1 FN | k=1 FP | mean CAMUS time (s, 1 thread) |
|---|---|---|---|---|---|---|
| n25 true-gt | default | 0.166 | 0.035 | 0.108 | 0.098 | 15.1 |
| n25 true-gt | swap | 0.221 | 0.098 | 0.098 | 0.106 | 15.5 |
| n25 true-gt | tqmc | 0.175 | 0.046 | 0.110 | 0.099 | 15.7 |
| n25 true-gt | true_major | 0.136 | 0.000 | 0.047 | 0.065 | 15.5 |

## Selecting the base tree by the CAMUS objective (no truth used)

For each replicate, among the candidate 1-reticulation networks (default, tqmc, wastral, swap), pick the one maximizing the number of filtered (t=0.5) gene-tree quartets displayed. Ties -> default.

| condition | n | default err | selected err | mean diff | W/T/L | p | oracle-best err |
|---|---|---|---|---|---|---|---|

## Threshold t: per-condition best fixed t and per-replicate oracle (upper bound for adaptive t)

| condition | t=0 | t=0.2 | t=0.3 | t=0.5 | t=0.7 | t=0.8 | per-rep oracle t |
|---|---|---|---|---|---|---|---|

