# CAMUS pilot: aggregated results

Error = (FN+FP)/2 of softwired clusters (equal to PhyloNet `CmpNets -m cluster`).
Paired by replicate; diff = variant - default (negative = better).
W/T/L = variant better / tie (|diff| < 1e-9) / worse. p = two-sided Wilcoxon signed-rank (zero differences dropped). Training condition (used to pick z settings): n25 FastTree (reps 20-39).

## Reproduction against the published CAMUS networks (V2 data, t = 0.5)

| condition | reps | our default = published (k=1 network identical) | published FN/FP (k=1) | ours FN/FP (k=1) |
|---|---|---|---|---|
| n15 iqtree | 7 | 7/7 | 0.159 / 0.094 | 0.159 / 0.094 |
| n25 fasttree | 20 | 20/20 | 0.110 / 0.082 | 0.110 / 0.082 |

## Variants vs default: 1-reticulation network (paper's protocol)

| variant | condition | n | default err | variant err | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|---|
| swap | n15 iqtree | 3 | 0.113 | 0.113 | +0.0000 | 0/3/0 | 1 |
| swap | n25 fasttree (train) | 20 | 0.096 | 0.098 | +0.0017 | 4/12/4 | 0.742 |
| swap | **held-out pooled** | 3 | 0.113 | 0.113 | +0.0000 | 0/3/0 | 1 |
| t0 | n15 iqtree | 7 | 0.126 | 0.213 | +0.0865 | 2/1/4 | 0.156 |
| t0 | n25 fasttree (train) | 20 | 0.096 | 0.192 | +0.0960 | 0/0/20 | 8.84e-05 |
| t0 | **held-out pooled** | 7 | 0.126 | 0.213 | +0.0865 | 2/1/4 | 0.156 |
| t0.2 | n15 iqtree | 6 | 0.121 | 0.191 | +0.0708 | 1/2/3 | 0.25 |
| t0.2 | n25 fasttree (train) | 20 | 0.096 | 0.095 | -0.0010 | 7/9/4 | 0.814 |
| t0.2 | **held-out pooled** | 6 | 0.121 | 0.191 | +0.0708 | 1/2/3 | 0.25 |
| t0.3 | n15 iqtree | 6 | 0.121 | 0.169 | +0.0486 | 1/3/2 | 0.5 |
| t0.3 | n25 fasttree (train) | 20 | 0.096 | 0.103 | +0.0070 | 3/13/4 | 0.578 |
| t0.3 | **held-out pooled** | 6 | 0.121 | 0.169 | +0.0486 | 1/3/2 | 0.5 |
| t0.7 | n15 iqtree | 6 | 0.121 | 0.139 | +0.0184 | 0/5/1 | 1 |
| t0.7 | n25 fasttree (train) | 20 | 0.096 | 0.099 | +0.0025 | 0/18/2 | 0.5 |
| t0.7 | **held-out pooled** | 6 | 0.121 | 0.139 | +0.0184 | 0/5/1 | 1 |
| t0.8 | n15 iqtree | 6 | 0.121 | 0.167 | +0.0465 | 0/3/3 | 0.25 |
| t0.8 | n25 fasttree (train) | 20 | 0.096 | 0.103 | +0.0070 | 0/17/3 | 0.25 |
| t0.8 | **held-out pooled** | 6 | 0.121 | 0.167 | +0.0465 | 0/3/3 | 0.25 |
| tqmc | n15 iqtree | 3 | 0.113 | 0.113 | +0.0000 | 0/3/0 | 1 |
| tqmc | n25 fasttree (train) | 20 | 0.096 | 0.095 | -0.0008 | 4/12/4 | 0.641 |
| tqmc | **held-out pooled** | 3 | 0.113 | 0.113 | +0.0000 | 0/3/0 | 1 |
| true_major | n15 iqtree | 3 | 0.113 | 0.034 | -0.0789 | 2/1/0 | 0.5 |
| true_major | n25 fasttree (train) | 20 | 0.096 | 0.055 | -0.0408 | 17/1/2 | 0.0033 |
| true_major | **held-out pooled** | 3 | 0.113 | 0.034 | -0.0789 | 2/1/0 | 0.5 |
| wastral | n15 iqtree | 3 | 0.113 | 0.187 | +0.0742 | 0/0/3 | 0.25 |
| wastral | n25 fasttree (train) | 20 | 0.096 | 0.110 | +0.0138 | 2/12/6 | 0.461 |
| wastral | **held-out pooled** | 3 | 0.113 | 0.187 | +0.0742 | 0/0/3 | 0.25 |
| z2and | n15 iqtree | 5 | 0.138 | 0.138 | +0.0000 | 0/5/0 | 1 |
| z2and | n25 fasttree (train) | 20 | 0.096 | 0.096 | +0.0000 | 0/20/0 | 1 |
| z2and | **held-out pooled** | 5 | 0.138 | 0.138 | +0.0000 | 0/5/0 | 1 |
| z3and | n15 iqtree | 5 | 0.138 | 0.138 | +0.0000 | 0/5/0 | 1 |
| z3and | n25 fasttree (train) | 20 | 0.096 | 0.096 | +0.0000 | 0/20/0 | 1 |
| z3and | **held-out pooled** | 5 | 0.138 | 0.138 | +0.0000 | 0/5/0 | 1 |
| z3only | n15 iqtree | 4 | 0.159 | 0.216 | +0.0574 | 1/1/2 | 0.5 |
| z3only | n25 fasttree (train) | 20 | 0.096 | 0.087 | -0.0096 | 8/9/3 | 0.451 |
| z3only | **held-out pooled** | 4 | 0.159 | 0.216 | +0.0574 | 1/1/2 | 0.5 |
| z3or | n15 iqtree | 4 | 0.159 | 0.216 | +0.0574 | 1/1/2 | 0.5 |
| z3or | n25 fasttree (train) | 20 | 0.096 | 0.087 | -0.0096 | 8/9/3 | 0.451 |
| z3or | **held-out pooled** | 4 | 0.159 | 0.216 | +0.0574 | 1/1/2 | 0.5 |
| z5only | n15 iqtree | 4 | 0.159 | 0.216 | +0.0574 | 1/1/2 | 0.5 |
| z5only | n25 fasttree (train) | 20 | 0.096 | 0.088 | -0.0084 | 5/15/0 | 0.0625 |
| z5only | **held-out pooled** | 4 | 0.159 | 0.216 | +0.0574 | 1/1/2 | 0.5 |

## Variants vs default: network with the true number of reticulations

| variant | condition | n | default err | variant err | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|---|
| swap | n15 iqtree | 3 | 0.113 | 0.113 | +0.0000 | 0/3/0 | 1 |
| swap | n25 fasttree (train) | 20 | 0.094 | 0.099 | +0.0049 | 4/11/5 | 0.426 |
| swap | **held-out pooled** | 3 | 0.113 | 0.113 | +0.0000 | 0/3/0 | 1 |
| t0 | n15 iqtree | 7 | 0.126 | 0.213 | +0.0865 | 2/1/4 | 0.156 |
| t0 | n25 fasttree (train) | 20 | 0.094 | 0.195 | +0.1010 | 0/0/20 | 8.84e-05 |
| t0 | **held-out pooled** | 7 | 0.126 | 0.213 | +0.0865 | 2/1/4 | 0.156 |
| t0.2 | n15 iqtree | 6 | 0.121 | 0.191 | +0.0708 | 1/2/3 | 0.25 |
| t0.2 | n25 fasttree (train) | 20 | 0.094 | 0.093 | -0.0010 | 7/9/4 | 0.814 |
| t0.2 | **held-out pooled** | 6 | 0.121 | 0.191 | +0.0708 | 1/2/3 | 0.25 |
| t0.3 | n15 iqtree | 6 | 0.121 | 0.169 | +0.0486 | 1/3/2 | 0.5 |
| t0.3 | n25 fasttree (train) | 20 | 0.094 | 0.101 | +0.0070 | 3/13/4 | 0.578 |
| t0.3 | **held-out pooled** | 6 | 0.121 | 0.169 | +0.0486 | 1/3/2 | 0.5 |
| t0.7 | n15 iqtree | 6 | 0.121 | 0.139 | +0.0184 | 0/5/1 | 1 |
| t0.7 | n25 fasttree (train) | 20 | 0.094 | 0.096 | +0.0025 | 0/18/2 | 0.5 |
| t0.7 | **held-out pooled** | 6 | 0.121 | 0.139 | +0.0184 | 0/5/1 | 1 |
| t0.8 | n15 iqtree | 6 | 0.121 | 0.167 | +0.0465 | 0/3/3 | 0.25 |
| t0.8 | n25 fasttree (train) | 20 | 0.094 | 0.101 | +0.0070 | 0/17/3 | 0.25 |
| t0.8 | **held-out pooled** | 6 | 0.121 | 0.167 | +0.0465 | 0/3/3 | 0.25 |
| tqmc | n15 iqtree | 3 | 0.113 | 0.113 | +0.0000 | 0/3/0 | 1 |
| tqmc | n25 fasttree (train) | 20 | 0.094 | 0.091 | -0.0025 | 4/12/4 | 0.742 |
| tqmc | **held-out pooled** | 3 | 0.113 | 0.113 | +0.0000 | 0/3/0 | 1 |
| true_major | n15 iqtree | 3 | 0.113 | 0.034 | -0.0789 | 2/1/0 | 0.5 |
| true_major | n25 fasttree (train) | 20 | 0.094 | 0.053 | -0.0410 | 17/1/2 | 0.0033 |
| true_major | **held-out pooled** | 3 | 0.113 | 0.034 | -0.0789 | 2/1/0 | 0.5 |
| wastral | n15 iqtree | 3 | 0.113 | 0.187 | +0.0742 | 0/0/3 | 0.25 |
| wastral | n25 fasttree (train) | 20 | 0.094 | 0.106 | +0.0121 | 2/13/5 | 0.578 |
| wastral | **held-out pooled** | 3 | 0.113 | 0.187 | +0.0742 | 0/0/3 | 0.25 |
| z2and | n15 iqtree | 5 | 0.138 | 0.138 | +0.0000 | 0/5/0 | 1 |
| z2and | n25 fasttree (train) | 20 | 0.094 | 0.094 | +0.0000 | 0/20/0 | 1 |
| z2and | **held-out pooled** | 5 | 0.138 | 0.138 | +0.0000 | 0/5/0 | 1 |
| z3and | n15 iqtree | 5 | 0.138 | 0.138 | +0.0000 | 0/5/0 | 1 |
| z3and | n25 fasttree (train) | 20 | 0.094 | 0.094 | +0.0000 | 0/20/0 | 1 |
| z3and | **held-out pooled** | 5 | 0.138 | 0.138 | +0.0000 | 0/5/0 | 1 |
| z3only | n15 iqtree | 4 | 0.159 | 0.216 | +0.0574 | 1/1/2 | 0.5 |
| z3only | n25 fasttree (train) | 20 | 0.094 | 0.082 | -0.0114 | 9/8/3 | 0.33 |
| z3only | **held-out pooled** | 4 | 0.159 | 0.216 | +0.0574 | 1/1/2 | 0.5 |
| z3or | n15 iqtree | 4 | 0.159 | 0.216 | +0.0574 | 1/1/2 | 0.5 |
| z3or | n25 fasttree (train) | 20 | 0.094 | 0.082 | -0.0114 | 9/8/3 | 0.33 |
| z3or | **held-out pooled** | 4 | 0.159 | 0.216 | +0.0574 | 1/1/2 | 0.5 |
| z5only | n15 iqtree | 4 | 0.159 | 0.216 | +0.0574 | 1/1/2 | 0.5 |
| z5only | n25 fasttree (train) | 20 | 0.094 | 0.085 | -0.0084 | 5/15/0 | 0.0625 |
| z5only | **held-out pooled** | 4 | 0.159 | 0.216 | +0.0574 | 1/1/2 | 0.5 |

## Base-tree error (softwired = tree clusters vs true network clusters) and runtime

| condition | variant | base FN | base FP | k=1 FN | k=1 FP | mean CAMUS time (s, 1 thread) |
|---|---|---|---|---|---|---|
| n15 iqtree | default | 0.227 | 0.041 | 0.159 | 0.094 | 1.1 |
| n15 iqtree | swap | 0.273 | 0.000 | 0.184 | 0.042 | 1.1 |
| n15 iqtree | tqmc | 0.306 | 0.048 | 0.184 | 0.042 | 1.2 |
| n15 iqtree | true_major | 0.273 | 0.000 | 0.068 | 0.000 | 1.1 |
| n15 iqtree | wastral | 0.374 | 0.143 | 0.237 | 0.137 | 1.1 |
| n25 fasttree | default | 0.176 | 0.046 | 0.110 | 0.082 | 7.8 |
| n25 fasttree | swap | 0.200 | 0.073 | 0.100 | 0.096 | 7.8 |
| n25 fasttree | tqmc | 0.181 | 0.052 | 0.109 | 0.082 | 7.9 |
| n25 fasttree | true_major | 0.136 | 0.000 | 0.047 | 0.064 | 7.9 |
| n25 fasttree | wastral | 0.181 | 0.052 | 0.130 | 0.090 | 8.0 |

## Selecting the base tree by the CAMUS objective (no truth used)

For each replicate, among the candidate 1-reticulation networks (default, tqmc, wastral, swap), pick the one maximizing the number of filtered (t=0.5) gene-tree quartets displayed. Ties -> default.

| condition | n | default err | selected err | mean diff | W/T/L | p | oracle-best err |
|---|---|---|---|---|---|---|---|
| n15 iqtree | 3 | 0.113 | 0.113 | +0.0000 | 0/3/0 | 1 | 0.113 |
| n25 fasttree | 20 | 0.096 | 0.098 | +0.0013 | 5/11/4 | 0.82 | 0.079 |

## Threshold t: per-condition best fixed t and per-replicate oracle (upper bound for adaptive t)

| condition | t=0 | t=0.2 | t=0.3 | t=0.5 | t=0.7 | t=0.8 | per-rep oracle t |
|---|---|---|---|---|---|---|---|
| n15 iqtree | 0.224 | 0.191 | 0.169 | 0.121 | 0.139 | 0.167 | 0.116 |
| n25 fasttree | 0.192 | 0.095 | 0.103 | 0.096 | 0.099 | 0.103 | 0.081 |

