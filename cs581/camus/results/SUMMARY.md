# CAMUS pilot: aggregated results

Error = (FN+FP)/2 of softwired clusters (equal to PhyloNet `CmpNets -m cluster`).
Paired by replicate; diff = variant - default (negative = better).
W/T/L = variant better / tie (|diff| < 1e-9) / worse. p = two-sided Wilcoxon signed-rank (zero differences dropped). Training condition (used to pick z settings): n25 FastTree (reps 20-39).

## Reproduction against the published CAMUS networks (V2 data, t = 0.5)

| condition | reps | our default = published (k=1 network identical) | published FN/FP (k=1) | ours FN/FP (k=1) |
|---|---|---|---|---|
| n15 iqtree | 20 | 20/20 | 0.159 / 0.112 | 0.159 / 0.112 |
| n25 fasttree | 20 | 20/20 | 0.110 / 0.082 | 0.110 / 0.082 |
| n25 iqtree | 20 | 10/20 | 0.131 / 0.116 | 0.137 / 0.110 |
| n50 fasttree | 12 | 12/12 | 0.113 / 0.079 | 0.113 / 0.079 |

## Variants vs default: 1-reticulation network (paper's protocol)

| variant | condition | n | default err | variant err | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|---|
| swap | n15 iqtree | 20 | 0.135 | 0.143 | +0.0080 | 2/15/3 | 0.438 |
| swap | n25 fasttree (train) | 20 | 0.096 | 0.098 | +0.0017 | 4/12/4 | 0.742 |
| swap | n25 iqtree | 20 | 0.124 | 0.119 | -0.0051 | 4/13/3 | 0.375 |
| swap | n50 fasttree | 8 | 0.084 | 0.076 | -0.0079 | 3/4/1 | 0.25 |
| swap | **held-out pooled** | 48 | 0.122 | 0.122 | -0.0001 | 9/32/7 | 0.597 |
| t0 | n15 iqtree | 20 | 0.135 | 0.218 | +0.0827 | 5/1/14 | 0.00484 |
| t0 | n25 fasttree (train) | 20 | 0.096 | 0.192 | +0.0960 | 0/0/20 | 8.84e-05 |
| t0 | n25 iqtree | 20 | 0.124 | 0.183 | +0.0596 | 3/2/15 | 0.000534 |
| t0 | n50 fasttree | 12 | 0.096 | 0.172 | +0.0760 | 0/0/12 | 0.000488 |
| t0 | **held-out pooled** | 52 | 0.122 | 0.194 | +0.0722 | 8/3/41 | 2.64e-07 |
| t0.2 | n15 iqtree | 20 | 0.135 | 0.150 | +0.0142 | 4/9/7 | 0.465 |
| t0.2 | n25 fasttree (train) | 20 | 0.096 | 0.095 | -0.0010 | 7/9/4 | 0.814 |
| t0.2 | n25 iqtree | 20 | 0.124 | 0.115 | -0.0091 | 9/6/5 | 0.463 |
| t0.2 | n50 fasttree | 8 | 0.084 | 0.094 | +0.0100 | 1/5/2 | 0.5 |
| t0.2 | **held-out pooled** | 48 | 0.122 | 0.126 | +0.0038 | 14/20/14 | 0.711 |
| t0.3 | n15 iqtree | 20 | 0.135 | 0.151 | +0.0152 | 2/14/4 | 0.562 |
| t0.3 | n25 fasttree (train) | 20 | 0.096 | 0.103 | +0.0070 | 3/13/4 | 0.578 |
| t0.3 | n25 iqtree | 20 | 0.124 | 0.107 | -0.0170 | 8/11/1 | 0.0547 |
| t0.3 | n50 fasttree | 8 | 0.084 | 0.092 | +0.0080 | 0/7/1 | 1 |
| t0.3 | **held-out pooled** | 48 | 0.122 | 0.123 | +0.0006 | 10/32/6 | 0.782 |
| t0.7 | n15 iqtree | 20 | 0.135 | 0.160 | +0.0246 | 2/15/3 | 0.312 |
| t0.7 | n25 fasttree (train) | 20 | 0.096 | 0.099 | +0.0025 | 0/18/2 | 0.5 |
| t0.7 | n25 iqtree | 20 | 0.124 | 0.136 | +0.0123 | 3/11/6 | 0.164 |
| t0.7 | n50 fasttree | 8 | 0.084 | 0.095 | +0.0112 | 0/5/3 | 0.25 |
| t0.7 | **held-out pooled** | 48 | 0.122 | 0.139 | +0.0172 | 5/31/12 | 0.0569 |
| t0.8 | n15 iqtree | 20 | 0.135 | 0.180 | +0.0443 | 2/11/7 | 0.0742 |
| t0.8 | n25 fasttree (train) | 20 | 0.096 | 0.103 | +0.0070 | 0/17/3 | 0.25 |
| t0.8 | n25 iqtree | 20 | 0.124 | 0.132 | +0.0079 | 4/10/6 | 0.375 |
| t0.8 | n50 fasttree | 8 | 0.084 | 0.095 | +0.0112 | 0/5/3 | 0.25 |
| t0.8 | **held-out pooled** | 48 | 0.122 | 0.146 | +0.0236 | 6/26/16 | 0.0209 |
| tqmc | n15 iqtree | 20 | 0.135 | 0.133 | -0.0025 | 3/15/2 | 1 |
| tqmc | n25 fasttree (train) | 20 | 0.096 | 0.095 | -0.0008 | 4/12/4 | 0.641 |
| tqmc | n25 iqtree | 20 | 0.124 | 0.111 | -0.0132 | 5/12/3 | 0.195 |
| tqmc | n50 fasttree | 8 | 0.084 | 0.091 | +0.0066 | 3/3/2 | 0.812 |
| tqmc | **held-out pooled** | 48 | 0.122 | 0.117 | -0.0054 | 11/30/7 | 0.417 |
| true_major | n15 iqtree | 20 | 0.135 | 0.035 | -0.1007 | 16/3/1 | 0.000598 |
| true_major | n25 fasttree (train) | 20 | 0.096 | 0.055 | -0.0408 | 17/1/2 | 0.0033 |
| true_major | n25 iqtree | 20 | 0.124 | 0.054 | -0.0696 | 16/3/1 | 0.000107 |
| true_major | n50 fasttree | 8 | 0.084 | 0.045 | -0.0393 | 7/0/1 | 0.0156 |
| true_major | **held-out pooled** | 48 | 0.122 | 0.045 | -0.0775 | 39/6/3 | 5.95e-08 |
| wastral | n15 iqtree | 20 | 0.135 | 0.189 | +0.0538 | 3/5/12 | 0.0256 |
| wastral | n25 fasttree (train) | 20 | 0.096 | 0.110 | +0.0138 | 2/12/6 | 0.461 |
| wastral | n25 iqtree | 20 | 0.124 | 0.154 | +0.0307 | 4/9/7 | 0.123 |
| wastral | n50 fasttree | 8 | 0.084 | 0.101 | +0.0171 | 5/1/2 | 0.938 |
| wastral | **held-out pooled** | 48 | 0.122 | 0.160 | +0.0380 | 12/15/21 | 0.00437 |
| z2and | n15 iqtree | 20 | 0.135 | 0.123 | -0.0121 | 2/18/0 | 0.5 |
| z2and | n25 fasttree (train) | 20 | 0.096 | 0.096 | +0.0000 | 0/20/0 | 1 |
| z2and | n25 iqtree | 20 | 0.124 | 0.124 | +0.0000 | 0/20/0 | 1 |
| z2and | n50 fasttree | 8 | 0.084 | 0.084 | +0.0000 | 0/8/0 | 1 |
| z2and | **held-out pooled** | 48 | 0.122 | 0.117 | -0.0050 | 2/46/0 | 0.5 |
| z3and | n15 iqtree | 20 | 0.135 | 0.123 | -0.0121 | 2/18/0 | 0.5 |
| z3and | n25 fasttree (train) | 20 | 0.096 | 0.096 | +0.0000 | 0/20/0 | 1 |
| z3and | n25 iqtree | 20 | 0.124 | 0.112 | -0.0113 | 4/16/0 | 0.125 |
| z3and | n50 fasttree | 8 | 0.084 | 0.084 | +0.0000 | 0/8/0 | 1 |
| z3and | **held-out pooled** | 48 | 0.122 | 0.112 | -0.0097 | 6/42/0 | 0.0312 |
| z3only | n15 iqtree | 20 | 0.135 | 0.130 | -0.0050 | 6/8/6 | 0.677 |
| z3only | n25 fasttree (train) | 20 | 0.096 | 0.087 | -0.0096 | 8/9/3 | 0.451 |
| z3only | n25 iqtree | 20 | 0.124 | 0.097 | -0.0263 | 13/5/2 | 0.0353 |
| z3only | n50 fasttree | 8 | 0.084 | 0.095 | +0.0112 | 1/4/3 | 0.25 |
| z3only | **held-out pooled** | 48 | 0.122 | 0.111 | -0.0112 | 20/17/11 | 0.239 |
| z3or | n15 iqtree | 20 | 0.135 | 0.132 | -0.0030 | 5/9/6 | 0.831 |
| z3or | n25 fasttree (train) | 20 | 0.096 | 0.087 | -0.0096 | 8/9/3 | 0.451 |
| z3or | n25 iqtree | 20 | 0.124 | 0.102 | -0.0220 | 11/7/2 | 0.0681 |
| z3or | n50 fasttree | 8 | 0.084 | 0.095 | +0.0112 | 1/4/3 | 0.25 |
| z3or | **held-out pooled** | 48 | 0.122 | 0.113 | -0.0086 | 17/20/11 | 0.399 |
| z5only | n15 iqtree | 20 | 0.135 | 0.125 | -0.0104 | 6/10/4 | 0.322 |
| z5only | n25 fasttree (train) | 20 | 0.096 | 0.088 | -0.0084 | 5/15/0 | 0.0625 |
| z5only | n25 iqtree | 20 | 0.124 | 0.095 | -0.0283 | 11/8/1 | 0.00684 |
| z5only | n50 fasttree | 8 | 0.084 | 0.091 | +0.0066 | 0/7/1 | 1 |
| z5only | **held-out pooled** | 48 | 0.122 | 0.107 | -0.0150 | 17/25/6 | 0.0301 |

## Variants vs default: network with the true number of reticulations

| variant | condition | n | default err | variant err | mean diff | W/T/L | p |
|---|---|---|---|---|---|---|---|
| swap | n15 iqtree | 20 | 0.139 | 0.144 | +0.0045 | 2/16/2 | 0.875 |
| swap | n25 fasttree (train) | 20 | 0.094 | 0.099 | +0.0049 | 4/11/5 | 0.426 |
| swap | n25 iqtree | 20 | 0.124 | 0.119 | -0.0051 | 4/13/3 | 0.375 |
| swap | n50 fasttree | 8 | 0.078 | 0.070 | -0.0079 | 3/4/1 | 0.25 |
| swap | **held-out pooled** | 48 | 0.123 | 0.121 | -0.0016 | 9/33/6 | 0.33 |
| t0 | n15 iqtree | 20 | 0.139 | 0.209 | +0.0699 | 6/1/13 | 0.01 |
| t0 | n25 fasttree (train) | 20 | 0.094 | 0.195 | +0.1010 | 0/0/20 | 8.84e-05 |
| t0 | n25 iqtree | 20 | 0.124 | 0.183 | +0.0588 | 3/2/15 | 0.00042 |
| t0 | n50 fasttree | 12 | 0.092 | 0.179 | +0.0869 | 0/0/12 | 0.000488 |
| t0 | **held-out pooled** | 52 | 0.122 | 0.192 | +0.0696 | 9/3/40 | 3.26e-07 |
| t0.2 | n15 iqtree | 20 | 0.139 | 0.153 | +0.0138 | 4/9/7 | 0.465 |
| t0.2 | n25 fasttree (train) | 20 | 0.094 | 0.093 | -0.0010 | 7/9/4 | 0.814 |
| t0.2 | n25 iqtree | 20 | 0.124 | 0.113 | -0.0108 | 10/6/4 | 0.241 |
| t0.2 | n50 fasttree | 8 | 0.078 | 0.088 | +0.0100 | 1/5/2 | 0.5 |
| t0.2 | **held-out pooled** | 48 | 0.123 | 0.126 | +0.0029 | 15/20/13 | 0.867 |
| t0.3 | n15 iqtree | 20 | 0.139 | 0.155 | +0.0154 | 2/14/4 | 0.469 |
| t0.3 | n25 fasttree (train) | 20 | 0.094 | 0.101 | +0.0070 | 3/13/4 | 0.578 |
| t0.3 | n25 iqtree | 20 | 0.124 | 0.106 | -0.0178 | 9/10/1 | 0.0371 |
| t0.3 | n50 fasttree | 8 | 0.078 | 0.086 | +0.0080 | 0/7/1 | 1 |
| t0.3 | **held-out pooled** | 48 | 0.123 | 0.123 | +0.0003 | 11/31/6 | 0.67 |
| t0.7 | n15 iqtree | 20 | 0.139 | 0.164 | +0.0246 | 2/15/3 | 0.312 |
| t0.7 | n25 fasttree (train) | 20 | 0.094 | 0.096 | +0.0025 | 0/18/2 | 0.5 |
| t0.7 | n25 iqtree | 20 | 0.124 | 0.136 | +0.0123 | 3/11/6 | 0.164 |
| t0.7 | n50 fasttree | 8 | 0.078 | 0.087 | +0.0092 | 0/6/2 | 0.5 |
| t0.7 | **held-out pooled** | 48 | 0.123 | 0.139 | +0.0169 | 5/32/11 | 0.0739 |
| t0.8 | n15 iqtree | 20 | 0.139 | 0.183 | +0.0440 | 2/11/7 | 0.0742 |
| t0.8 | n25 fasttree (train) | 20 | 0.094 | 0.101 | +0.0070 | 0/17/3 | 0.25 |
| t0.8 | n25 iqtree | 20 | 0.124 | 0.132 | +0.0079 | 4/10/6 | 0.375 |
| t0.8 | n50 fasttree | 8 | 0.078 | 0.098 | +0.0205 | 0/5/3 | 0.25 |
| t0.8 | **held-out pooled** | 48 | 0.123 | 0.148 | +0.0250 | 6/26/16 | 0.0172 |
| tqmc | n15 iqtree | 20 | 0.139 | 0.132 | -0.0076 | 4/15/1 | 0.312 |
| tqmc | n25 fasttree (train) | 20 | 0.094 | 0.091 | -0.0025 | 4/12/4 | 0.742 |
| tqmc | n25 iqtree | 20 | 0.124 | 0.111 | -0.0132 | 5/12/3 | 0.195 |
| tqmc | n50 fasttree | 8 | 0.078 | 0.083 | +0.0054 | 4/2/2 | 1 |
| tqmc | **held-out pooled** | 48 | 0.123 | 0.115 | -0.0078 | 13/29/6 | 0.156 |
| true_major | n15 iqtree | 20 | 0.139 | 0.036 | -0.1031 | 16/3/1 | 0.000502 |
| true_major | n25 fasttree (train) | 20 | 0.094 | 0.053 | -0.0410 | 17/1/2 | 0.0033 |
| true_major | n25 iqtree | 20 | 0.124 | 0.054 | -0.0697 | 17/2/1 | 7.63e-05 |
| true_major | n50 fasttree | 8 | 0.078 | 0.036 | -0.0422 | 7/0/1 | 0.0156 |
| true_major | **held-out pooled** | 48 | 0.123 | 0.044 | -0.0790 | 40/5/3 | 3.93e-08 |
| wastral | n15 iqtree | 20 | 0.139 | 0.186 | +0.0468 | 3/5/12 | 0.0413 |
| wastral | n25 fasttree (train) | 20 | 0.094 | 0.106 | +0.0121 | 2/13/5 | 0.578 |
| wastral | n25 iqtree | 20 | 0.124 | 0.154 | +0.0307 | 4/9/7 | 0.123 |
| wastral | n50 fasttree | 8 | 0.078 | 0.086 | +0.0084 | 5/1/2 | 0.938 |
| wastral | **held-out pooled** | 48 | 0.123 | 0.156 | +0.0337 | 12/15/21 | 0.00643 |
| z2and | n15 iqtree | 20 | 0.139 | 0.127 | -0.0121 | 2/18/0 | 0.5 |
| z2and | n25 fasttree (train) | 20 | 0.094 | 0.094 | +0.0000 | 0/20/0 | 1 |
| z2and | n25 iqtree | 20 | 0.124 | 0.124 | +0.0000 | 0/20/0 | 1 |
| z2and | n50 fasttree | 8 | 0.078 | 0.078 | +0.0000 | 0/8/0 | 1 |
| z2and | **held-out pooled** | 48 | 0.123 | 0.118 | -0.0050 | 2/46/0 | 0.5 |
| z3and | n15 iqtree | 20 | 0.139 | 0.127 | -0.0121 | 2/18/0 | 0.5 |
| z3and | n25 fasttree (train) | 20 | 0.094 | 0.094 | +0.0000 | 0/20/0 | 1 |
| z3and | n25 iqtree | 20 | 0.124 | 0.112 | -0.0113 | 4/16/0 | 0.125 |
| z3and | n50 fasttree | 8 | 0.078 | 0.078 | +0.0000 | 0/8/0 | 1 |
| z3and | **held-out pooled** | 48 | 0.123 | 0.113 | -0.0097 | 6/42/0 | 0.0312 |
| z3only | n15 iqtree | 20 | 0.139 | 0.132 | -0.0074 | 7/7/6 | 0.588 |
| z3only | n25 fasttree (train) | 20 | 0.094 | 0.082 | -0.0114 | 9/8/3 | 0.33 |
| z3only | n25 iqtree | 20 | 0.124 | 0.091 | -0.0330 | 14/5/1 | 0.00336 |
| z3only | n50 fasttree | 8 | 0.078 | 0.089 | +0.0112 | 1/4/3 | 0.25 |
| z3only | **held-out pooled** | 48 | 0.123 | 0.108 | -0.0150 | 22/16/10 | 0.0803 |
| z3or | n15 iqtree | 20 | 0.139 | 0.134 | -0.0054 | 6/8/6 | 0.733 |
| z3or | n25 fasttree (train) | 20 | 0.094 | 0.082 | -0.0114 | 9/8/3 | 0.33 |
| z3or | n25 iqtree | 20 | 0.124 | 0.094 | -0.0294 | 12/7/1 | 0.00806 |
| z3or | n50 fasttree | 8 | 0.078 | 0.089 | +0.0112 | 1/4/3 | 0.25 |
| z3or | **held-out pooled** | 48 | 0.123 | 0.110 | -0.0126 | 19/19/10 | 0.156 |
| z5only | n15 iqtree | 20 | 0.139 | 0.129 | -0.0102 | 6/10/4 | 0.322 |
| z5only | n25 fasttree (train) | 20 | 0.094 | 0.085 | -0.0084 | 5/15/0 | 0.0625 |
| z5only | n25 iqtree | 20 | 0.124 | 0.093 | -0.0308 | 12/8/0 | 0.000488 |
| z5only | n50 fasttree | 8 | 0.078 | 0.084 | +0.0066 | 0/7/1 | 1 |
| z5only | **held-out pooled** | 48 | 0.123 | 0.107 | -0.0160 | 18/25/5 | 0.0135 |

## Base-tree error (softwired = tree clusters vs true network clusters) and runtime

| condition | variant | base FN | base FP | k=1 FN | k=1 FP | mean CAMUS time (s, 1 thread) |
|---|---|---|---|---|---|---|
| n15 iqtree | default | 0.242 | 0.043 | 0.159 | 0.112 | 1.2 |
| n15 iqtree | swap | 0.304 | 0.118 | 0.154 | 0.133 | 1.2 |
| n15 iqtree | tqmc | 0.245 | 0.046 | 0.151 | 0.115 | 1.2 |
| n15 iqtree | true_major | 0.209 | 0.000 | 0.038 | 0.031 | 1.2 |
| n15 iqtree | wastral | 0.292 | 0.111 | 0.194 | 0.184 | 1.2 |
| n25 fasttree | default | 0.176 | 0.046 | 0.110 | 0.082 | 7.8 |
| n25 fasttree | swap | 0.200 | 0.073 | 0.100 | 0.096 | 7.8 |
| n25 fasttree | tqmc | 0.181 | 0.052 | 0.109 | 0.082 | 7.9 |
| n25 fasttree | true_major | 0.136 | 0.000 | 0.047 | 0.064 | 7.9 |
| n25 fasttree | wastral | 0.181 | 0.052 | 0.130 | 0.090 | 8.0 |
| n25 iqtree | default | 0.197 | 0.044 | 0.137 | 0.110 | 7.7 |
| n25 iqtree | swap | 0.265 | 0.119 | 0.121 | 0.116 | 7.8 |
| n25 iqtree | tqmc | 0.204 | 0.052 | 0.119 | 0.102 | 7.9 |
| n25 iqtree | true_major | 0.161 | 0.000 | 0.047 | 0.062 | 7.9 |
| n25 iqtree | wastral | 0.231 | 0.087 | 0.163 | 0.146 | 7.8 |
| n50 fasttree | default | 0.157 | 0.048 | 0.113 | 0.079 | 189.3 |
| n50 fasttree | swap | 0.158 | 0.056 | 0.085 | 0.067 | 186.1 |
| n50 fasttree | tqmc | 0.150 | 0.048 | 0.104 | 0.077 | 190.1 |
| n50 fasttree | true_major | 0.107 | 0.000 | 0.051 | 0.039 | 183.3 |
| n50 fasttree | wastral | 0.140 | 0.038 | 0.106 | 0.097 | 187.8 |

## Selecting the base tree by the CAMUS objective (no truth used)

For each replicate, among the candidate 1-reticulation networks (default, tqmc, wastral, swap), pick the one maximizing the number of filtered (t=0.5) gene-tree quartets displayed. Ties -> default.

| condition | n | default err | selected err | mean diff | W/T/L | p | oracle-best err |
|---|---|---|---|---|---|---|---|
| n15 iqtree | 20 | 0.135 | 0.134 | -0.0010 | 4/13/3 | 0.938 | 0.112 |
| n25 fasttree | 20 | 0.096 | 0.098 | +0.0013 | 5/11/4 | 0.82 | 0.079 |
| n25 iqtree | 20 | 0.124 | 0.115 | -0.0090 | 5/13/2 | 0.219 | 0.098 |
| n50 fasttree | 8 | 0.084 | 0.076 | -0.0080 | 4/3/1 | 0.125 | 0.062 |

## Threshold t: per-condition best fixed t and per-replicate oracle (upper bound for adaptive t)

| condition | t=0 | t=0.2 | t=0.3 | t=0.5 | t=0.7 | t=0.8 | per-rep oracle t |
|---|---|---|---|---|---|---|---|
| n15 iqtree | 0.218 | 0.150 | 0.151 | 0.135 | 0.160 | 0.180 | 0.112 |
| n25 fasttree | 0.192 | 0.095 | 0.103 | 0.096 | 0.099 | 0.103 | 0.081 |
| n25 iqtree | 0.183 | 0.115 | 0.107 | 0.124 | 0.136 | 0.132 | 0.094 |
| n50 fasttree | 0.165 | 0.094 | 0.092 | 0.084 | 0.095 | 0.095 | 0.083 |

