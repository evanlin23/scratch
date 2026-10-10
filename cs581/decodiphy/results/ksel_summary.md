runs: 1338 (train 750, test 588); seeds <= 1000000000; runs where the authors' solver crashed (excluded): 12
tuned on train trees (['1kp', 'bees', 'beetles', 'birds-jarvis', 'hemipteroid']): bic lam=2, ratio tau=0.9, paper min_p=0.005, adj min_p=0.005
logistic stop-rule coefficients [log10 min p, rel. ybar drop, log10 loss ratio, adjacent]: [-1.89, -8.82, 0.39, 1.93] intercept -2.99

### TEST trees, noise0 (n=196)
| rule | exact-k acc | mean abs(k_hat-k) | bias | mean Jaccard | dJac vs paper | Wilcoxon p | wins/losses |
|---|---|---|---|---|---|---|---|
| paper(min_p=0.01) | 0.735 | 0.36 | -0.34 | 0.918 | +0.000 | nan | 0/0 |
| bic(lam=2) | 0.668 | 0.56 | +0.55 | 0.882 | -0.036 | 0.0067 | 37/48 |
| ratio(tau=0.9) | 0.964 | 0.04 | -0.01 | 0.962 | +0.044 | 4.1e-08 | 45/4 |
| paper(min_p tuned=0.005) | 0.883 | 0.12 | -0.10 | 0.951 | +0.033 | 6.3e-06 | 29/3 |
| adj (paper+adjacency stop) | 0.714 | 0.41 | -0.40 | 0.912 | -0.006 | 0.074 | 2/7 |
| adj(min_p tuned=0.005) | 0.847 | 0.19 | -0.18 | 0.945 | +0.027 | 0.0008 | 30/9 |
| learned (logistic, 4 features) | 0.893 | 0.20 | -0.18 | 0.953 | +0.035 | 1.6e-06 | 37/4 |
| oracle(true k) | 1.000 | 0.00 | +0.00 | 0.958 | +0.040 | 2.2e-07 | 44/8 |

### TEST trees, noise1 (n=196)
| rule | exact-k acc | mean abs(k_hat-k) | bias | mean Jaccard | dJac vs paper | Wilcoxon p | wins/losses |
|---|---|---|---|---|---|---|---|
| paper(min_p=0.01) | 0.153 | 1.73 | +1.36 | 0.454 | +0.000 | nan | 0/0 |
| bic(lam=2) | 0.005 | 2.83 | +2.83 | 0.419 | -0.035 | 0.00093 | 40/93 |
| ratio(tau=0.9) | 0.000 | 2.98 | +2.98 | 0.409 | -0.045 | 5.6e-05 | 38/98 |
| paper(min_p tuned=0.005) | 0.077 | 2.19 | +2.11 | 0.442 | -0.012 | 0.14 | 28/41 |
| adj (paper+adjacency stop) | 0.173 | 1.62 | +1.09 | 0.479 | +0.024 | 2.3e-05 | 25/3 |
| adj(min_p tuned=0.005) | 0.117 | 1.98 | +1.73 | 0.475 | +0.021 | 0.017 | 54/37 |
| learned (logistic, 4 features) | 0.270 | 1.19 | +0.47 | 0.489 | +0.035 | 0.00067 | 69/30 |
| oracle(true k) | 1.000 | 0.00 | +0.00 | 0.503 | +0.048 | 0.11 | 83/79 |

### TEST trees, noise2 (n=196)
| rule | exact-k acc | mean abs(k_hat-k) | bias | mean Jaccard | dJac vs paper | Wilcoxon p | wins/losses |
|---|---|---|---|---|---|---|---|
| paper(min_p=0.01) | 0.082 | 2.05 | +1.40 | 0.307 | +0.000 | nan | 0/0 |
| bic(lam=2) | 0.005 | 2.78 | +2.70 | 0.289 | -0.018 | 0.088 | 44/66 |
| ratio(tau=0.9) | 0.000 | 3.00 | +3.00 | 0.277 | -0.031 | 0.00034 | 33/76 |
| paper(min_p tuned=0.005) | 0.056 | 2.39 | +2.13 | 0.304 | -0.003 | 0.79 | 24/32 |
| adj (paper+adjacency stop) | 0.082 | 2.05 | +1.38 | 0.309 | +0.002 | nan | 2/0 |
| adj(min_p tuned=0.005) | 0.061 | 2.36 | +2.09 | 0.308 | +0.001 | 0.76 | 26/32 |
| learned (logistic, 4 features) | 0.112 | 1.79 | +0.79 | 0.305 | -0.003 | 0.33 | 25/27 |
| oracle(true k) | 1.000 | 0.00 | +0.00 | 0.256 | -0.051 | 1.3e-07 | 45/121 |

### TEST trees, all (n=588)
| rule | exact-k acc | mean abs(k_hat-k) | bias | mean Jaccard | dJac vs paper | Wilcoxon p | wins/losses |
|---|---|---|---|---|---|---|---|
| paper(min_p=0.01) | 0.323 | 1.38 | +0.81 | 0.560 | +0.000 | nan | 0/0 |
| bic(lam=2) | 0.226 | 2.05 | +2.03 | 0.530 | -0.030 | 5.6e-06 | 121/207 |
| ratio(tau=0.9) | 0.321 | 2.01 | +1.99 | 0.549 | -0.010 | 0.066 | 116/178 |
| paper(min_p tuned=0.005) | 0.338 | 1.57 | +1.38 | 0.566 | +0.006 | 0.11 | 81/76 |
| adj (paper+adjacency stop) | 0.323 | 1.36 | +0.69 | 0.567 | +0.007 | 0.0014 | 29/10 |
| adj(min_p tuned=0.005) | 0.342 | 1.51 | +1.21 | 0.576 | +0.016 | 0.00014 | 110/78 |
| learned (logistic, 4 features) | 0.425 | 1.06 | +0.36 | 0.582 | +0.022 | 2e-06 | 131/61 |
| oracle(true k) | 1.000 | 0.00 | +0.00 | 0.572 | +0.012 | 0.89 | 172/208 |

### ALL trees, noise0 (n=446)
| rule | exact-k acc | mean abs(k_hat-k) | bias | mean Jaccard | dJac vs paper | Wilcoxon p | wins/losses |
|---|---|---|---|---|---|---|---|
| paper(min_p=0.01) | 0.747 | 0.34 | -0.33 | 0.922 | +0.000 | nan | 0/0 |
| bic(lam=2) | 0.713 | 0.50 | +0.48 | 0.892 | -0.030 | 0.00026 | 80/105 |
| ratio(tau=0.9) | 0.962 | 0.04 | -0.02 | 0.961 | +0.039 | 7.9e-15 | 93/13 |
| paper(min_p tuned=0.005) | 0.874 | 0.15 | -0.14 | 0.950 | +0.028 | 7.6e-11 | 58/3 |
| adj (paper+adjacency stop) | 0.709 | 0.39 | -0.39 | 0.912 | -0.010 | 0.00041 | 3/20 |
| adj(min_p tuned=0.005) | 0.825 | 0.22 | -0.21 | 0.940 | +0.018 | 0.00044 | 58/22 |
| learned (logistic, 4 features) | 0.868 | 0.29 | -0.27 | 0.933 | +0.011 | 0.001 | 65/21 |
| oracle(true k) | 1.000 | 0.00 | +0.00 | 0.957 | +0.035 | 1.3e-12 | 91/22 |

### ALL trees, noise1 (n=446)
| rule | exact-k acc | mean abs(k_hat-k) | bias | mean Jaccard | dJac vs paper | Wilcoxon p | wins/losses |
|---|---|---|---|---|---|---|---|
| paper(min_p=0.01) | 0.143 | 1.70 | +1.29 | 0.520 | +0.000 | nan | 0/0 |
| bic(lam=2) | 0.011 | 2.82 | +2.78 | 0.453 | -0.067 | 7.8e-20 | 64/245 |
| ratio(tau=0.9) | 0.007 | 2.92 | +2.90 | 0.449 | -0.072 | 2.6e-22 | 59/256 |
| paper(min_p tuned=0.005) | 0.087 | 2.11 | +2.01 | 0.502 | -0.018 | 6.7e-05 | 46/109 |
| adj (paper+adjacency stop) | 0.191 | 1.51 | +0.88 | 0.550 | +0.029 | 7.7e-10 | 71/16 |
| adj(min_p tuned=0.005) | 0.152 | 1.83 | +1.43 | 0.539 | +0.018 | 0.026 | 115/99 |
| learned (logistic, 4 features) | 0.265 | 1.20 | +0.20 | 0.554 | +0.034 | 4.3e-06 | 154/81 |
| oracle(true k) | 1.000 | 0.00 | +0.00 | 0.582 | +0.062 | 4e-05 | 207/167 |

### ALL trees, noise2 (n=446)
| rule | exact-k acc | mean abs(k_hat-k) | bias | mean Jaccard | dJac vs paper | Wilcoxon p | wins/losses |
|---|---|---|---|---|---|---|---|
| paper(min_p=0.01) | 0.087 | 2.11 | +1.41 | 0.328 | +0.000 | nan | 0/0 |
| bic(lam=2) | 0.020 | 2.71 | +2.62 | 0.325 | -0.004 | 0.6 | 111/132 |
| ratio(tau=0.9) | 0.007 | 2.95 | +2.94 | 0.312 | -0.017 | 0.00094 | 84/157 |
| paper(min_p tuned=0.005) | 0.056 | 2.43 | +2.20 | 0.331 | +0.003 | 0.45 | 62/75 |
| adj (paper+adjacency stop) | 0.105 | 2.04 | +1.22 | 0.341 | +0.013 | 7.7e-05 | 32/8 |
| adj(min_p tuned=0.005) | 0.076 | 2.31 | +1.88 | 0.348 | +0.020 | 0.00025 | 90/69 |
| learned (logistic, 4 features) | 0.132 | 1.82 | +0.52 | 0.323 | -0.006 | 0.24 | 87/79 |
| oracle(true k) | 1.000 | 0.00 | +0.00 | 0.287 | -0.042 | 4.3e-10 | 128/257 |

### ALL trees, all (n=1338)
| rule | exact-k acc | mean abs(k_hat-k) | bias | mean Jaccard | dJac vs paper | Wilcoxon p | wins/losses |
|---|---|---|---|---|---|---|---|
| paper(min_p=0.01) | 0.326 | 1.38 | +0.79 | 0.590 | +0.000 | nan | 0/0 |
| bic(lam=2) | 0.248 | 2.01 | +1.96 | 0.557 | -0.033 | 6.4e-17 | 255/482 |
| ratio(tau=0.9) | 0.325 | 1.97 | +1.94 | 0.574 | -0.016 | 5.6e-07 | 236/426 |
| paper(min_p tuned=0.005) | 0.339 | 1.57 | +1.36 | 0.595 | +0.004 | 0.15 | 166/187 |
| adj (paper+adjacency stop) | 0.335 | 1.32 | +0.57 | 0.601 | +0.011 | 3e-07 | 106/44 |
| adj(min_p tuned=0.005) | 0.351 | 1.46 | +1.03 | 0.609 | +0.019 | 4.2e-08 | 263/190 |
| learned (logistic, 4 features) | 0.422 | 1.10 | +0.15 | 0.603 | +0.013 | 9.1e-06 | 306/181 |
| oracle(true k) | 1.000 | 0.00 | +0.00 | 0.609 | +0.018 | 0.13 | 426/446 |

### Leave-one-tree-out: learned rule vs paper rule on ALL trees
| noise | n | paper exact-k | LOTO exact-k | paper mean abs err | LOTO mean abs err | paper Jaccard | LOTO Jaccard | dJac | Wilcoxon p | wins/losses | sign-test p (exact-k) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| noise0 | 446 | 0.747 | 0.843 | 0.34 | 0.36 | 0.922 | 0.922 | +0.000 | 0.29 | 62/31 | 1.3e-06 (61 vs 18) |
| noise1 | 446 | 0.143 | 0.291 | 1.70 | 1.15 | 0.520 | 0.562 | +0.042 | 2.6e-07 | 178/93 | 1.6e-09 (94 vs 28) |
| noise2 | 446 | 0.087 | 0.148 | 2.11 | 1.83 | 0.328 | 0.313 | -0.016 | 0.012 | 100/114 | 0.00036 (41 vs 14) |
| all | 1338 | 0.326 | 0.428 | 1.38 | 1.11 | 0.590 | 0.599 | +0.009 | 0.0029 | 340/238 | 5e-18 (196 vs 60) |

### TEST trees, noisy only (noise1+noise2): exact-k accuracy by true k
| true k | n | paper(min_p=0.01) | bic(lam=2) | ratio(tau=0.9) | paper(min_p tuned=0.005) | adj (paper+adjacency stop) | adj(min_p tuned=0.005) | learned (logistic, 4 features) | oracle(true k) |
|---|---|---|---|---|---|---|---|---|---|
| 2 | 80 | 0.05 | 0.00 | 0.00 | 0.03 | 0.05 | 0.03 | 0.11 | 1.00 |
| 3 | 80 | 0.09 | 0.00 | 0.00 | 0.01 | 0.15 | 0.07 | 0.17 | 1.00 |
| 5 | 80 | 0.07 | 0.00 | 0.00 | 0.05 | 0.09 | 0.07 | 0.20 | 1.00 |
| 7 | 78 | 0.19 | 0.00 | 0.00 | 0.09 | 0.18 | 0.12 | 0.33 | 1.00 |
| 10 | 74 | 0.19 | 0.03 | 0.00 | 0.16 | 0.18 | 0.16 | 0.14 | 1.00 |

true placements with an adjacent pair: 189/1338 (the authors' generator excludes near placements)
