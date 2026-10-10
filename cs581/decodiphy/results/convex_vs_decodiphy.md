runs: 135; tau tuned on train trees: 0.2

| set | noise | n | DecoDiPhy (paper k) EMD | DecoDiPhy (oracle k) EMD | convex max-ybar EMD | dEMD (convex - paper) | Wilcoxon p | convex better / worse | median support (m>0.01) vs true 2k |
|---|---|---|---|---|---|---|---|---|---|
| test | noise0 | 20 | 0.0004 | 0.0000 | 0.0002 | -0.0001 | 0.07 | 3/17 | 9 vs 10 |
| test | noise1 | 20 | 0.0217 | 0.0222 | 0.0226 | +0.0009 | 0.0064 | 4/16 | 14 vs 10 |
| test | noise2 | 20 | 0.0413 | 0.0430 | 0.0441 | +0.0028 | 0.011 | 3/17 | 14 vs 10 |
| all | noise0 | 45 | 0.0003 | 0.0000 | 0.0005 | +0.0001 | 0.0035 | 7/38 | 9 vs 10 |
| all | noise1 | 45 | 0.0185 | 0.0180 | 0.0193 | +0.0008 | 3.6e-05 | 7/38 | 12 vs 10 |
| all | noise2 | 45 | 0.0360 | 0.0379 | 0.0388 | +0.0028 | 3e-07 | 8/37 | 15 vs 10 |

mean runtime per instance: convex (2 + 4 solves, all taus) 1.71s; DecoDiPhy greedy up to chosen k+1 2.75s
tau=0.0: mean EMD all=0.0195
tau=0.05: mean EMD all=0.0195
tau=0.2: mean EMD all=0.0195
tau=0.5: mean EMD all=0.0195
tau=1.0: mean EMD all=0.0195
