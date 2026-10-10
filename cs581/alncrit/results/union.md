Replicates in both sets: 20 (1000L1/R0, 1000L1/R1, 1000L1/R2, 1000L2/R0, 1000L2/R1, 1000L3/R0, 1000L3/R3, 1000M2/R0, 1000M3/R0, 1000M3/R1, 1000M3/R3, 1000M4/R0, 1000M4/R1, 1000S1/R0, 1000S1/R1, 1000S1/R2, 1000S2/R0, 1000S2/R1, 1000S3/R0, RNASim/R0); alignments per replicate: MAGUS, MAGUS(Slow), MAGUS(Slow)-rerun, MAGUS-rerun, PASTA(1), PASTA(3), PASTA(3)+GCM, TRUE, slow-soft-m3

#### Union, incl. TRUE (n=180)

| model | coefs | cluster-robust p | within R² |
|---|---|---|---|
| SPFN | 0.183 | 3.3e-05 | 0.344 |
| SPFP | 0.189 | 6.4e-05 | 0.347 |
| log_comp | -2.204 | 0.062 | 0.064 |
| TC_err | 0.015 | 6.8e-05 | 0.182 |
| SPFN + SPFP | 0.022, 0.166 | 0.88, 0.37 | 0.347 |
| SPFN + SPFP + log_comp | 0.039, 0.146, -0.270 | 0.79, 0.41, 0.6 | 0.348 |

within-replicate correlation of SPFN and SPFP: r = 0.994; SPFN vs log compression r = -0.367

#### Union, estimated only (n=160)

| model | coefs | cluster-robust p | within R² |
|---|---|---|---|
| SPFN | 0.256 | 0.0017 | 0.271 |
| SPFP | 0.252 | 0.0019 | 0.276 |
| log_comp | -1.771 | 0.054 | 0.053 |
| TC_err | 0.085 | 0.024 | 0.203 |
| SPFN + SPFP | 0.057, 0.197 | 0.82, 0.42 | 0.276 |
| SPFN + SPFP + log_comp | 0.032, 0.256, 0.941 | 0.89, 0.26, 0.18 | 0.286 |

within-replicate correlation of SPFN and SPFP: r = 0.984; SPFN vs log compression r = -0.567

