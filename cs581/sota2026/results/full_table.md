## Full datasets: average error (SPFN+SPFP)/2, %, and wall-clock seconds (4 threads)

MAGUS(pub)/PASTA(pub) = FastSP on the authors' published alignment of the same replicate; their seconds are the paper's own timing (different hardware). MAGUS(4c) = earlier rerun of MAGUS(Fast) with the paper's flags on this 4-core machine type (cs581/experiments/runs/<rep>/prep.json; random backbones, so not identical to the published alignment).

| dataset | MAGUS(pub) | MAGUS(4c) | PASTA(pub) | famsa | twilight-1 | twilight | mafft-parttree | mafft-auto | mafft-linsi | muscle5 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1000L1_R0 | 8.1 (716s) | 7.5 (1870s) | 8.3 (1844s) | 33.1 (48s) | 97.9 (17s) | 98.1 (136s) | 98.0 (39s) | 99.0 (18s) |  |  |
| 1000L2_R0 | 3.8 (673s) | 4.5 (1500s) | 6.3 (1615s) | 37.9 (223s) | 97.8 (20s) | 98.1 (243s) | 98.0 (39s) | 99.3 (20s) |  |  |
| 1000L3_R0 | 12.7 (873s) | 11.3 (1737s) | 15.8 (2486s) | 63.6 (344s) | 98.0 (21s) | 98.2 (253s) | 98.3 (39s) | 99.3 (20s) |  |  |
| RNASim_R0 | 9.6 (2227s) | 9.9 (888s) | 10.1 (3102s) | fail rc=-9 (357s) | 14.4 (21s) | 11.1 (230s) | 23.7 (60s) | 22.0 (26s) |  |  |
| BBA0039_R0 | 4.5 (157s) | 4.7 (244s) | 4.3 (462s) | 5.5 (56s) | 6.9 (8s) | 8.0 (100s) | 6.7 (2s) | 5.5 (1s) |  |  |
| BBA0067_R0 | 25.6 (371s) | 26.3 (1051s) | 25.7 (373s) | 32.0 (87s) | 42.9 (6s) | 37.3 (54s) | 38.4 (2s) | 33.2 (13s) | 26.1 (367s) | timeout (900s) |
| BBA0081_R0 | 57.0 (1012s) |  | 74.1 (548s) | 72.2 (96s) | 75.9 (6s) | 72.3 (45s) | 75.7 (2s) | 68.3 (45s) | 67.4 (588s) |  |
| BBA0101_R0 | 28.5 (477s) | 29.2 (829s) | 29.2 (381s) | 35.2 (88s) | 42.4 (7s) | 39.7 (65s) | 39.7 (3s) | 36.5 (14s) |  |  |
| BBA0117_R0 | 12.1 (14s) |  | 12.1 (35s) | 15.2 (11s) | 24.5 (0s) | 20.4 (45s) | 25.1 (0s) | 16.4 (1s) |  |  |
| BBA0134_R0 | 18.7 (845s) |  | 20.4 (573s) | 33.2 (83s) | 32.8 (6s) | 28.3 (68s) | 29.8 (2s) | 25.2 (13s) |  |  |
| BBA0154_R0 | 21.1 (562s) | 21.7 (1106s) | 22.5 (307s) | 25.3 (77s) | 25.4 (6s) | 26.3 (55s) | 25.7 (2s) | 25.1 (7s) |  |  |
| BBA0190_R0 | 23.1 (1117s) | 23.4 (2177s) | 23.9 (879s) | 33.7 (97s) | 26.0 (16s) | 30.4 (80s) | 26.9 (4s) | 26.5 (18s) |  |  |
| 16S.T_R0 | 9.9 (2918s) |  | 12.9 (6979s) | 17.6 (123s) | 12.7 (78s) |  | 18.8 (215s) | 19.4 (241s) |  |  |
| 1000M1_R0 | 11.5 (674s) |  | 13.0 (1993s) | 50.4 (773s) | 97.6 (30s) |  | 98.1 (39s) | 99.2 (25s) |  |  |
| 1000M2_R0 | 8.3 (968s) | 8.2 (1493s) | 13.7 (2007s) | 43.6 (471s) | 97.9 (30s) |  | 98.0 (38s) | 99.0 (23s) |  |  |
| 1000M3_R0 | 4.1 (586s) | 4.5 (1145s) | 5.6 (1328s) | 15.2 (254s) | 88.2 (22s) |  | 83.6 (33s) | 71.8 (17s) |  |  |
| 1000M4_R0 | 1.1 (408s) | 1.2 (992s) | 1.5 (1205s) | 2.5 (27s) | 3.3 (8s) |  | 5.5 (21s) | 5.4 (6s) |  |  |
| 1000S1_R0 | 8.5 (685s) | 9.8 (1563s) | 11.7 (1769s) | 26.5 (39s) | 96.8 (15s) |  | 97.1 (37s) | 98.8 (16s) |  |  |
| 1000S2_R0 | 3.5 (476s) | 4.7 (994s) | 5.7 (1530s) | 13.7 (37s) | 96.7 (15s) |  | 96.8 (37s) | 98.9 (16s) |  |  |
| 1000S3_R0 | 4.5 (533s) | 4.5 (1338s) | 4.9 (1343s) | 24.7 (38s) | 97.2 (16s) |  | 97.7 (38s) | 99.0 (19s) |  |  |

## Paired comparison vs published MAGUS(Fast), same replicate (Δ = tool − MAGUS, points; W/T/L = tool better/tie(±0.05)/worse; Wilcoxon signed-rank)

| data | tool | n | mean Δ | W/T/L | p |
|---|---|---|---|---|---|
| ROSE | famsa | 10 | +24.5 | 0/0/10 | 0.002 |
| ROSE | twilight-1 | 10 | +80.5 | 0/0/10 | 0.002 |
| ROSE | twilight | 3 | +90.0 | 0/0/3 | – |
| ROSE | mafft-parttree | 10 | +80.5 | 0/0/10 | 0.002 |
| ROSE | mafft-auto | 10 | +80.4 | 0/0/10 | 0.002 |
| ROSE | MAGUS(4c) | 9 | +0.2 | 2/2/5 | 0.57 |
| ROSE | PASTA(pub) | 10 | +2.1 | 0/0/10 | 0.002 |
| RNASim | twilight-1 | 1 | +4.8 | 0/0/1 | – |
| RNASim | twilight | 1 | +1.5 | 0/0/1 | – |
| RNASim | mafft-parttree | 1 | +14.1 | 0/0/1 | – |
| RNASim | mafft-auto | 1 | +12.4 | 0/0/1 | – |
| RNASim | MAGUS(4c) | 1 | +0.3 | 0/0/1 | – |
| RNASim | PASTA(pub) | 1 | +0.5 | 0/0/1 | – |
| BAliBASE | famsa | 8 | +7.7 | 0/0/8 | 0.0078 |
| BAliBASE | twilight-1 | 8 | +10.8 | 0/0/8 | 0.0078 |
| BAliBASE | twilight | 8 | +9.0 | 0/0/8 | 0.0078 |
| BAliBASE | mafft-parttree | 8 | +9.7 | 0/0/8 | 0.0078 |
| BAliBASE | mafft-auto | 8 | +5.7 | 0/0/8 | 0.0078 |
| BAliBASE | mafft-linsi | 2 | +5.4 | 0/0/2 | – |
| BAliBASE | MAGUS(4c) | 5 | +0.5 | 0/0/5 | 0.062 |
| BAliBASE | PASTA(pub) | 8 | +2.7 | 2/0/6 | 0.055 |
| 16S | famsa | 1 | +7.7 | 0/0/1 | – |
| 16S | twilight-1 | 1 | +2.7 | 0/0/1 | – |
| 16S | mafft-parttree | 1 | +8.9 | 0/0/1 | – |
| 16S | mafft-auto | 1 | +9.4 | 0/0/1 | – |
| 16S | PASTA(pub) | 1 | +3.0 | 0/0/1 | – |
| all | famsa | 19 | +16.5 | 0/0/19 | 3.8e-06 |
| all | twilight-1 | 20 | +45.0 | 0/0/20 | 1.9e-06 |
| all | twilight | 12 | +28.6 | 0/0/12 | 0.00049 |
| all | mafft-parttree | 20 | +45.3 | 0/0/20 | 1.9e-06 |
| all | mafft-auto | 20 | +43.6 | 0/0/20 | 1.9e-06 |
| all | mafft-linsi | 2 | +5.4 | 0/0/2 | – |
| all | MAGUS(4c) | 15 | +0.3 | 2/2/11 | 0.064 |
| all | PASTA(pub) | 20 | +2.3 | 2/0/18 | 1.9e-05 |

## Diagnostic: same aligner, true tree as guide tree (error %, seconds)

| dataset | MAGUS(pub) | famsa | famsa-truetree | twilight-1 | twilight | twilight-truetree |
|---|---|---|---|---|---|---|
| 1000L1_R0 | 8.1 (716s) | 33.1 (48s) | 16.3 (133s) | 97.9 (17s) | 98.1 (136s) | 29.4 (16s) |
| 1000L2_R0 | 3.8 (673s) | 37.9 (223s) | 14.5 (125s) | 97.8 (20s) | 98.1 (243s) | 34.9 (13s) |
| 1000L3_R0 | 12.7 (873s) | 63.6 (344s) | 27.2 (177s) | 98.0 (21s) | 98.2 (253s) | 40.9 (13s) |
| RNASim_R0 | 9.6 (2227s) | fail rc=-9 | fail rc=-9 | 14.4 (21s) | 11.1 (230s) | 10.1 (36s) |
| 1000M1_R0 | 11.5 (674s) | 50.4 (773s) | 21.0 (141s) | 97.6 (30s) |  | 33.9 (13s) |
| 1000M2_R0 | 8.3 (968s) | 43.6 (471s) | 22.4 (136s) | 97.9 (30s) |  | 23.9 (12s) |
| 1000M3_R0 | 4.1 (586s) | 15.2 (254s) | 7.2 (252s) | 88.2 (22s) |  | 8.8 (20s) |
| 1000M4_R0 | 1.1 (408s) | 2.5 (27s) | 2.3 (241s) | 3.3 (8s) |  | 1.1 (17s) |
| 1000S1_R0 | 8.5 (685s) | 26.5 (39s) | 13.8 (253s) | 96.8 (15s) |  | 24.5 (24s) |
| 1000S2_R0 | 3.5 (476s) | 13.7 (37s) | 7.4 (228s) | 96.7 (15s) |  | 11.7 (23s) |
| 1000S3_R0 | 4.5 (533s) | 24.7 (38s) | 10.2 (250s) | 97.2 (16s) |  | 22.9 (23s) |
