## Subset level: Δ error vs MAFFT L-INS-i (MAGUS's command, rerun here), points; negative = more accurate than L-INS-i

| kind | tool | n | mean err | mean Δ | W/T/L | p | mean s (1 core) |
|---|---|---|---|---|---|---|---|
| magus40 | cached-linsi | 18 | 3.6 | +0.04 | 0/17/1 | 0.11 | – |
| magus40 | mafft-linsi-magus | 18 | 3.5 | +0.00 | 0/18/0 | – | 11 |
| magus40 | mafft-ginsi | 18 | 3.5 | -0.06 | 4/6/8 | 0.36 | 11 |
| magus40 | mafft-einsi | 18 | 6.5 | +2.92 | 2/4/12 | 0.0031 | 14 |
| magus40 | famsa | 18 | 16.7 | +13.20 | 0/1/17 | 7.6e-06 | 1 |
| magus40 | twilight-1 | 18 | 43.7 | +40.12 | 0/0/18 | 7.6e-06 | 2 |
| magus40 | twilight | 18 | 38.2 | +34.69 | 0/0/18 | 7.6e-06 | 8 |
| magus40 | clustalo | 18 | 8.7 | +5.12 | 0/0/18 | 7.6e-06 | 3 |
| magus40 | mafft-auto | 18 | 3.9 | +0.36 | 0/11/7 | 0.012 | 11 |
| rand40 | mafft-linsi-magus | 10 | 31.2 | +0.00 | 0/10/0 | – | 24 |
| rand40 | mafft-ginsi | 10 | 34.1 | +2.89 | 5/0/5 | 0.32 | 29 |
| rand40 | mafft-einsi | 10 | 59.4 | +28.13 | 1/0/9 | 0.0039 | 51 |
| rand40 | famsa | 10 | 83.2 | +51.93 | 0/0/10 | 0.002 | 1 |
| rand40 | twilight-1 | 10 | 86.7 | +55.45 | 0/0/10 | 0.002 | 3 |
| rand40 | twilight | 10 | 86.3 | +55.05 | 0/0/10 | 0.002 | 10 |
| rand40 | clustalo | 10 | 72.6 | +41.35 | 0/0/10 | 0.002 | 4 |
| rand40 | mafft-auto | 10 | 41.6 | +10.32 | 2/1/7 | 0.037 | 35 |

### magus40 by data type (mean Δ vs L-INS-i, points)

| tool | ROSE | RNASim | 16S | BAliBASE |
|---|---|---|---|---|
| cached-linsi | +0.04 (n=18) |  |  |  |
| mafft-linsi-magus | +0.00 (n=18) |  |  |  |
| mafft-ginsi | -0.06 (n=18) |  |  |  |
| mafft-einsi | +2.92 (n=18) |  |  |  |
| famsa | +13.20 (n=18) |  |  |  |
| twilight-1 | +40.12 (n=18) |  |  |  |
| twilight | +34.69 (n=18) |  |  |  |
| clustalo | +5.12 (n=18) |  |  |  |
| mafft-auto | +0.36 (n=18) |  |  |  |
