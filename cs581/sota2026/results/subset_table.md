## Subset level: Δ error vs MAFFT L-INS-i (MAGUS's command, rerun here), points; negative = more accurate than L-INS-i

| kind | tool | n | mean err | mean Δ | W/T/L | p | mean s (1 core) |
|---|---|---|---|---|---|---|---|
| magus40 | cached-linsi | 30 | 6.3 | +0.01 | 2/26/2 | 0.39 | – |
| magus40 | mafft-linsi-magus | 30 | 6.2 | +0.00 | 0/30/0 | – | 8 |
| magus40 | mafft-ginsi | 30 | 6.2 | -0.05 | 8/12/10 | 0.72 | 8 |
| magus40 | mafft-einsi | 30 | 8.1 | +1.88 | 4/9/17 | 0.00069 | 10 |
| magus40 | famsa | 30 | 16.6 | +10.36 | 1/1/28 | 9.3e-09 | 2 |
| magus40 | twilight-1 | 30 | 31.3 | +25.02 | 2/0/28 | 9.3e-09 | 2 |
| magus40 | twilight | 30 | 28.2 | +21.95 | 1/0/29 | 3.7e-09 | 7 |
| magus40 | clustalo | 30 | 9.5 | +3.28 | 5/0/25 | 6.9e-06 | 3 |
| magus40 | mafft-auto | 30 | 6.5 | +0.26 | 2/17/11 | 0.0038 | 9 |
| magus40 | muscle5 | 16 | 9.3 | +0.76 | 6/2/8 | 0.27 | 96 |
| rand40 | mafft-linsi-magus | 16 | 27.6 | +0.00 | 0/16/0 | – | 20 |
| rand40 | mafft-ginsi | 16 | 29.3 | +1.65 | 10/0/6 | 0.9 | 24 |
| rand40 | mafft-einsi | 16 | 45.2 | +17.59 | 3/2/11 | 0.0042 | 37 |
| rand40 | famsa | 16 | 63.9 | +36.29 | 0/0/16 | 3.1e-05 | 2 |
| rand40 | twilight-1 | 16 | 64.7 | +37.02 | 0/0/16 | 3.1e-05 | 3 |
| rand40 | twilight | 16 | 64.8 | +37.17 | 0/0/16 | 3.1e-05 | 11 |
| rand40 | clustalo | 16 | 54.6 | +26.95 | 1/0/15 | 6.1e-05 | 4 |
| rand40 | mafft-auto | 16 | 34.0 | +6.41 | 6/1/9 | 0.12 | 28 |
| bb200 | cached-linsi | 1 | 30.8 | -0.27 | 1/0/0 | – | – |
| bb200 | mafft-linsi-magus | 1 | 31.0 | +0.00 | 0/1/0 | – | 377 |
| bb200 | mafft-ginsi | 1 | 33.4 | +2.34 | 0/0/1 | – | 350 |
| bb200 | famsa | 1 | 94.1 | +63.06 | 0/0/1 | – | 4 |

### magus40 by data type (mean Δ vs L-INS-i, points)

| tool | ROSE | RNASim | 16S | BAliBASE |
|---|---|---|---|---|
| cached-linsi | +0.04 (n=18) | +0.00 (n=2) | +0.02 (n=2) | -0.05 (n=8) |
| mafft-linsi-magus | +0.00 (n=18) | +0.00 (n=2) | +0.00 (n=2) | +0.00 (n=8) |
| mafft-ginsi | -0.06 (n=18) | -0.14 (n=2) | -0.02 (n=2) | -0.00 (n=8) |
| mafft-einsi | +2.92 (n=18) | +0.01 (n=2) | +0.09 (n=2) | +0.46 (n=8) |
| famsa | +13.20 (n=18) | +31.68 (n=2) | +0.17 (n=2) | +1.19 (n=8) |
| twilight-1 | +40.12 (n=18) | +1.20 (n=2) | +0.38 (n=2) | +3.18 (n=8) |
| twilight | +34.69 (n=18) | +0.88 (n=2) | +0.43 (n=2) | +3.92 (n=8) |
| clustalo | +5.12 (n=18) | +0.54 (n=2) | +3.41 (n=2) | -0.20 (n=8) |
| mafft-auto | +0.36 (n=18) | +0.04 (n=2) | +0.00 (n=2) | +0.15 (n=8) |
| muscle5 | +3.25 (n=4) | -0.67 (n=2) | -0.00 (n=2) | +0.06 (n=8) |
