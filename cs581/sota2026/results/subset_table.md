## Subset level: Δ error vs MAFFT L-INS-i (MAGUS's command, rerun here), points; negative = more accurate than L-INS-i

| kind | tool | n | mean err | mean Δ | W/T/L | p | mean s (1 core) |
|---|---|---|---|---|---|---|---|
| magus40 | cached-linsi | 8 | 2.3 | +0.00 | 0/8/0 | 1 | – |
| magus40 | mafft-linsi-magus | 8 | 2.3 | +0.00 | 0/8/0 | – | 11 |
| magus40 | mafft-ginsi | 8 | 2.6 | +0.33 | 1/1/6 | 0.2 | 11 |
| magus40 | mafft-einsi | 8 | 5.2 | +2.91 | 1/1/6 | 0.023 | 12 |
| magus40 | famsa | 8 | 18.9 | +16.63 | 0/0/8 | 0.0078 | 1 |
| magus40 | twilight-1 | 8 | 53.1 | +50.87 | 0/0/8 | 0.0078 | 2 |
| magus40 | twilight | 8 | 46.0 | +43.79 | 0/0/8 | 0.0078 | 8 |
| magus40 | clustalo | 8 | 7.8 | +5.58 | 0/0/8 | 0.0078 | 3 |
| magus40 | mafft-auto | 8 | 2.7 | +0.48 | 0/3/5 | 0.062 | 11 |
| rand40 | mafft-linsi-magus | 5 | 45.6 | +0.00 | 0/5/0 | – | 28 |
| rand40 | mafft-ginsi | 4 | 52.6 | +4.67 | 2/0/2 | – | 43 |
| rand40 | mafft-einsi | 4 | 76.6 | +28.72 | 0/0/4 | – | 68 |
| rand40 | famsa | 4 | 96.8 | +48.90 | 0/0/4 | – | 1 |
| rand40 | twilight-1 | 4 | 96.7 | +48.86 | 0/0/4 | – | 3 |
| rand40 | twilight | 4 | 96.6 | +48.69 | 0/0/4 | – | 11 |
| rand40 | clustalo | 4 | 91.0 | +43.14 | 0/0/4 | – | 5 |
| rand40 | mafft-auto | 4 | 58.6 | +10.72 | 1/0/3 | – | 48 |

### magus40 by data type (mean Δ vs L-INS-i, points)

| tool | ROSE | RNASim | 16S | BAliBASE |
|---|---|---|---|---|
| cached-linsi | +0.00 (n=8) |  |  |  |
| mafft-linsi-magus | +0.00 (n=8) |  |  |  |
| mafft-ginsi | +0.33 (n=8) |  |  |  |
| mafft-einsi | +2.91 (n=8) |  |  |  |
| famsa | +16.63 (n=8) |  |  |  |
| twilight-1 | +50.87 (n=8) |  |  |  |
| twilight | +43.79 (n=8) |  |  |  |
| clustalo | +5.58 (n=8) |  |  |  |
| mafft-auto | +0.48 (n=8) |  |  |  |
