## Full datasets: average error (SPFN+SPFP)/2, %, and wall-clock seconds (4 threads)

MAGUS(pub)/PASTA(pub) = FastSP on the authors' published alignment of the same replicate; their seconds are the paper's own timing (different hardware). `magus` = our rerun with the paper's flags.

| dataset | MAGUS(pub) | PASTA(pub) | famsa | twilight-1 | twilight | mafft-parttree | mafft-auto |
|---|---|---|---|---|---|---|---|
| 1000L1_R0 | 8.1 (716s) | 8.3 (1844s) | 33.1 (48s) | 97.9 (17s) | 98.1 (136s) | 98.0 (39s) | 99.0 (18s) |
| 1000L2_R0 | 3.8 (673s) | 6.3 (1615s) | 37.9 (223s) | 97.8 (20s) | 98.1 (243s) | 98.0 (39s) | 99.3 (20s) |
| 1000L3_R0 | 12.7 (873s) | 15.8 (2486s) | 63.6 (344s) | 98.0 (21s) | 98.2 (253s) | 98.3 (39s) | 99.3 (20s) |

## Paired comparison vs published MAGUS(Fast), same replicate (Δ = tool − MAGUS, points; W/T/L = tool better/tie(±0.05)/worse; Wilcoxon signed-rank)

| data | tool | n | mean Δ | W/T/L | p |
|---|---|---|---|---|---|
| ROSE | famsa | 3 | +36.7 | 0/0/3 | – |
| ROSE | twilight-1 | 3 | +89.7 | 0/0/3 | – |
| ROSE | twilight | 3 | +90.0 | 0/0/3 | – |
| ROSE | mafft-parttree | 3 | +89.9 | 0/0/3 | – |
| ROSE | mafft-auto | 3 | +91.0 | 0/0/3 | – |
| ROSE | PASTA(pub) | 3 | +1.9 | 0/0/3 | – |
| all | famsa | 3 | +36.7 | 0/0/3 | – |
| all | twilight-1 | 3 | +89.7 | 0/0/3 | – |
| all | twilight | 3 | +90.0 | 0/0/3 | – |
| all | mafft-parttree | 3 | +89.9 | 0/0/3 | – |
| all | mafft-auto | 3 | +91.0 | 0/0/3 | – |
| all | PASTA(pub) | 3 | +1.9 | 0/0/3 | – |
