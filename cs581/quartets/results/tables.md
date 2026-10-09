## Mean error (normalized RF, %) per method; held-out replicates 06-20

| data | genes | condition | n | astral4 | wastral | treeqmc | wqfm | astrid | ls[astral4] | ls[treeqmc] | ls[wqfm] | ls[astrid] | harvest | capminor | vote | reweight |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ils | 50-raxml | ILS 2M-1E-7 | 15 | 18.0 | - | 16.5 | 11.6 | 13.0 | 17.4 | 17.4 | 17.4 | 17.4 | 17.4 | 17.4 | 17.9 | 17.9 |
| ils | 50-raxml | ILS 500K-1E-6 | 15 | 34.8 | - | 32.2 | - | - | - | - | - | - | - | - | - | - |
| ils | 50-raxml | ILS 500K-1E-7 | 15 | 31.9 | - | 33.9 | - | - | - | - | - | - | - | - | - | - |

## Mean normalized quartet score (x100) per method; held-out replicates

| data | genes | condition | n | astral4 | wastral | treeqmc | wqfm | astrid | ls[astral4] | ls[treeqmc] | ls[wqfm] | ls[astrid] | harvest | capminor | vote | reweight |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ils | 50-raxml | ILS 2M-1E-7 | 15 | 61.296 | - | 61.200 | 61.291 | 61.261 | 61.334 | 61.334 | 61.334 | 61.334 | 61.334 | 61.322 | 61.317 | 61.254 |
| ils | 50-raxml | ILS 500K-1E-6 | 15 | 47.073 | - | 46.949 | - | - | - | - | - | - | - | - | - | - |
| ils | 50-raxml | ILS 500K-1E-7 | 15 | 47.270 | - | 47.112 | - | - | - | - | - | - | - | - | - | - |

## Paired comparisons, pooled over conditions within each dataset/gene count (held-out reps)

Δerr = mean(candidate − baseline) normalized RF in percentage points (negative = candidate better). W/T/L = candidate better / tie (identical RF) / worse. p = two-sided Wilcoxon signed-rank (ties dropped). Δscore = mean normalized quartet score difference (x100); ↑/↓ = #reps where candidate score is higher/lower.

| data | genes | candidate | baseline | n | Δerr (pp) | W/T/L | p | Δscore x100 | ↑/↓ |
|---|---|---|---|---|---|---|---|---|---|
| ils | 50-raxml | ls[astral4] | astral4 | 9 | +0.00 | 0/9/0 | n/a | +0.0000 | 0/0 |
| ils | 50-raxml | ls[treeqmc] | treeqmc | 9 | +2.90 | 1/4/4 | 0.25 | +0.0776 | 7/0 |
| ils | 50-raxml | ls[wqfm] | wqfm | 9 | +5.80 | 0/4/5 | 0.062 | +0.0434 | 6/0 |
| ils | 50-raxml | ls[astrid] | astrid | 9 | +4.35 | 1/2/6 | 0.2 | +0.0735 | 8/0 |
| ils | 50-raxml | harvest | astral4 | 9 | +0.00 | 0/9/0 | n/a | +0.0000 | 0/0 |
| ils | 50-raxml | capminor | astral4 | 9 | +0.00 | 1/7/1 | 1 | -0.0127 | 0/3 |
| ils | 50-raxml | vote | astral4 | 9 | +0.48 | 1/6/2 | 1 | -0.0169 | 0/4 |
| ils | 50-raxml | reweight | astral4 | 9 | +0.48 | 3/3/3 | 1 | -0.0803 | 0/7 |
| ils | 50-raxml | treeqmc | astral4 | 45 | -0.68 | 14/19/12 | 0.68 | -0.1259 | 0/39 |
| ils | 50-raxml | wqfm | astral4 | 9 | -5.80 | 5/4/0 | 0.062 | -0.0434 | 0/6 |
| ils | 50-raxml | astrid | astral4 | 9 | -4.35 | 6/2/1 | 0.2 | -0.0735 | 0/8 |

## Per-condition paired comparisons for the HGT-aware variants (held-out reps)

| data | genes | condition | candidate | n | Δerr (pp) | W/T/L | p |
|---|---|---|---|---|---|---|---|
| ils | 50-raxml | ILS 2M-1E-7 | capminor | 9 | +0.00 | 1/7/1 | 1 |
| ils | 50-raxml | ILS 2M-1E-7 | vote | 9 | +0.48 | 1/6/2 | 1 |
| ils | 50-raxml | ILS 2M-1E-7 | reweight | 9 | +0.48 | 3/3/3 | 1 |
| ils | 50-raxml | ILS 2M-1E-7 | harvest | 9 | +0.00 | 0/9/0 | n/a |
| ils | 50-raxml | ILS 2M-1E-7 | ls[treeqmc] | 9 | +2.90 | 1/4/4 | 0.25 |

## Development replicates 01-05 (pooled; not used for claims)

| data | genes | candidate | baseline | n | Δerr (pp) | W/T/L |
|---|---|---|---|---|---|---|
| hgt | 200-est | ls[astral4] | astral4 | 1 | +0.00 | 0/1/0 |
| hgt | 200-est | ls[treeqmc] | treeqmc | 1 | +0.00 | 0/1/0 |
| hgt | 200-est | ls[wqfm] | wqfm | 1 | +0.00 | 0/1/0 |
| hgt | 200-est | ls[astrid] | astrid | 1 | +0.00 | 0/1/0 |
| hgt | 200-est | harvest | astral4 | 1 | +0.00 | 0/1/0 |
| hgt | 200-est | capminor | astral4 | 1 | +0.00 | 0/1/0 |
| hgt | 200-est | vote | astral4 | 1 | +0.00 | 0/1/0 |
| hgt | 200-est | reweight | astral4 | 1 | +0.00 | 0/1/0 |
| hgt | 50-est | ls[astral4] | astral4 | 8 | +0.00 | 0/8/0 |
| hgt | 50-est | ls[treeqmc] | treeqmc | 8 | -0.26 | 2/5/1 |
| hgt | 50-est | ls[wqfm] | wqfm | 8 | -0.52 | 2/5/1 |
| hgt | 50-est | ls[astrid] | astrid | 7 | -1.79 | 3/3/1 |
| hgt | 50-est | harvest | astral4 | 7 | +0.00 | 0/7/0 |
| hgt | 50-est | capminor | astral4 | 4 | -0.52 | 1/3/0 |
| hgt | 50-est | vote | astral4 | 4 | -0.52 | 1/3/0 |
| hgt | 50-est | reweight | astral4 | 4 | +0.52 | 0/3/1 |
| ils | 200-raxml | ls[astral4] | astral4 | 1 | +0.00 | 0/1/0 |
| ils | 200-raxml | ls[treeqmc] | treeqmc | 1 | +8.70 | 0/0/1 |
| ils | 200-raxml | ls[wqfm] | wqfm | 1 | +8.70 | 0/0/1 |
| ils | 200-raxml | ls[astrid] | astrid | 1 | -8.70 | 1/0/0 |
| ils | 200-raxml | harvest | astral4 | 1 | +0.00 | 0/1/0 |
| ils | 200-raxml | capminor | astral4 | 1 | +0.00 | 0/1/0 |

## Mean runtime (s, one thread) per method, all replicates

| data | genes | astral4 | wastral | treeqmc | wqfm | astrid | ls[astral4] | ls[treeqmc] | ls[wqfm] | ls[astrid] | harvest | capminor | vote | reweight |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hgt | 200-est | 4.1 | 6.0 | 2.4 | 3.2 | 0.5 | 11.5 | 11.3 | 11.4 | 11.2 | 4.1 | 11.4 | 11.1 | 8.1 |
| hgt | 50-est | 1.1 | 1.5 | 1.7 | 2.3 | 0.2 | 13.4 | 17.3 | 16.3 | 15.9 | 1.3 | 15.3 | 16.6 | 3.3 |
| ils | 200-raxml | 1.9 | - | 1.2 | 2.6 | 0.2 | 0.3 | 0.3 | 0.4 | 0.6 | 1.4 | 0.3 | - | - |
| ils | 50-raxml | 0.4 | - | 0.6 | 0.9 | 0.1 | 0.2 | 0.3 | 0.3 | 0.3 | 0.3 | 0.2 | 0.2 | 0.4 |

Replicates loaded: 70
