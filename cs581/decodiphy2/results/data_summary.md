Runs analysed: 4171 (authors' stored E1 runs, 12 trees, k in {2,3,5,7,10}, seeds 1-10, noise0 + 6 noisy settings)

## 1. How often are configurations non-identifiable?

| noise | k | runs | TRUE: adjacent | TRUE: kernel (claw etc.) | FIT: adjacent | FIT: kernel | FIT: boundary x | FIT: class width > 0 | FIT: any flag |
|---|---|---|---|---|---|---|---|---|---|
| all | all | 4171 | 12.1% | 4.0% | 2.7% | 2.2% | 4.6% | 4.7% | 8.9% |
| noise0 | 2 | 120 | 2.5% | 0.0% | 2.5% | 0.0% | 0.0% | 2.5% | 2.5% |
| noise0 | 3 | 120 | 5.0% | 0.0% | 5.8% | 0.8% | 1.7% | 5.8% | 6.7% |
| noise0 | 5 | 120 | 7.5% | 0.8% | 5.8% | 2.5% | 3.3% | 8.3% | 10.0% |
| noise0 | 7 | 119 | 15.1% | 5.0% | 10.9% | 6.7% | 3.4% | 16.0% | 17.6% |
| noise0 | 10 | 117 | 30.8% | 14.5% | 16.2% | 9.4% | 2.6% | 24.8% | 26.5% |
| noise0 | all | 596 | 12.1% | 4.0% | 8.2% | 3.9% | 2.2% | 11.4% | 12.6% |
| noise1 | 2 | 360 | 2.5% | 0.0% | 0.8% | 0.0% | 1.4% | 0.8% | 2.2% |
| noise1 | 3 | 360 | 5.0% | 0.0% | 0.8% | 0.0% | 1.7% | 0.8% | 2.5% |
| noise1 | 5 | 359 | 7.5% | 0.8% | 0.6% | 0.6% | 3.9% | 1.1% | 4.7% |
| noise1 | 7 | 357 | 15.1% | 5.0% | 3.4% | 1.4% | 11.5% | 4.8% | 16.0% |
| noise1 | 10 | 351 | 30.8% | 14.5% | 3.1% | 6.0% | 13.4% | 8.8% | 21.4% |
| noise1 | all | 1787 | 12.1% | 4.0% | 1.7% | 1.6% | 6.3% | 3.2% | 9.3% |
| noise2 | 2 | 360 | 2.5% | 0.0% | 0.6% | 0.0% | 2.2% | 0.6% | 2.8% |
| noise2 | 3 | 360 | 5.0% | 0.0% | 1.1% | 0.0% | 2.2% | 1.1% | 3.3% |
| noise2 | 5 | 360 | 7.5% | 0.8% | 1.1% | 1.4% | 3.6% | 2.5% | 5.8% |
| noise2 | 7 | 357 | 15.1% | 5.0% | 2.8% | 3.4% | 3.4% | 5.0% | 7.8% |
| noise2 | 10 | 351 | 30.8% | 14.5% | 3.7% | 6.8% | 7.4% | 10.3% | 16.8% |
| noise2 | all | 1788 | 12.1% | 4.0% | 1.8% | 2.3% | 3.7% | 3.9% | 7.3% |

## 2. Fits whose equivalence class is non-trivial: 195 of 4171 (4.7%)

Max abundance range within the class (max_q hi_q - lo_q): median 0.100, mean 0.125, 90th pct 0.231, max 0.823.
Share of flagged fits with range > 0.05: 75.4%; > 0.10: 49.2%.

### Accuracy: DecoDiPhy's point vs the canonical class representative (same loss, same d)

| subset | n | EMD DecoDiPhy | EMD canonical | Wilcoxon p | canonical better / worse | EMD range over class vertices (mean max-min) |
|---|---|---|---|---|---|---|
| all flagged | 195 | 0.00444 | 0.00435 | 0.2 | 111/84 | 0.00320 |
| flagged, noise0 | 68 | 0.00349 | 0.00336 | 0.63 | 37/31 | 0.00317 |
| flagged, noise1 | 58 | 0.00556 | 0.00543 | 0.3 | 34/24 | 0.00343 |
| flagged, noise2 | 69 | 0.00444 | 0.00443 | 0.48 | 40/29 | 0.00304 |
| flagged, adjacent fits | 113 | 0.00372 | 0.00369 | 0.34 | 64/49 | 0.00325 |
| flagged, non-adjacent (kernel) fits | 82 | 0.00544 | 0.00526 | 0.39 | 47/35 | 0.00313 |

### Flagged fits with exactly the true edge set: 73

Abundance L1 error: DecoDiPhy 0.1341, canonical 0.1359, best point in class (lower bound) 0.0011.
canonical better/worse: 37/36, Wilcoxon p = 0.77.
True p inside the reported per-placement intervals [lo, hi] (+-0.001): 86.3%; DecoDiPhy's point equal to the true p (+-0.001): 0.0%.
- noise0: n=41, L1 DecoDiPhy 0.1189 / canonical 0.1200 / class lb 0.0009; truth in intervals 87.8%, point exact 0.0%
- noise1: n=10, L1 DecoDiPhy 0.2151 / canonical 0.1925 / class lb 0.0004; truth in intervals 90.0%, point exact 0.0%
- noise2: n=22, L1 DecoDiPhy 0.1255 / canonical 0.1399 / class lb 0.0018; truth in intervals 81.8%, point exact 0.0%

## 3. Whole-benchmark effect of reporting the canonical representative instead of DecoDiPhy's point

- noise0: n=596, mean EMD 0.00086 -> 0.00085 (-1.82%)
- noise1: n=1787, mean EMD 0.00402 -> 0.00402 (-0.10%)
- noise2: n=1788, mean EMD 0.00249 -> 0.00249 (-0.03%)
- all: n=4171, mean EMD 0.00291 -> 0.00291 (-0.15%)

## 4. The flag predicts abundance error (runs with exactly the true edge set)

| noise | class trivial: n | abundance L1 mean (median) | class non-trivial: n | abundance L1 mean (median) | share of total L1 error in flagged runs |
|---|---|---|---|---|---|
| noise0 | 460 | 0.0003 (0.0000) | 41 | 0.1189 (0.0962) | 97% |
| noise1 | 1129 | 0.0044 (0.0007) | 10 | 0.2151 (0.1920) | 30% |
| noise2 | 1130 | 0.0040 (0.0006) | 22 | 0.1255 (0.1002) | 38% |
