# Pre-registered held-out test on DISCO data (see PREREG_disco_heldout.md)


## Set A (reps 05-10: default, gdl_1e-9_1, gdl_1e-9_05): n = 18 runs

| comparison | mean FN (Pro) | mean FN (base) | mean diff | W/T/L | p | Holm-adjusted p |
|---|---|---|---|---|---|---|
| astrid-pro vs astrid-multi | 0.0471 | 0.0646 | -0.0176 | 12/5/1 | 0.00672 | 0.0202 |
| astrid-pro vs astrid-disco | 0.0471 | 0.0425 | +0.0045 | 4/5/9 | 0.131 | 0.263 |
| astrid-pro vs astral-pro | 0.0471 | 0.0454 | +0.0017 | 5/6/7 | 0.554 | 0.554 |
| (secondary) astrid-pro-w vs astrid-multi | 0.0459 | 0.0646 | -0.0187 | 10/8/0 | 0.00503 | - |
| (secondary) astrid-pro-w vs astrid-disco | 0.0459 | 0.0425 | +0.0034 | 3/10/5 | 0.325 | - |
| (secondary) astrid-pro-w vs astral-pro | 0.0459 | 0.0454 | +0.0006 | 7/5/6 | 0.726 | - |
| (secondary) astrid-pro-min vs astrid-multi | 0.0482 | 0.0646 | -0.0164 | 10/7/1 | 0.0143 | - |
| (secondary) astrid-pro-min vs astrid-disco | 0.0482 | 0.0425 | +0.0057 | 4/7/7 | 0.108 | - |
| (secondary) astrid-pro-min vs astral-pro | 0.0482 | 0.0454 | +0.0028 | 7/3/8 | 0.442 | - |

Per-condition mean FN rate:

| condition | n | astrid-multi | astrid-pro | astrid-pro-w | astrid-pro-min | astrid-disco | astral-pro |
|---|---|---|---|---|---|---|---|
| default | 6 | 0.071 | 0.066 | 0.068 | 0.065 | 0.063 | 0.063 |
| gdl_1e-9_05 | 6 | 0.094 | 0.051 | 0.048 | 0.063 | 0.046 | 0.056 |
| gdl_1e-9_1 | 6 | 0.029 | 0.024 | 0.022 | 0.017 | 0.019 | 0.017 |

## Set B (new conditions, reps 01-05): n = 15 runs

| comparison | mean FN (Pro) | mean FN (base) | mean diff | W/T/L | p | Holm-adjusted p |
|---|---|---|---|---|---|---|
| astrid-pro vs astrid-multi | 0.0517 | 0.0517 | -0.0000 | 4/7/4 | 0.574 | 0.658 |
| astrid-pro vs astrid-disco | 0.0517 | 0.0483 | +0.0034 | 1/8/6 | 0.219 | 0.658 |
| astrid-pro vs astral-pro | 0.0517 | 0.0476 | +0.0041 | 3/5/7 | 0.325 | 0.658 |
| (secondary) astrid-pro-w vs astrid-multi | 0.0483 | 0.0517 | -0.0034 | 6/7/2 | 0.183 | - |
| (secondary) astrid-pro-w vs astrid-disco | 0.0483 | 0.0483 | -0.0000 | 2/8/5 | 0.865 | - |
| (secondary) astrid-pro-w vs astral-pro | 0.0483 | 0.0476 | +0.0007 | 4/6/5 | 0.856 | - |
| (secondary) astrid-pro-min vs astrid-multi | 0.0490 | 0.0517 | -0.0027 | 5/8/2 | 0.201 | - |
| (secondary) astrid-pro-min vs astrid-disco | 0.0490 | 0.0483 | +0.0007 | 3/7/5 | 0.725 | - |
| (secondary) astrid-pro-min vs astral-pro | 0.0490 | 0.0476 | +0.0014 | 4/7/4 | 0.618 | - |

Per-condition mean FN rate:

| condition | n | astrid-multi | astrid-pro | astrid-pro-w | astrid-pro-min | astrid-disco | astral-pro |
|---|---|---|---|---|---|---|---|
| gdl_5e-10_05 | 5 | 0.076 | 0.071 | 0.069 | 0.071 | 0.073 | 0.073 |
| ils_1e4 | 5 | 0.027 | 0.031 | 0.024 | 0.027 | 0.027 | 0.022 |
| ils_2e8 | 5 | 0.053 | 0.053 | 0.051 | 0.049 | 0.045 | 0.047 |

## Set A ∪ B (primary): n = 33 runs

| comparison | mean FN (Pro) | mean FN (base) | mean diff | W/T/L | p | Holm-adjusted p |
|---|---|---|---|---|---|---|
| astrid-pro vs astrid-multi | 0.0492 | 0.0588 | -0.0096 | 16/12/5 | 0.00702 | 0.0211 |
| astrid-pro vs astrid-disco | 0.0492 | 0.0451 | +0.0040 | 5/13/15 | 0.09 | 0.18 |
| astrid-pro vs astral-pro | 0.0492 | 0.0464 | +0.0028 | 8/11/14 | 0.302 | 0.302 |
| (secondary) astrid-pro-w vs astrid-multi | 0.0470 | 0.0588 | -0.0118 | 16/15/2 | 0.00146 | - |
| (secondary) astrid-pro-w vs astrid-disco | 0.0470 | 0.0451 | +0.0019 | 5/18/10 | 0.474 | - |
| (secondary) astrid-pro-w vs astral-pro | 0.0470 | 0.0464 | +0.0006 | 11/11/11 | 0.744 | - |
| (secondary) astrid-pro-min vs astrid-multi | 0.0485 | 0.0588 | -0.0102 | 15/15/3 | 0.00426 | - |
| (secondary) astrid-pro-min vs astrid-disco | 0.0485 | 0.0451 | +0.0034 | 7/14/12 | 0.157 | - |
| (secondary) astrid-pro-min vs astral-pro | 0.0485 | 0.0464 | +0.0022 | 11/10/12 | 0.368 | - |

Per-condition mean FN rate:

| condition | n | astrid-multi | astrid-pro | astrid-pro-w | astrid-pro-min | astrid-disco | astral-pro |
|---|---|---|---|---|---|---|---|
| default | 6 | 0.071 | 0.066 | 0.068 | 0.065 | 0.063 | 0.063 |
| gdl_1e-9_05 | 6 | 0.094 | 0.051 | 0.048 | 0.063 | 0.046 | 0.056 |
| gdl_1e-9_1 | 6 | 0.029 | 0.024 | 0.022 | 0.017 | 0.019 | 0.017 |
| gdl_5e-10_05 | 5 | 0.076 | 0.071 | 0.069 | 0.071 | 0.073 | 0.073 |
| ils_1e4 | 5 | 0.027 | 0.031 | 0.024 | 0.027 | 0.027 | 0.022 |
| ils_2e8 | 5 | 0.053 | 0.053 | 0.051 | 0.049 | 0.045 | 0.047 |
