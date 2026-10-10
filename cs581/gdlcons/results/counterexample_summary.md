# Adversarial 4-taxon pools: fraction of datasets with a wrong species tree

True gene trees; disjoint blocks of K families from one pool per configuration.

## cand2 (hidden-paralog tagging)

| method | 100 | 500 | 1000 | 2000 | 5000 | 10000 |
|---|---|---|---|---|---|---|
| ASTRAL-Pro, true tags | 24/40 | 4/40 |  | 0/20 |  | 0/4 |
| ASTRAL-Pro3 (own root+tags) | 38/40 | 40/40 |  | 20/20 |  | 4/4 |
| ASTRAL-Pro, overlap tags (true root) | 35/40 | 39/40 |  | 20/20 |  | 4/4 |
| ASTRAL-multi (astral4) | 26/40 | 24/40 |  | 13/20 |  | 2/4 |
| astrid-multi |  | 14/20 |  | 14/20 |  | 3/4 |
| astrid-disco |  | 20/20 |  | 20/20 |  | 4/4 |
| astral-disco |  | 20/20 |  | 20/20 |  | 4/4 |
| fastmulrfs |  | 10/20 |  | 12/20 |  | 3/4 |
| stag |  | 0/20 |  | 0/20 |  | 0/4 |

## cand3 (random re-rooting)

| method | 100 | 500 | 1000 | 2000 | 5000 | 10000 |
|---|---|---|---|---|---|---|
| ASTRAL-Pro, true tags |  |  |  | 0/20 |  | 0/4 |
| ASTRAL-Pro3 (own root+tags) |  |  |  | 6/20 |  | 0/4 |
| ASTRAL-Pro, overlap tags (true root) |  |  |  | 6/20 |  | 0/4 |
| ASTRAL-Pro, random re-root p=0.3 |  |  |  | 13/20 |  | 4/4 |
| ASTRAL-Pro, random re-root p=1 |  |  |  | 20/20 |  | 4/4 |
| ASTRAL-Pro, root on an A leaf |  |  |  | 20/20 |  | 4/4 |

## cand4 (ASTRAL-Pro's own rooting)

| method | 100 | 500 | 1000 | 2000 | 5000 | 10000 |
|---|---|---|---|---|---|---|
| ASTRAL-Pro, true tags |  |  | 0/20 |  | 0/4 |  |
| ASTRAL-Pro3 (own root+tags) |  |  | 20/20 |  | 4/4 |  |
| ASTRAL-Pro, overlap tags (true root) |  |  | 3/20 |  | 0/4 |  |
| ASTRAL-Pro, random re-root p=0.3 |  |  | 5/20 |  | 0/4 |  |
| ASTRAL-Pro, root on an A leaf |  |  | 20/20 |  | 4/4 |  |
| astrid-multi |  |  | 1/20 |  | 0/4 |  |
| astrid-disco |  |  | 20/20 |  | 4/4 |  |
| astral-disco |  |  | 20/20 |  | 4/4 |  |
| fastmulrfs |  |  | 0/20 |  | 0/4 |  |
| stag |  |  | 20/20 |  | 4/4 |  |

## cand1 (only naive flips fail)

| method | 100 | 500 | 1000 | 2000 | 5000 | 10000 |
|---|---|---|---|---|---|---|
| ASTRAL-Pro, true tags | 0/40 | 0/40 |  | 0/20 |  | 0/4 |
| ASTRAL-Pro3 (own root+tags) | 4/40 | 0/40 |  | 0/20 |  | 0/4 |
| ASTRAL-Pro, naive D->S q=0.02 | 11/40 | 2/40 |  | 0/20 |  | 0/4 |
| ASTRAL-Pro, naive D->S q=0.1 | 34/40 | 40/40 |  | 20/20 |  | 4/4 |
| astrid-multi |  | 15/20 |  | 11/18 |  |  |
| astrid-disco |  | 1/20 |  | 0/18 |  |  |
| astral-disco |  | 1/20 |  | 0/18 |  |  |
| fastmulrfs |  | 18/20 |  | 18/18 |  |  |
| stag |  | 0/20 |  | 0/17 |  |  |

