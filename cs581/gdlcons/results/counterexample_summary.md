# Adversarial 4-taxon pools: fraction of datasets with a wrong species tree

True gene trees; disjoint blocks of K families from one pool per configuration.

## cand2

| method | 100 | 500 | 2000 | 10000 |
|---|---|---|---|---|
| ASTRAL-Pro, true tags | 24/40 | 4/40 | 0/20 | 0/4 |
| ASTRAL-Pro3 (own root+tags) | 38/40 | 40/40 | 20/20 | 4/4 |
| ASTRAL-Pro, overlap tags (true root) | 35/40 | 39/40 | 20/20 | 4/4 |
| ASTRAL-multi (astral4) | 26/40 | 24/40 | 13/20 | 2/4 |
| astrid-multi |  | 14/20 | 7/12 |  |
| astrid-disco |  | 20/20 | 12/12 |  |
| astral-disco |  | 20/20 | 12/12 |  |
| fastmulrfs |  | 10/20 | 9/12 |  |
| stag |  | 0/20 | 0/11 |  |

## cand1

| method | 100 | 500 | 2000 | 10000 |
|---|---|---|---|---|
| ASTRAL-Pro, true tags | 0/40 | 0/40 | 0/20 | 0/4 |
| ASTRAL-Pro3 (own root+tags) | 4/40 | 0/40 | 0/20 | 0/4 |
| ASTRAL-Pro, naive D->S q=0.02 | 11/40 | 2/40 | 0/20 | 0/4 |
| ASTRAL-Pro, naive D->S q=0.1 | 34/40 | 40/40 | 20/20 | 4/4 |

