# CS581 ML project scouting: validation, pilots, recommendation

Branch `claude/cs581-ml`. Everything here can be regenerated from `cs581/ml/code/`:

| file | what it does |
|---|---|
| `fetch.py` | streams the inputs from the MAGUS paper's `Results.zip` (IDB doi:10.13012/B2IDB-2643961_V1) |
| `validate_park.py` | rescores the published Park et al. 2021 trees (IDB doi:10.13012/B2IDB-7008049_V1) |
| `runtrees.py` | runs FastTree, IQ-TREE 3 and RAxML-NG, and scores each tree with `treeerr.py` (dendropy) |
| `colsupport.py`, `make_variants.py` | masking pilot |
| `frag.py` | fragment-aware pilot |
| `dnc.py` | GTM pipeline (written, not run) |
| `make_tables.py` | regenerates `results/tables.md` (every number below) |

The literature review is in `literature.md`.

Software and setup:
- FastTree 2.1.11 (double precision), IQ-TREE 3.1.4, RAxML-NG 2.0.3.
- Every job used 1 thread, with 4 jobs at a time on a 4-core VM.
- CPU time is the user+sys time of the tool process.
- FN is the missing-branch rate vs the true tree: the fraction of true internal edges absent from the estimate.

## 1. What was validated (our numbers vs. published ones)

**Park, Zaharias & Warnow 2021 (*Algorithms* 14:148).** These are the benchmarks behind the ML lecture's "RF error / runtime of ML heuristics" slides. The authors deposited inputs *and* output trees (IDB-7008049). Rescoring their trees with our scorer reproduces their tables:

| dataset (paper table) | FastTree | IQ-TREE 2 | RAxML-NG | GTM (IQ-TREE start) |
|---|---|---|---|---|
| RNASim1000, 5 reps (Table 3): paper | 14.9 | 15.1 | 15.1 | 14.4 |
| RNASim1000: published trees rescored by us | **14.9** | **15.1** | **15.1** | **14.4** |
| RNASim1000: our own fresh runs | **14.9** (FastTree) | 15.6 (IQ-TREE 3 `--fast`) | – | – |
| 1000M1-HF, 5 reps (Table 5): paper | 50.9 | 30.2 | 24.9 | 28.4 |
| 1000M1-HF: published trees rescored by us | 51.4 | **30.2** | **24.9** | **28.4** |
| 1000M1-HF: our own fresh runs | 48.9 | 37.5 (`--fast`) | see §2.1 | – |

Notes on this table:
- 1000M1-HF is ROSE 1000M1 with 500 of the 1000 sequences cut to fragments (about 25% of the median length).
- The input is not deposited as a single file. We rebuilt it from their per-subset alignments, which are full-width, and checked it against the ROSE originals: all 1000 sequences match, and 500 are identical.
- The published RAxML-NG trees are the `lastTree` of a run stopped at the 24 h cap (RAxML-NG 1.0.1, 20 starts, 2 threads).

**PASTA's own tree (`pasta.tre`, MAGUS paper data).** PASTA's final tree is FastTree `-gtr -gamma -fastest` on its final alignment (PASTA's default config). We re-ran that command on the published `pasta_align.txt`:

| condition | published `pasta.tre` FN (mean of R0–R2) | our FastTree `-fastest` FN |
|---|---|---|
| 1000M2 | 11.2 | 10.2 |
| 1000M3 | 7.2 | 7.2 |
| 1000L1 | 12.9 | 11.7 |
| RNASim 1K | 15.2 | 15.3 |

- Per-replicate differences are up to 2.5 points; RF between the two trees is 0.5–7.7%.
- So the setting is reproduced, but not bit-for-bit. Likely reasons: FastTree 2.1.10 vs 2.1.11, and PASTA seeding FastTree with its previous tree.
- Per-replicate table: `results/pasta_tre_check.md`.

**MAGUS-paper data, FastTree / IQ-TREE `--fast` on true / MAGUS / PASTA alignments (R0–R2):** see `results/tables.md`. The MAGUS paper reports alignment error only, so there is no published tree error to match. The relevant published claim is UPP 2015, Table 2: on these ROSE models, FastTree on PASTA alignments is about 1.3 points of FN worse than on the true alignment. We see the same small gap.

| FN, FastTree | 1000M2 | 1000M3 | 1000L1 | RNASim 1K |
|---|---|---|---|---|
| true alignment | 9.4 | 8.0 | 9.3 | 15.1 |
| MAGUS (`gcm.txt`) | 10.8 | 7.7 | 11.2 | 16.1 |
| PASTA | 10.2 | 7.2 | 12.0 | 15.0 |

Two takeaways for project choice:
1. The pipeline and data are validated against published numbers, so an apples-to-apples project is feasible.
2. On full-length ROSE/RNASim data the whole spread is small. The true alignment is only about 1 point better than MAGUS, and RAxML-NG is barely better than FastTree (§2.2). The big, unsolved gap is on **fragmentary data**: on 1000M1-HF, FastTree gets 49–51% FN, while RAxML-NG reaches about 22–25%.

## 2. Pilots

PILOT_TABLES

## 3. Recommendation

RECOMMENDATION
