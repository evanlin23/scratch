# Weighted ASTRID-Pro (wASTRID-Pro): pilot report

*CS581 pilot, 2026-10-10/11. Branch `claude/cs581-wapro`. The pilot was stopped early at the orchestrator's request,
because the student chose a different project. **Incomplete parts are marked ⚠.***

- Code: `code/`. Pre-registration: `PREREG.md`, pushed before any result.
- Frozen variant choice: `results/CHOSEN.md`. Interim look: `results/INTERIM.md`.
- Full tables: `results/tables.md`. Raw runs: `results/*.jsonl`.

## 0. Verdict

**Unclear, leaning not promising, for "faster AND more accurate than the best existing methods".**

**What weighting buys.**
- wASTRID-Pro (support-weighted ASTRID-Pro) is significantly but only slightly more accurate than unweighted
  ASTRID-Pro on held-out data: −0.0025 FN rate, 38/36/22, p = 0.018.
- It is also slightly ahead of ASTRID-DISCO (−0.0030, p = 0.019, Holm 0.072).
- These gains are below the pre-registered 0.005 margin, so by the pre-registered rule the verdict is "unclear".

**Where it stands against the other baselines.**
- It ties Asteroid, DISCO+wASTRID and wQFM-GDL.
- It clearly beats ASTRAL-Pro3 and DISCO+wASTRAL (−0.017 to −0.018, p < 1e-5).

**Speed.** It keeps ASTRID-Pro's speed: 0.1 s CPU per 100 genes and about 1 s per 1000 genes. That is 10× faster than
ASTRID-DISCO / Asteroid and 300–2000× faster than ASTRAL-Pro3 / wQFM-GDL.

**Warning sign on the 1000-gene DISCO data.**
- wQFM-GDL is far more accurate there: 0.029 vs 0.054, on only 6 paired datasets (⚠ incomplete).
- DISCO+wASTRID is also better there: +0.012 against us, p = 0.06.
- So "more accurate than the best method" is not supported where gene trees are better (100 bp, 1000 genes).
- The gain is concentrated in the noisy FastMulRFS gene trees (25 bp, 100 genes).

## 1. Prior art and novelty (`prior_art.md`)

- **No GDL-aware weighted ASTRID exists.**
  - wASTRID (Liu & Warnow 2023) lists GDL / DISCO+wASTRID as future work.
  - Its code says ASTRID-multi is not implemented.
- **ASTRAL-Pro3 (ASTER v1.25.3.8) has no weighting option.**
  - The ASTER 2025 paper (Zhang, Nielsen & Mirarab, MBE 42:msaf172) says weighted multi-copy quartet scoring is
    unsolved.
  - wASTRAL (`wastral --mode 1`) is single-copy only.
- **Other methods do not weight by support.**
  - wQFM-GDL uses no support or length.
  - Asteroid has only `--use-gene-bl` (raw lengths, no tagging) and length-threshold contraction.
- **DISCO + wASTRAL and DISCO + wASTRID** appear in no benchmark I could find. Both are run here as baselines.
- **Conclusion:** support/length weighting on orthologous-pair, speciation-node distances is new. It is a small
  increment on ASTRID-Pro + wASTRID-s.

## 2. Methods

**Gene trees.** The published gene trees carry no support, so all methods use the same re-estimated gene trees.
- Tool: FastTree 2 `-nt -gtr -gamma`, which gives SH-like local support plus branch lengths.
- **FastMulRFS / ASTRAL-Pro S100 data.**
  - 100 species; DL rate 1e-10, 2e-10 or 5e-10 × Ne 1e7 or 5e7; reps 01–10.
  - First 100 genes, first 25 or 100 sites.
  - 120 datasets.
- **DISCO data.**
  - 1000 families, 100 bp; conditions `default`, `gdl_1e-10_05`, `gdl_5e-10_05`, `ils_2e8`; reps 01–04.
  - 16 datasets.
- Scripts: `code/estimate_fmrfs.py` and `code/estimate_disco.py`.

**wASTRID-Pro** (`code/wapro.cpp`, an extension of the round-2 `apro.cpp`).
- **Unchanged from ASTRID-Pro:**
  - MinDup rooting.
  - Species-overlap tagging.
  - Orthologous pairs only.
  - Speciation nodes on the path, gene-tree root counted.
  - Per-gene mean, then mean over genes.
  - FastME 2.
- **New:** lengths and supports are carried through rooting. The root edge is split: half its length to each side,
  and the same support on both halves.
- **New:** a weight per counted node. Variants:
  - `-W sup` (chosen): each counted speciation node weighs the support of the path edge below it. The LCA weighs the
    mean of its two path edges.
  - `-W supe`: edge sum. On single-copy trees this equals wASTRID-s; verified against a brute-force computation,
    max difference 4e-15.
  - `-W len -N diam`: wASTRID-pl analogue.
  - Contraction `-c τ`.
  - Per-gene scale normalisation `-N mean`.
  - Per-pair pooling `-A pair`.
  - Gene weights `-G msup`.
- **Check:** with unit weights the output matrix is byte-identical to `apro`.

**Baselines** (`code/bench.py`; versions and flags):
- ASTRAL-Pro3 v1.25.3.8, `-u 0`. No weighting option exists.
- DISCO v1.4.1 + wASTRAL v1.25.3.8, `--mode 1` (hybrid).
- DISCO + wASTRID-s: internode v0.0.7, `-m support -b 0-1`. Two compatibility patches were needed to build it:
  bindgen 0.69, and one `usize` cast.
- ASTRID-DISCO: internode `-m internode`.
- Asteroid at GitHub HEAD 496d1b4: default and `--use-gene-bl`. `-bl` was best on TRAIN, so it is the designated
  Asteroid baseline.
- wQFM-GDL v1.0.3, `-t`.
- Unweighted ASTRID-Pro.

**Protocol** (`PREREG.md`).
- **TRAIN:** FastMulRFS reps 01–03 (36 datasets) + DISCO `default` / `gdl_1e-10_05` reps 01–02 (4 datasets).
- **HELD-OUT:** FastMulRFS reps 04–10 (84 datasets) + 12 DISCO datasets.
- **Metric:** FN rate. Two-sided Wilcoxon, bootstrap CI, Holm over baselines.

## 3. TRAIN: variant selection (n = 40)

| variant | flags | mean FN |
|---|---|---|
| **S (chosen)** | `-W sup` | **0.0706** |
| U (unweighted) | – | 0.0716 |
| G | `-G msup` | 0.0714 |
| S+G | `-W sup -G msup` | 0.0714 |
| S2 | `-W sup -p 2` | 0.0716 |
| SE | `-W supe` | 0.0719 |
| C5 / C3 / C1 | `-c 0.5` / `0.3` / `0.1` | 0.0721 / 0.0745 / 0.0745 |
| N | `-N mean` | 0.0734 |
| P | `-A pair` | 0.0783 |
| L | `-W len -N diam` | 0.0801 |

**What each lever did:**
- Support weighting helps a little.
- Contraction, normalisation, per-pair pooling and length weighting all hurt.

**Baselines on TRAIN:**

| baseline | mean FN |
|---|---|
| ASTRID-DISCO | 0.0706 |
| DISCO+wASTRID | 0.0706 |
| Asteroid-bl | 0.0721 |
| wQFM-GDL | 0.0729 |
| Asteroid | 0.0732 |
| ASTRAL-Pro3 | 0.0853 |
| DISCO+wASTRAL | 0.0888 |

## 4. HELD-OUT (evaluated once with the frozen choice)

wASTRID-Pro minus each method. Negative means wASTRID-Pro is better.

| vs | n | mean diff | 95% CI | W/T/L | p | Holm p |
|---|---|---|---|---|---|---|
| **ASTRID-Pro (best baseline by held-out mean, 0.0739): PRIMARY** | 96 | **−0.0025** | [−0.0047, −0.0002] | 38/36/22 | **0.018** | 0.072 |
| ASTRID-DISCO | 96 | −0.0030 | [−0.0059, 0.0000] | 37/29/30 | 0.019 | 0.072 |
| Asteroid (default) | 96 | −0.0034 | [−0.0075, +0.0005] | 40/23/33 | 0.13 | – |
| Asteroid-bl (designated) | 96 | −0.0073 | [−0.0129, −0.0017] | 43/21/32 | 0.004 | 0.021 |
| DISCO+wASTRID | 96 | −0.0035 | [−0.0077, +0.0006] | 38/30/28 | 0.11 | 0.22 |
| wQFM-GDL ⚠ | 90 | −0.0045 | [−0.0107, +0.0017] | 40/15/35 | 0.41 | 0.41 |
| ASTRAL-Pro3 ⚠ | 90 | −0.0171 | [−0.0241, −0.0102] | 56/13/21 | 7e-6 | 4e-5 |
| DISCO+wASTRAL ⚠ | 90 | −0.0183 | [−0.0253, −0.0119] | 54/15/21 | 6e-7 | 4e-6 |

**Mean held-out FN rate:**

| method | all | FastMulRFS 25 bp | FastMulRFS 100 bp | DISCO |
|---|---|---|---|---|
| wASTRID-Pro | 0.0715 | 0.105 | 0.043 | 0.054 |
| ASTRID-Pro | 0.0739 | 0.108 | 0.045 | 0.055 |
| ASTRID-DISCO | 0.0745 | 0.111 | 0.044 | 0.054 |
| Asteroid | 0.0749 | 0.112 | 0.044 | 0.053 |
| DISCO+wASTRID | 0.0749 | 0.116 | 0.043 | **0.042** |
| wQFM-GDL | 0.0760 | 0.117 | **0.042** | **0.029** (6 datasets) |
| ASTRAL-Pro3 | 0.0886 | 0.136 | 0.048 | 0.048 (6) |

**Subsets.**
- **Excluding the 9 interim-seen datasets:**
  - vs ASTRID-Pro −0.0031, p = 0.010.
  - vs ASTRID-DISCO −0.0034, p = 0.013.
  - vs wQFM-GDL −0.0066, p = 0.16.
- **DISCO only (n = 12):**
  - vs ASTRID-Pro −0.0017, p = 0.62.
  - vs DISCO+wASTRID **+0.0119**, p = 0.06.
  - vs wQFM-GDL **+0.0085** (n = 6, 1/1/4).
- **FastMulRFS only (n = 84):**
  - vs ASTRID-Pro −0.0026, p = 0.026.
  - vs DISCO+wASTRID −0.0056, p = 0.027.
  - vs wQFM-GDL −0.0054, p = 0.31.

## 5. Runtime (mean seconds, wall / CPU, 1 thread; TRAIN and HELD-OUT pooled)

| method | FastMulRFS (100 genes) | DISCO (1000 genes) |
|---|---|---|
| wASTRID-Pro | 0.25 / 0.08 | 2.1 / 0.9 |
| ASTRID-Pro | 0.15 / 0.08 | 2.1 / 0.9 |
| ASTRID-DISCO | 0.51 / 0.38 | 16 / 8.0 |
| DISCO+wASTRID | 0.51 / 0.38 | 18 / 8.1 |
| Asteroid | 1.4 / 1.1 | 20 / 9.6 |
| ASTRAL-Pro3 | 7.2 / 6.1 | 445 / 287 |
| DISCO+wASTRAL | 12 / 10 | 727 / 446 |
| wQFM-GDL | 52 / 47 | 1658 / 1865 (multi-threaded Java) |

⚠ **Runtime caveats.**
- These were measured on a shared, overloaded 4-core machine: 5–7 jobs on 4 cores, so wall time is often more than
  CPU time.
- The pre-registered clean 1-thread vs 4-thread timings were **not run**.
- The CPU ratios are large enough to stand anyway. Weighting adds no measurable cost.

## 6. What is incomplete ⚠

- **wQFM-GDL, ASTRAL-Pro3 and DISCO+wASTRAL on DISCO held-out.**
  - They finished on 6 of 12 datasets: `default` 03–04, `gdl_1e-10_05` 03–04, and `gdl_5e-10_05` 01–02, in the
    pre-registered drop order.
  - None ran on `ils_2e8`.
- **TRAIN gaps.** wQFM-GDL is missing on 2 train DISCO datasets, and DISCO+wASTRAL on 1. These affect only the
  descriptive train baseline numbers.
- **Interim look.** The orchestrator asked for one interim held-out look at 9 FastMulRFS rep-04 datasets
  (`results/INTERIM.md`). The conclusion holds with and without them.
- **Not done:**
  - IQ-TREE aBayes gene trees. wASTRID found aBayes best; it is possible that weighting gains more with better
    support values.
  - Larger DISCO replication.
  - 4-thread timings.

## 7. If continued: weeks 1–4, and risks

- **Week 1: aBayes support.** IQ-TREE `-abayes -fast` costs about 0.75 s per gene. Re-run TRAIN with it: is the
  gain from weighting larger with aBayes than with SH-like support?
- **Week 2: close the DISCO gap.** Finish wQFM-GDL / ASTRAL-Pro3 on all DISCO held-out datasets, then work out why
  DISCO+wASTRID beats both ASTRID-Pro variants on the 1000-gene data:
  - support weighting after DISCO's orthology split, vs
  - our in-tree tagging.
  A hybrid could weight the DISCO-split trees.
- **Week 3:** theory. With exact supports (s = 1 on true edges), the β-consistency theorem from round 2 carries over
  unchanged. Bound the bias added by support < 1 on true edges.
- **Week 4:** write-up.

**Risks.**
1. **The effect size is tiny** (0.0025 FN rate, ≈ 0.25 edges per 97). That is below the pre-registered meaningful
   margin.
2. **The quartet methods may win on better gene trees.** wQFM-GDL and wASTRAL-family methods appear stronger once
   gene trees are better (DISCO, 1000 genes). A distance method may be fundamentally capped there.
3. **The novelty is incremental:** ASTRID-Pro + wASTRID-s.
