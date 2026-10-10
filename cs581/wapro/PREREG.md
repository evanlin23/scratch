# Pre-registration: weighted ASTRID-Pro (written and pushed before any benchmark result was looked at)

*2026-10-10. Branch `claude/cs581-wapro`.*

**Disclosure.** Before writing this, I ran every method once on a single smoke-test input to check the plumbing:
FastMulRFS dl 1e-10, Ne 1e7, rep 01, 25 bp. That input is in TRAIN. No variant was chosen from it.

## Gene trees (identical for every method)

The published gene trees have no branch support, so I re-estimate them from the published alignments. Tool:
**FastTree 2.1 `-nt -gtr -gamma`**, which gives SH-like local support in [0, 1] on internal edges, plus branch lengths.

- **FastMulRFS data** (Molloy & Warnow 2020, the ASTRAL-Pro S100 data).
  - 100 species. DL rate 1e-10, 2e-10 or 5e-10; Ne 1e7 or 5e7; reps 01–10.
  - The first 100 genes of `genes-with-gt3species.txt`. These are the same genes as the published 100-gene RAxML sets.
  - Each gene is truncated to its first 25 or 100 sites, as in the paper.
  - 120 datasets in total.
- **DISCO data** (Willson et al. 2022).
  - 100 species, all 1000 families, 100 bp alignments.
  - Conditions: `default`, `gdl_1e-10_05`, `gdl_5e-10_05` and `ils_2e8`; reps 01–04.
  - 16 datasets in total.

## Split

| set | FastMulRFS | DISCO | datasets |
|---|---|---|---|
| **TRAIN** (choosing the variant, and the baselines' settings) | reps 01–03, all 6 conditions × {25, 100} bp | `default`, `gdl_1e-10_05`, reps 01–02 | 36 + 4 = 40 |
| **HELD-OUT** (evaluated once, with the frozen choice) | reps 04–10, all conditions × {25, 100} bp | `default`, `gdl_1e-10_05` reps 03–04; `gdl_5e-10_05`, `ils_2e8` reps 01–04 | 84 + 12 = 96 |

If compute runs short, I drop DISCO held-out datasets in this order: `ils_2e8` reps 04, 03, …, then `gdl_5e-10_05`
reps 04, 03, …. I report every drop.

## Candidate variants (TRAIN only)

All variants use ASTRID-Pro with its usual steps: MinDup rooting, species-overlap tags, orthologous pairs only,
gene-tree root counted, FastME 2 (BalME + NNI + SPR). Flags are those of `code/wapro`. Support s is FastTree's
SH-like value, mapped to [0, 1].

| id | flags | idea |
|---|---|---|
| U | (none) | unweighted ASTRID-Pro (= previous round) |
| S | `-W sup` | each counted speciation node weighs the support of the path edge below it; the LCA weighs the mean of its two path edges |
| SE | `-W supe` | edge form: sum of support over path edges that enter a speciation node. Equals wASTRID-s on single-copy trees |
| S2 | `-W sup -p 2` | support squared |
| L | `-W len -N diam` | wASTRID-pl analogue: patristic length over orthologous pairs, per gene divided by the tree diameter |
| C1, C3, C5 | `-c 0.1`, `-c 0.3`, `-c 0.5` | lever (a): contract edges with support < τ before rooting and tagging |
| N | `-N mean` | lever (b): divide each gene's distances by its mean (scale normalisation for missing species / gene size) |
| P | `-A pair` | lever (c): pool orthologous pairs over genes (per-pair) instead of per-gene averaging |
| G | `-G msup` | gene weight = mean internal support |

**Combination round (still TRAIN).**
- Take the best weighting among {S, SE, S2, L}, or U if none beats U.
- Add each lever from {C1/C3/C5 (best τ), N, P, G} that improved mean TRAIN FN when used alone. Add them greedily,
  in order of their solo gain, and keep a lever only if it lowers mean TRAIN FN.

**Selection rule.** Lowest mean TRAIN FN rate. If two variants are within 0.001 of each other, pick the one with
fewer flags. The chosen variant is called **wASTRID-Pro**. I write it into `results/CHOSEN.md` and commit that file
*before* running the held-out set.

## Baselines (best recommended setting; run on TRAIN and HELD-OUT)

- **ASTRAL-Pro3** (ASTER v1.25.3.8). It has no weighting option for multi-copy input; I checked the source and the
  ASTER 2025 paper.
- **DISCO + wASTRAL**: DISCO v1.4.1, then `wastral --mode 1` (hybrid weighting; support range auto-detected).
- **DISCO + wASTRID-s**: internode v0.0.7, `-m support -b 0-1`.
- **ASTRID-DISCO**: DISCO, then ASTRID (internode `-m internode`).
- **Asteroid**: default setting (missing-data correction on), and `--use-gene-bl`. Its best setting on TRAIN (lowest
  mean FN) is the one used on HELD-OUT.
- **wQFM-GDL** v1.0.3, tree mode (`-t`). It uses no support.
- **unweighted ASTRID-Pro** (U).

## Hypotheses and tests (HELD-OUT only)

- **Metric.** Species-tree FN rate. Paired by (dataset, condition, replicate, sequence length).
- **Primary.**
  - Comparison: wASTRID-Pro vs **the best baseline**, i.e. the baseline with the lowest mean HELD-OUT FN rate.
    Choosing it this way is conservative: it works against us.
  - Test: two-sided Wilcoxon signed-rank (zero differences dropped), α = 0.05.
  - Reported with the mean difference, a 95% bootstrap CI, and W/T/L.
- **Secondary.** wASTRID-Pro vs each baseline, with Holm correction over the 7 baselines. Also per-dataset tables,
  and runtime.
- **Runtime.**
  - Wall-clock and CPU time (sum over child processes), with 1 thread on every held-out dataset.
  - With 4 threads (where the tool supports it) on a subset: all DISCO held-out datasets, plus FastMulRFS rep 04.

**Verdict rule** (beating the best existing methods):
- **promising**: primary p < 0.05 in favour of wASTRID-Pro AND mean improvement ≥ 0.005 FN rate.
- **unclear**: improvement in the right direction but not significant, or significant but < 0.005.
- **not promising**: otherwise.
