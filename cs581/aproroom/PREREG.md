# Pre-registration: where do GDL species-tree methods still have headroom? (aproroom)

*Written 2026-10-10, before any experiment below was run (only a smoke test on one FastMulRFS replicate of true
gene trees, where every method had FN = 0). Analysis code: `code/analyze.py`.*

## Data
- **FastMulRFS / ASTRAL-Pro S100 data**, Molloy & Warnow 2020, Illinois Data Bank doi:10.13012/B2IDB-5721322_V1.
  100 species; DL rate 1e-10, 2e-10, 5e-10 × Ne 1e7, 5e7; 10 replicates. True gene trees (SimPhy), true locus
  trees, RAxML gene trees from 25 and 100 bp. First 100 genes of `genes-with-gt3species.txt`.
- **DISCO data**, Willson et al. 2022, Illinois Data Bank doi:10.13012/B2IDB-4050038_V1. 100 species (default,
  gdl_*, ils_*, missing_1000, gtrees_10000_l1) and 1000 species (species_1000). True gene trees + locus trees, and
  gene trees estimated from 50/100/500 bp.

## Levels of the error decomposition (Q1)
- **L0** true gene tree, **true root and true D/S tags** (from SimPhy locus trees, `code/truetag.py`).
  Only tag-aware pipelines can use L0: ASTRID-Pro-tt, ASTRID-DISCO-tt, DISCO-ASTRAL-tt (DISCO decomposition at
  the true duplications). ASTRAL-Pro3, Asteroid and wQFM-GDL do their own rooting/tagging (or none) and have no
  L0 input.
- **L1** true gene trees, estimated root/tags (each method's own).
- **L2** estimated gene trees (FastMulRFS: RAxML 100 bp and 25 bp; DISCO: 100 bp, plus 50/500 bp on default).

Components per method (mean FN rate, paired by replicate):
- method + finite-sample error = L0 (or L1 for methods without an L0 input);
- tagging/rooting error = L1 − L0;
- gene-tree-estimation error (GTEE) = L2 − L1.

## Main comparisons (decided now)
1. **Decomposition.** For each method, L2 − L1 and L1 − L0 (Wilcoxon, two-sided). Hypothesis: GTEE dominates
   (> 70% of L2 error) everywhere; tagging error is < 0.005 FN on average.
2. **Between-method gap per level.** For each level, ASTRID-Pro vs each of ASTRID-DISCO, Asteroid, ASTRAL-Pro3,
   wQFM-GDL: mean diff, W/T/L, Wilcoxon, Holm over the 4 within each level. Hypothesis: gaps are ~0 at L1 and
   largest at L2 with short sequences (25 bp).
3. **Headroom.** Per replicate, the oracle best method (minimum FN over the 5 main methods) vs the best single
   method; and the L1 error of the best method as a floor for what better handling of GTEE could reach.
   This is the bound reported as "how much a new method can gain at most" on these data.
4. **Regimes (Q2).** On DISCO estimated trees: per condition (high-GDL, high-ILS, missing data, 1000 species,
   few vs many genes), ASTRID-Pro vs each method, W/T/L and Wilcoxon (uncorrected; reported as exploratory where
   n < 6).
5. **Hybrid (Q4).** ASTRAL-Pro3 seeded with the ASTRID-Pro tree as a guide (`-g`), with reduced search
   (`-r 1 -s 0`, "apro3-guide-fast") and with default search ("apro3-guide"), vs default ASTRAL-Pro3 and vs
   ASTRID-Pro. Control: reduced search without guide ("apro3-fast"). Outcome: FN rate and runtime. Success if
   apro3-guide-fast is not worse than ASTRAL-Pro3 (mean diff ≤ +0.002, p > 0.05) at ≥ 3× lower runtime.
6. **Speed/memory (Q3).** Wall time and peak RSS at 1 and 4 threads; scaling in species (50 → 1000) and genes
   (100 → 10,000). Descriptive only.

Statistics: FN rate (= RF rate for binary trees), paired by (dataset, condition, replicate, level), two-sided
Wilcoxon signed-rank (zero differences dropped, Pratt not used), W/T/L. Failed runs (crash/timeout) are reported
and excluded pairwise.
