# Pre-registration (written before any benchmark run of this round, 2026-10-10)

**Method under test (fixed by the theory, no tuning):** ASTRID-Pro = `apro -M pro` (MinDup rooting with
fewest-loss tie-break, species-overlap tagging, orthologous pairs only, speciation nodes only, gene-tree
root counted, per-gene mean then mean over genes) + FastME 2 (BalME, NNI+SPR). Secondary: ASTRID-Pro-S
(first pass = ASTRID-Pro, rooted by fewest total duplications; survival-reweighted second pass).

**Caveat on "held-out".** The prior pilot already ran an ASTRID-Pro variant on the FastMulRFS data and on
DISCO conditions default, gdl_1e-9_1, gdl_1e-9_05, gdl_5e-10_05, ils_1e4, ils_2e8. Nothing is tuned here;
the remaining DISCO conditions (gdl_1e-10_*, gdl_5e-10_0, gdl_1e-9_0, missing_1000, species_1000) are new.

**Metric.** Species-tree FN rate. Paired by (dataset, condition, replicate, sequence length, #genes).

**Primary comparisons** (two-sided Wilcoxon signed-rank on all simulated estimated-gene-tree runs pooled,
Holm over the four): ASTRID-Pro vs ASTRAL-Pro3, vs ASTRID-multi, vs ASTRID-DISCO, vs Asteroid.

**"Matches ASTRAL-Pro"** means: mean paired difference (ASTRID-Pro - ASTRAL-Pro3) <= +0.005 and the upper
end of its 95% bootstrap CI <= +0.010 FN rate.

**Go/kill.** Go if ASTRID-Pro matches ASTRAL-Pro3 (above) and is >= 10x faster at 1000 genes x 100 taxa.
Secondary (descriptive, no correction): DISCO+ASTRAL, FastMulRFS, wQFM-GDL, DupLoss-2, ASTRID-Pro-S;
per-condition tables; runtime scaling.
