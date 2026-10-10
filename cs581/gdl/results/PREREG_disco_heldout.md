# Pre-registration: DISCO-data held-out test (written 2026-10-10 ~01:55 UTC, before inspecting held-out scores)

**Status of the data when this was written.**
- The runs for replicates 05–10 of `default`, `gdl_1e-9_1` and `gdl_1e-9_05` had already been computed by the background workers. Their scores had **not been looked at**; only their (condition, replicate) keys were listed.
- I had looked at replicates 01–04. These are the **dev** set, and they motivated the variant choice below.

**Fixed in advance.**
- **Primary method.** `astrid-pro` = mean over orthologous pairs, speciation nodes only, inferred root and tags, root not counted, unweighted.
- **Primary comparisons.** `astrid-pro` vs `astrid-multi`, `astrid-disco` and `astral-pro`.
- **Metric.** Species-tree FN rate. Paired by (condition, replicate).
- **Test.** Two-sided Wilcoxon signed-rank (zero differences dropped). Holm correction across the 3 primary comparisons at α = 0.05. Report the mean difference and W/T/L, where tie = equal FN count.
- **Held-out set A.** Replicates 05–10 of `default` (1000 genes), `gdl_1e-9_1` (1000 genes) and `gdl_1e-9_05` (first 100 genes, because families average about 1000 leaves). That is 18 runs.
- **Held-out set B.** New conditions that have not been run yet, replicates 01–05, 1000 genes from 100 bp:
  - `gdl_5e-10_05` (loss/dup 0.5, supercritical);
  - `ils_1e4` (very low ILS);
  - `ils_2e8` (high ILS).
  That is 15 runs. Its scores are reported pooled with A, and also separately.
- **Secondary, no correction.** `astrid-pro-w` (support-weighted) and `astrid-pro-min` against the same baselines.
- **Success criterion.** "Confirmed" if `astrid-pro` beats `astrid-multi` after Holm correction on A ∪ B. Otherwise the DISCO-data gain stays exploratory.
