# Pre-registration: MAGUS / PASTA with different base (subset) aligners

Written 2026-10-10, before any base-method result was computed.

## Question

MAGUS aligns its 25 subsets with MAFFT L-INS-i and merges them with GCM. Does another subset aligner
(MUSCLE5, FAMSA2, ProbCons, Clustal Omega; Prank on nucleotides) give a more accurate MAGUS alignment?
Instructor projects (a) "MAGUS with different base methods", (b) "Prank/ProbCons within PASTA or MAGUS",
(c) "Regressive vs PASTA and MAGUS".

## Design (paired, merge-only)

For each dataset, one MAGUS run (paper flags: 25 subsets, PASTA-style decomposition, 10 backbones x 200
sequences aligned by L-INS-i, MCL, minclusters) provides the subsets and the backbones. Each subset is
re-aligned from its unaligned sequences by each base method; then the same GCM merge
(`gcmx.run_magus -s SUBSETS -b BACKBONES`, same backbones) is run. Only the subset aligner changes.
Control: `merge-linsi` = the subsets re-aligned with MAGUS's own L-INS-i command (same as MAGUS).

## Primary analysis

- Data: all protein datasets that complete for all methods (BAliBASE RV100 length-filtered sets,
  AliSim SIMMOD/SIMHIGH replicates, HomFam families subsampled to 2,000 sequences; one MAGUS draw each).
- Metric: SP error = (SPFN + SPFP) / 2 of the final merged alignment (FastSP; HomFam scored on the seed
  sequences only).
- Candidates (proteins): MUSCLE5 (`-align`), FAMSA2 (default), ProbCons (default), Clustal Omega (default).
- Test: for each candidate, paired difference (candidate minus L-INS-i) over the protein datasets,
  two-sided Wilcoxon signed-rank; Holm correction over the 4 candidates.
- "Best base method" = the candidate with the lowest mean SP error over the protein datasets.
  **Primary claim**: best candidate vs L-INS-i inside MAGUS, paired Wilcoxon, Holm-adjusted p < 0.05 and
  mean difference < 0. Because "best" is picked on the same data, the Holm adjustment over all four is
  what is reported for it.
- W/T/L with a tie band of |d| < 0.05 SP-error points (percentage points).

## Secondary (descriptive, no multiplicity control)

- SPFN, SPFP, TC; per data family (BAliBASE / simulated / HomFam); nucleotide sets (ROSE 1000M2, 1000L1,
  RNASim 1000, 16S.M) with Prank added, MUSCLE5, Clustal Omega, FAMSA (nucleotide mode as available).
- Subset-level accuracy (each method's subset alignments scored against the reference restricted to the
  subset), to explain the merge results.
- Runtime: subset-alignment wall and CPU time per method (1 thread per subset, 4 subsets in parallel),
  merge wall time.
- Trees: FastTree `-lg -gamma` on the simulated protein sets for MAGUS (L-INS-i subsets) and the best variant;
  nRF vs the true tree.
- PASTA 1.8.3 with its built-in aligner choices on a few datasets.
- Whole-data baselines: MUSCLE5, FAMSA2, regressive T-Coffee, where they finish in ~30 min.

## Decision rule for the verdict

"Promising" if the primary test succeeds, or if a candidate is at least as accurate as L-INS-i (no
significant loss) while much faster, giving a speed/accuracy story. "Not promising" if every candidate is
worse or equal at no speed gain. "Unclear" otherwise.
