# Pre-registration: replacing MCL in MAGUS's GCM

Written 2026-10-10 ~22:55 UTC, before any grid result. Already seen at that point: MAGUS control merges
(reproduction checks) on 1000M2, 1000L1, 1000L2, BBA0101, BBA0134; a debugging probe on the training set
BBA0101 only (`results/probe_BBA0101.jsonl`: raw:mcl:4, es4:mcl:4 with a slightly different filter,
raw:leidmod:20, raw:cc:4, raw:agglo:1, run with Leiden iterated to convergence); and cluster statistics (no
accuracy) of Leiden-CPM / Louvain / Leiden-modularity / LPA on BBA0101 to place the resolution grids.

## Fixed inputs

Merge-only. Each replicate = one cached (or fresh, for AliSim) MAGUS draw: its 25 L-INS-i subset alignments and
10 L-INS-i backbones (200 sequences) are fixed; only the clustering step of GCM changes. MAGUS's own
`purgeDuplicateClusters` + `purgeClusterViolations` + `minclusters` trace and the alignment writer are unchanged
(`--graphtraceoptimize false`, i.e. the paper's flags). Baseline = `raw:mcl:4` (= MAGUS; reproduces the
cached score exactly).

Graphs: `raw` (MAGUS's) and `es4` (every non-self edge contributed by fewer than 4 of the 10 backbones deleted;
gcmgen's `linsi#es4`; the trace also sees the filtered graph).

## Split (gcmgen's, `cs581/gcmgen/SPLIT.md`)

- **Training (hyperparameters chosen here):** BBA0101, BBA0134, BBA0067, BBA0039, SIMMOD_R1, SIMHIGH_R1 (protein);
  1000M2, 1000L1, 1000L2, 16S.M (DNA/RNA).
- **Held out (evaluated once):** BBA0154, BBA0190, SIMMOD_R2, SIMHIGH_R2 (protein); 1000S1, 1000L3, 1000M3,
  1000M4, 1000S2, 1000S3, RNASim (R0), RNASim_R1 (DNA/RNA).

## Training grid (31 variants + baseline, every training replicate)

| family | raw graph | es4 graph |
|---|---|---|
| MCL inflation I | 2, 3, 6 (4 = baseline) | 2, 3, 4, 6 |
| Leiden, modularity (RBConfiguration), resolution γ | 10, 30, 100 | 10, 30, 100 |
| Leiden, CPM on w/sqrt(s_a s_b), γ | 0.003, 0.01, 0.02 | 0.003, 0.01, 0.02 |
| Louvain (igraph multilevel, modularity), γ | 30 | 30 |
| connected components of cross-subset edges with support ≥ t | 6, 8, 10 | — (identical to raw for t ≥ 4) |
| greedy constrained agglomeration, support ≥ t | 1, 2, 4 | 1 |
| label propagation | lpa | lpa |

Leiden/Louvain/LPA use seed 1 (leidenalg default 2 iterations). Merge wall limit 30 min; a variant that times out
or fails on any training replicate is not eligible.

## Selection rule (training only)

Δ = variant error − MAGUS error, error = (SPFN + SPFP)/2 × 100 (FastSP). For each family × graph, the parameter
with the lowest mean Δ over the 10 training replicates (each replicate weighted equally) is selected
(≤ 14 selected variants). Ties (within 0.01) go to the parameter closer to MAGUS's setting / the cheaper one.

## Held-out evaluation (once)

Run the selected variants plus `raw:mcl:4` and `es4:mcl:4` on the 12 held-out replicates. Reported with
proteins and DNA/RNA separately: mean Δ, W/T/L (tie band 0.05), two-sided Wilcoxon signed-rank, ΔSPFN / ΔSPFP,
merge wall time, cluster statistics.

Primary questions:
1. **Alone:** does the training-best variant on the raw graph beat `raw:mcl:4` on held-out?
2. **Combined:** does the training-best variant on the es4 graph beat `es4:mcl:4` (and MAGUS)?
3. **Mechanism:** is the es4 gain (es4:mcl:4 − raw:mcl:4) smaller for a clustering method other than MCL
   (i.e. is the low-support damage specific to MCL)? Answered from training + held-out without further tuning.

Exploratory only (labelled as such): per data type best variants, Infomap, FastTree nRF on AliSim.
