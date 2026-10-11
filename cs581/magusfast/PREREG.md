# Pre-registration: MAGUS fast settings (written 2026-10-11, before any result of this pilot)

## Datasets

| role | dataset | type | source |
|---|---|---|---|
| train | BBA0101_R0 | protein (BAliBASE RV100, 322 seqs) | `cs581/data/balibase_clean/RV100_BBA0101.fasta` |
| train | BBA0190_R0 | protein (BAliBASE RV100, 343 seqs) | `cs581/data/balibase_clean/RV100_BBA0190.fasta` |
| train | SIMHIGH_R1 | protein (AliSim LG+G4, 1000 seqs) | regenerated as in branch claude/cs581-protbench |
| train | 1000M2_R0 | DNA (ROSE) | MAGUS paper data |
| train | 1000L1_R0 | DNA (ROSE) | MAGUS paper data |
| train | RNASim1000_R0 | RNA (RNASim) | MAGUS paper data |
| **held out** | BBA0067_R0 | protein (BAliBASE RV100, 274 seqs) | `cs581/data/balibase_clean/RV100_BBA0067.fasta` |
| **held out** | SIMHIGH_R2 | protein (AliSim, 1000 seqs) | regenerated as in claude/cs581-protbench |
| **held out** | 1000M3_R0 | DNA (ROSE, a condition not in training) | MAGUS paper data |
| **held out** | RNASim1000_R1 | RNA (RNASim) | MAGUS paper data |

Held-out data are not looked at until the setting is chosen (their default MAGUS runs may be queued
earlier for scheduling reasons, but no number from them is read before the choice is committed).

## Settings and the support-threshold rule

Grid: number of backbones N in {10, 6, 4, 3} x backbone size s in {200, 100} x support pruning {off, on}
(16 settings; N=10, s=200, off = MAGUS default). MAGUS flags otherwise the paper's (`-r N -m s`, 25 subsets,
FastTree guide tree on a 300-sequence skeleton, MCL, minclusters trace, no optimization).

**Support pruning rule (fixed here):** with N backbones, keep a cross-node graph edge only if it is
supported (the two residues share a column) by at least K(N) = ceil(0.4 N) backbones:
K(10) = 4 (the previously tested `#es4`), K(6) = 3, K(4) = 2, K(3) = 2.

Exact engineering speed-ups (vectorized graph build, MCL with 4 threads) are verified to give the
identical graph / clustering and are allowed in every candidate setting.

## Accuracy measure and selection rule (training data only)

- SP error = (SPFN + SPFP) / 2 x 100 (FastSP), on the full reference alignment.
- Accuracy per setting is measured **paired**: merge-only reruns on the default MAGUS run's own 25
  subalignments; N size-200 backbones = the first N of MAGUS's own 10; size-100 backbones = for each of
  those backbone sequence sets, a random half per subset (4 of 8, seed fixed) realigned with MAGUS's
  L-INS-i command. Delta = setting error - error of the default merge (10 x 200, no pruning) on the same subsets.
- Time per setting = measured guide-tree stage + subset/backbone pool time scaled by measured MAFFT CPU
  (subsets + that setting's backbones) + measured merge-only wall time of that setting.
- **Choice:** among the 16 settings, those with Delta <= +0.25 points on **every** training dataset are
  eligible; the chosen setting is the eligible one with the smallest geometric-mean predicted wall time
  ratio (setting / default) over the training datasets. Ties within 2% go to more backbones.
  If only the default-sized settings are eligible, the verdict on the 2x goal is "not promising" and the
  fastest setting within +0.5 on every training dataset is reported as a secondary, explicitly
  post hoc, option.

## Held-out evaluation

For each held-out dataset: MAGUS default end to end, and the chosen setting end to end (real wall
times, 4 cores, idle machine), plus the paired merge-only Delta on the default run's subsets.
"Equal accuracy" = mean paired Delta <= +0.25 and no held-out dataset worse than +0.5 points.
Speed-up = default wall / chosen-setting wall (geometric mean). Tree error: FastTree
(`-nt -gtr -gamma` for DNA, `-lg -gamma` for protein) nRF vs the true tree on ROSE and AliSim sets.
W/T/L use a 0.05-point tie band.
