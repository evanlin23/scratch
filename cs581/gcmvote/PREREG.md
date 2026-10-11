# Pre-registration: a reference-free vote model for GCM edge support

AI-assisted (Claude), exploration code for CS581 project.

Written and pushed before any variant alignment was scored. Seen before writing: the control merge on BBA0101
(`magus` = 29.19 % error, reproduces the cached MAGUS run) and the fitted mixture parameters on BBA0101 (no
accuracy, no reference). Nothing else.

## Model (code/vote.py)

For every cross-subset GCM edge (a, b): k = number of backbones contributing ≥ 1 residue pair, n = number of
backbones holding a residue of a's subset in column a **and** a residue of b's subset in column b. Only edges
with k ≥ 1 exist, so both components are zero-truncated. Two-component mixture fitted by EM per run, no
reference: binomial (k ~ Bin(n, p1) true / Bin(n, p0) false, weight π) and beta-binomial (each component
BetaBin(n, mean, concentration)). Posterior w = P(true | k, n). Component 1 = the one with the larger mean.

## Variants (all merge-only on MAGUS's own subsets; paired against `magus`, MAGUS's merge through the same code)

| name | rule |
|---|---|
| `magus` | MAGUS's graph (control) |
| `es4` | delete edges with k < 4 (gcmgen's hand-picked filter) |
| `hard` / `hard-bb` | keep edges with w > 0.5 (binomial / beta-binomial) |
| `soft` / `soft-bb` | edge weight × w |
| `soft2`, `soft4` | edge weight × w^γ, binomial (γ grid {2, 4}; γ chosen on training only) |
| `fracF` | delete edges with k < ceil(F·B), F ∈ {0.2, 0.3, 0.4, 0.5} (sensitivity only) |
| `hard+mask` | tree variant: `hard`, then remove final columns whose cross-subset evidence (unfiltered edges whose two nodes ended in that column) has weight-weighted mean posterior < 0.5 |

## Data and split

gcmgen's fixed split (cs581/gcmgen/SPLIT.md), one MAGUS draw each (25 subsets, 10 L-INS-i backbones × 200,
MCL + minclusters):
- **training:** BBA0101, BBA0134, BBA0067, BBA0039, SIMMOD_R1, SIMHIGH_R1, 1000M2, 1000L1, 1000L2, 16S.M
- **held out:** BBA0154, BBA0190, SIMMOD_R2, SIMHIGH_R2, 1000L3, 1000M3, 1000S1, 1000S2, RNASim, 1000M4,
  1000S3, 1000M2_R1, 1000L1_R1, RNASim_R1 (BBA0081 / BBA0117 dropped as in gcmgen)
The AliSim sets are regenerated with gcmgen/gcmtrees's seeds and get fresh MAGUS draws here (cached DNA/RNA
and BAliBASE MAGUS inputs are reused). If a draw does not finish in the time budget it is reported as not run.

## Selection (training only)

The "best vote variant" is chosen among {`hard`, `soft`, `soft2`, `soft4`, `hard-bb`, `soft-bb`} as the
lowest mean Δ error (variant − magus, error = (SPFN+SPFP)/2 × 100) over the 10 training replicates (proteins and
DNA/RNA pooled); ties within 0.05 → the smaller worst-case Δ. Held-out replicates are scored with all
pre-declared variants but are never used to choose.

## Endpoints

- **Primary:** selected vote variant vs `magus` on the held-out replicates: mean Δ error, W/T/L (tie band
  0.05), two-sided Wilcoxon signed-rank (scipy default). Also vs `es4` on the same replicates.
  Reported for all held-out and separately for proteins and DNA/RNA.
- **Calibration:** reliability table per replicate: posterior bins vs measured precision (fraction of an edge's
  backbone residue-pair units that are true homologies, the gcmgen definition; and edge-level: share of edges
  whose units are mostly true). Effective cutoff = smallest k with w > 0.5, per exposure n.
- **Sensitivity:** B ∈ {5, 10, 20} (5 = MAGUS's backbones 1–5; 20 = MAGUS's 10 + 10 new L-INS-i backbones on new
  random 200-sequence sets, 8 per subset), backbone size 200; `frac0.2/0.3/0.4/0.5` vs `hard` (and the selected
  variant), on the training replicates (proteins first), as time allows.
- **Trees:** FastTree 2 `-lg -gamma` nRF vs the true tree on SIMHIGH R1–R4 (helpers: R5–R12) for the true
  alignment, `magus`, `es4`, the selected vote variant, and `hard+mask` (masked). Endpoint: Δ nRF vs `magus`,
  paired, Wilcoxon, pooled with the helpers' rows when available. "Improves trees" only if p < 0.05 with a
  negative mean.
- Verdict rule: "the vote model is a principled replacement for es4" if the selected variant is not worse than
  es4 on held-out (mean Δ vs es4 ≤ +0.1 and no data type worse by > 0.25) while needing no hand-picked k.
