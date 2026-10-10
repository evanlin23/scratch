# Pre-registration: backbone → place → polish on fragmentary data (CS581 pilot `fragml2`)

Written and pushed **before any tree in this pilot was estimated or scored** (2026-10-10).
Code: `code/pipe.py` (arms), `code/prep.py` (data), `code/summarize.py` (tables, written later; it
implements only what is stated here, and any deviation is listed under "Amendments" below).

## Question
Can a backbone-then-place-then-polish pipeline give RAxML-NG-level (or better) trees on fragmentary
data at a fraction of the cost, robustly? Does the EPA-ng premasking/rate-scaler fix
(`cs581/epang/code/epa-ng-fix.patch`) matter inside such a pipeline?

## Data (true alignments for every method)
- **M1HF** (primary): ROSE 1000M1 with 500/1000 sequences cut to one fragment each (~25% of median).
  R0–R4 = the exact Park, Zaharias & Warnow 2021 inputs (IDB-7008049); R5–R9 = ROSE 1000M1 R5–R9
  (MAGUS Datasets.zip, IDB-2643961) fragmented with `make_frag.py` (same protocol). n = 10.
- **RNASimHF** (secondary): RNASim 1K R0–R4, same fragmentation protocol. n = 5.
- **RNASim10KHF** (EPA-ng question; > 2,000 backbone tips): RNASim 10K R0–R1, same protocol
  (5,000 fragments, ~5,040-tip backbone). Analysed on the true alignment restricted to columns
  with < 95% gaps among the full-length sequences (memory; as in the epang pilot), for all arms.
- The Smirnov & Warnow 2021 Dryad deposit (doi:10.5061/dryad.95x69p8h8) returned 403 from this VM,
  so it is not used.

## Arms
Baselines: `base_fasttree`, `base_iqfast` (IQ-TREE 3 `--fast`), `base_raxmlng` (RAxML-NG 2.0.3, GTR+G,
one parsimony start, 1 thread), `base_iqtree` (IQ-TREE 3 default; lowest priority, run if time
allows); published Park et al. trees (RAxML-NG 24 h, IQ-TREE 2, GTM) for M1HF R0–R4.

Pipelines (backbone = sequences with ungapped length ≥ τ × median; BB ∈ {ft = FastTree, iqf = IQ-TREE `--fast`}):
- (A) `constr_<BB>_<τ>`: RAxML-NG (1 parsimony start) on all sequences, backbone tree as `--tree-constraint`.
  Grid: BB ∈ {ft, iqf} × τ ∈ {0.5, 0.75}.
- (B) `place_<BB>_<τ>_fix_<polish>`: RAxML-NG `--evaluate` on the backbone → patched EPA-ng places all
  other sequences → `gappa examine graft --fully-resolve` (best placement) → polish ∈
  {graft (none), rxfast (RAxML-NG fast mode from the grafted tree), iqfast (IQ-TREE `--fast -t`),
  ft (FastTree `-intree`), rxfull (default RAxML-NG search from the grafted tree)}.
- (C) the same with stock EPA-ng 0.3.8 (`_stock_`).

**Primary pipeline (declared now, single):** `place_ft_0.5_fix_rxfast`.

## Primary analysis
M1HF, n = 10: `place_ft_0.5_fix_rxfast` vs `base_raxmlng` on FN, paired per replicate, two-sided
Wilcoxon signed-rank (scipy, exact), W/T/L with a 0.5-point tie band; plus the per-replicate CPU-time
ratio (pipeline / RAxML-NG, all pipeline steps included), reported as median and range.
Interpretation fixed in advance: "RAxML-NG-level at a fraction of the cost" = mean ΔFN ≤ +0.5 points
and median CPU ratio ≤ 0.5. Better = mean ΔFN < −0.5 with p < 0.05.

## Secondary (descriptive; Holm correction across the pipeline arms on M1HF where p-values are quoted)
- All other arms vs `base_raxmlng`, same statistics, on M1HF and RNASimHF.
- FP rate, lnL of every final tree re-evaluated with RAxML-NG `--evaluate` under one fixed model
  (that replicate's `base_raxmlng` GTR+G parameters, `--opt-model off`, branch lengths optimised),
  CPU, wall, peak RSS.
- (B) vs (C): on 1000-taxon data the backbone has ≈ 500 tips (< 2,000), so the bug cannot trigger;
  prediction: identical placements. On RNASim10KHF: graft FN and graft + FastTree-polish FN, patched vs
  stock (n = 2, descriptive), vs `base_fasttree` on all 10K sequences.

## Priorities (compute: 4 cores, ~4 h)
1. M1HF R0–R9: base_fasttree, base_iqfast, base_raxmlng, primary pipeline, constr_ft_0.5, graft arms.
2. M1HF: rest of the A/B grid, (C) stock graft; RNASim10KHF R0.
3. RNASimHF R0–R4: base_fasttree, base_raxmlng, primary, constr_ft_0.5, place_ft_0.5_fix_ft.
4. base_iqtree, RNASim10KHF R1, anything left.
Arms not finished in time are reported as not run.

## Amendments
(none yet)
