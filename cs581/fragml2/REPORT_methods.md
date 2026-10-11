## 2. Prior art: what is known and what is not

Full notes, with citations, are in `prior_art.md`.

**Known.**
- *Backbone tree plus independent placement of fragments, without polishing.*
  - Smirnov & Warnow 2021 (*Syst Biol* 70:268, doi:10.1093/sysbio/syaa058) tested this. On 1000M1-HF, UPP+pplacer gives FN 48.8% vs UPP+RAxML-NG 37.0%, on estimated alignments. Their conclusion: placement is clearly worse than "align, then ML".
  - SEPP/UPP, SCAMPP and BSCAMPP place fragments but measure placement error, not tree error after polishing.
- *Grafting placements into a tree* is standard (`guppy tog`, `gappa examine graft`).
- *Placement or insertion followed by tree search* is known in other settings:
  - UShER + matOptimize (parsimony SPR; near-complete SARS-CoV-2 genomes; Thornlow et al. 2022);
  - PUmPER (RAxML-Light, 2014);
  - uDance (Balaban et al. 2024, *Nat Biotechnol* 42:768): APPLES-2 placement, then local unconstrained RAxML per partition and ASTRAL stitching. uDance **filters** short sequences instead of keeping fragments.
  - INC/Constrained-INC and GTM, the Warnow-lab divide-and-conquer methods. Park, Zaharias & Warnow 2021 benchmark these on 1000M1-HF.
- *Incomplete constraint trees* in RAxML-NG (`--tree-constraint`) and IQ-TREE (`-g`) are standard features. Constrained search around a reference tree goes back to EPA-era practice (Berger et al. 2011).

**Has "backbone tree + constrained joint ML of the fragments" been benchmarked on these fragmentary datasets, or compared with uDance or SEPP-style placement?**

As far as we can find, **no**:
- Smirnov & Warnow compare placement with align-then-ML, not constrained ML or graft-then-polish.
- Park et al. compare GTM, TreeMerge, Constrained-INC, IQ-TREE 2 and RAxML-NG, but not a backbone-constrained search.
- uDance's paper has no fragmentary-data benchmark of this kind.
- Two targeted web searches (2026-10-10) found no such comparison.

The only measurements we know of are this course's own earlier pilot (`cs581/ml`, branch `claude/cs581-ml`) and this one. **New here:**
1. A paired head-to-head of graft-only, graft+polish, constrained RAxML-NG, uDance and full RAxML-NG, with all steps' CPU time charged.
2. The effect of the EPA-ng rate-scaler/premasking bug inside a tree-building pipeline.

**Not new:**
- every building block;
- the idea of "place, then optimise";
- the finding that placement alone is worse than ML.

## 3. Data

All arms use the same true alignment of each replicate. The true tree is the reference.

| condition | replicates | taxa | fragments | source |
|---|---|---|---|---|
| **M1HF** (1000M1-HF) | R0–R4 | 1000 | 500 (mean ≈ 250 nt; median full length ≈ 1000 nt) | exact inputs of Park, Zaharias & Warnow 2021, Illinois Data Bank doi:10.13012/B2IDB-7008049_V1 (`1000M1_HF_Analysis.tar.gz`), rebuilt from their full-width subset alignments; ROSE true tree |
| M1HF | R5–R9 | 1000 | 500 | ROSE 1000M1 R5–R9 from the MAGUS paper's Datasets.zip (doi:10.13012/B2IDB-2643961_V1), fragmented by `code/make_frag.py` with the same protocol (50% of sequences → one contiguous fragment, length ~ N(0.25 × median, 60)); fragment lengths match R0–R4 (mean 250–256 nt) |
| RNASimHF | R0–R2 | 1000 | 500 | RNASim 1K (Datasets.zip), same protocol |
| RNASim10KHF | R0–R1 | 10,000 | 5,000 | RNASim 10K (Datasets.zip), same protocol. Columns with ≥ 95% gaps among the full-length sequences are masked (8,700 → 1,626 columns), as in the epang pilot, so that 5,000-tip EPA-ng fits in memory; all arms use the masked alignment |

The Smirnov & Warnow Dryad deposit has a dead DOI in the paper; the real one is doi:10.5061/dryad.95x69p8h8. It returned HTTP 403 from this VM, so it is not used.

## 4. Methods

`code/pipe.py` runs every arm on one replicate. Every tool uses 1 thread. CPU time is the user+sys time of all of an arm's steps, including the backbone tree, model evaluation, placement and grafting. Steps shared between arms are cached, but their cost is charged to every arm that uses them. Peak RSS is the maximum over the arm's steps.

**Baselines** (all sequences):
- FastTree 2.1.11 (`-gtr -gamma`);
- IQ-TREE 3.1.4 `--fast` (GTR+G);
- RAxML-NG 2.0.3 (GTR+G, `--tree pars{1}`, seed 1), the pre-registered reference;
- published Park et al. trees for M1HF R0–R4: RAxML-NG with 20 starts and a 24 h cap, IQ-TREE 2, and GTM with an IQ-TREE start;
- uDance v1.6.5 (commit d074760). Single-gene DNA run with the same backbone sequences and FastTree backbone tree, defaults otherwise, plus 3 local patches needed for a single-gene run (`udance/NOTES.md`).

**Backbone.** Sequences with ungapped length ≥ τ × median, where the median is over all sequences:
- τ = 0.5 keeps the 500 full-length sequences plus about 20 long fragments;
- τ = 0.75 keeps exactly the 500 full-length sequences.

**(A) Constrained ML, `constr_<BB>_<τ>`.**
1. Estimate the backbone tree with FastTree (`ft`) or IQ-TREE `--fast` (`iqf`).
2. Run RAxML-NG (1 parsimony start) on all sequences, with the backbone tree as a non-comprehensive `--tree-constraint`.

**(B) Place, then polish, `place_<BB>_<τ>_fix_<polish>`.**
1. Estimate the backbone tree.
2. Re-estimate its branch lengths and GTR+G parameters with RAxML-NG `--evaluate`.
3. Place each non-backbone sequence with patched EPA-ng 0.3.8 (`cs581/epang/code/epa-ng-fix.patch`, defaults otherwise).
4. Graft with `gappa examine graft --fully-resolve` (best-LWR placement).
5. Polish, one of:
   - `graft`: none. This is the Smirnov & Warnow-style row;
   - `rxfast`: RAxML-NG from the grafted tree in fast mode (`--opt-topology simplified --stop-rule kh-mult`). This is the pre-registered primary;
   - `ft`: FastTree `-intree`;
   - `iqfast`: IQ-TREE `-t … --fast`;
   - `rxfull`: a default RAxML-NG search.

**(C)** The same as (B), with stock EPA-ng 0.3.8 built from the same source tree.

**Metrics.**
- FN and FP vs the true tree (dendropy), as rates over internal edges.
- CPU, wall time and peak RSS.
- lnL of each final tree under one fixed model: that replicate's RAxML-NG GTR+G parameters, RAxML-NG `--evaluate --opt-model off`.

**Statistics.** Pre-registered (`PREREG.md`, pushed before any result except a single R0 smoke test of the graft arm):
- paired per replicate vs RAxML-NG;
- W/T/L with a 0.5-point tie band;
- two-sided Wilcoxon signed-rank test;
- Holm correction across the pipeline arms;
- CPU ratio (median and range).
