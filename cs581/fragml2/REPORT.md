# Backbone → place → polish on fragmentary data (CS581 pilot `fragml2`)

Branch `claude/cs581-fragml2`, directory `cs581/fragml2/`. The pre-registration (`PREREG.md`) was pushed before any results were seen. The run was **stopped early** at the orchestrator's request, so some pre-registered arms are incomplete; §6 lists them.

## 1. Summary

**Primary (pre-registered): 1000M1-HF, true alignment, vs RAxML-NG (1 parsimony start, 1 thread).** The pipeline is a FastTree tree on the full-length sequences, then patched EPA-ng placement of the fragments, a gappa graft, and a RAxML-NG fast-mode polish.

| | primary pipeline | RAxML-NG |
|---|---|---|
| mean FN | 23.7% (n=10) | 24.4% (n=8) |
| CPU | 12.8 min | 35.4 min |

- Over the 8 paired replicates the FN difference is **+0.11 points** (W/T/L 3/1/4, Wilcoxon p = 0.84).
- The CPU ratio is **0.34** (median; range 0.30–0.45).
- This meets the pre-registered bar for "RAxML-NG-level at a fraction of the cost" (ΔFN ≤ +0.5, ratio ≤ 0.5). It is **not better** than RAxML-NG.

**Secondary comparisons, 1000M1-HF.**

| method | FN | vs RAxML-NG | CPU |
|---|---|---|---|
| Constrained RAxML-NG (earlier pilot's method) | — | ΔFN +0.40 (3/0/5), p = 0.31 | 0.36× |
| Placement only, graft, no polish (Smirnov–Warnow style) | 33.1% | ΔFN +9.1 (0/0/8), Holm p = 0.03 | 1 CPU-min |
| uDance | about the same as graft-only (§5.2) | 5–6 points worse than RAxML-NG, the primary pipeline and constrained RAxML-NG on the same leaves | 22 CPU-min |
| IQ-TREE 3 `--fast` | 37.4% | | |
| FastTree | 47.2% | | |
| Published GTM | 28.4% | | |
| Published IQ-TREE 2 | 30.2% | | |
| Published 24 h RAxML-NG | 24.9% | | |

- The graft-only row shows that the polish is what recovers about 9 FN points.
- uDance keeps only 87–96% of the taxa, and its run failed on 1 of 5 replicates.
- The published rows are on R0–R4; the primary pipeline scores 25.0% on those same replicates.

**Does the EPA-ng fix matter?**
- **Not at 1,000 taxa.** The backbone has about 520 tips, below the 2,000-tip trigger. Placements were identical for 100% of fragments on 9/9 replicates.
- **Yes at 10,000 taxa** (RNASim10K-HF, 5,040-tip backbone, n = 2):
  - Stock EPA-ng puts 41–44% of fragments on a different edge, and is falsely confident: LWR > 0.99 for 94% of fragments vs 41–44% patched.
  - Graft FN is **47.7% stock vs 35.6% patched**, i.e. +12 points.
  - After a cheap FastTree polish the gap shrinks to **32.2% vs 30.9%** (+1.3).
  - For reference, FastTree on all 10K sequences gets 64.3% FN. The patched pipeline gets 30.9% at 1.6× FastTree's CPU.

**Verdict: unclear, leaning not promising, as a 4-week accuracy project.** Details in §7.
- The pipeline robustly matches RAxML-NG at about 1/3 of the CPU, but does not beat it.
- That is the same outcome the earlier constrained-ML pilot already had, so the speed-up is not new.
- **The most interesting result is at 10K taxa.** There, place-then-polish beats FastTree by about 33 points, and the EPA-ng bug costs 12 points before polish. That is the direction we would pursue, but it is a different, larger-scale project.

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
| RNASimHF | R0 (only R0 ran) | 1000 | 500 | RNASim 1K (Datasets.zip), same protocol |
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
- lnL under one fixed model was pre-registered (`code/lnl.py`) but **not computed** (stopped early).

**Statistics.** Pre-registered (`PREREG.md`, pushed before any result except a single R0 smoke test of the graft arm):
- paired per replicate vs RAxML-NG;
- W/T/L with a 0.5-point tie band;
- two-sided Wilcoxon signed-rank test;
- Holm correction across the pipeline arms;
- CPU ratio (median and range).

## 5. Results

All tables are generated by `code/summarize.py` → `results/tables.md`, `code/udance_score.py` → `results/udance.md` and `code/cmp_place.py` → `results/place_cmp.md`. Raw rows are in `results/runs.jsonl`.

### 5.1 1000M1-HF (primary; true alignment)

| method | n | mean FN | ΔFN vs RAxML-NG (W/T/L) | Wilcoxon p | Holm p | mean CPU min | CPU ratio, median [range] | peak RSS MB |
|---|---|---|---|---|---|---|---|---|
| FastTree 2 | 10 | 47.2% | +23.2 (0/0/8) | 0.008 | – | 1.7 | 0.05 [0.05–0.06] | 116 |
| IQ-TREE 3 `--fast` | 9 | 37.4% | +13.0 (0/0/8) | 0.008 | – | 7.0 | 0.19 [0.16–0.25] | 463 |
| **RAxML-NG (1 parsimony start)** | 8 | 24.4% | – | – | – | 35.4 | – | 138 |
| published RAxML-NG (20 starts, 24 h cap; R0–R4) | 5 | 24.9% | +0.16 (1/1/3) | 0.88 | – | – | – | – |
| published IQ-TREE 2 (R0–R4) | 5 | 30.2% | +5.4 (0/0/5) | 0.06 | – | – | – | – |
| published GTM, IQ-TREE start (R0–R4) | 5 | 28.4% | +3.6 (0/0/5) | 0.06 | – | – | – | – |
| **(B) primary: FastTree backbone τ=0.5 → patched EPA-ng → graft → RAxML-NG fast polish** | 10 | 23.7% | **+0.11 (3/1/4)** | **0.84** | 0.84 | **12.8** | **0.34 [0.30–0.45]** | 618 |
| (A) constrained RAxML-NG, FastTree backbone τ=0.5 | 9 | 24.5% | +0.40 (3/0/5) | 0.31 | 0.63 | 12.3 | 0.36 [0.33–0.37] | 128 |
| (B) graft only, no polish (Smirnov–Warnow style), patched EPA-ng | 9 | 33.1% | +9.1 (0/0/8) | 0.008 | 0.031 | 1.0 | 0.03 | 618 |
| (C) graft only, stock EPA-ng | 9 | 33.1% | +9.1 (0/0/8) | 0.008 | 0.031 | 1.0 | 0.03 | 618 |

How to read this table:
- n differs between rows because RAxML-NG on R8–R9 was stopped at wrap-up.
- Paired statistics use the replicates both arms finished. For the primary that is R0–R7, n = 8; p = 0.84 is the exact two-sided Wilcoxon on 8 pairs.
- The CPU time includes every step of an arm: backbone tree, RAxML-NG `--evaluate`, EPA-ng, graft and polish.

Per-replicate FN (%):

| method | R0 | R1 | R2 | R3 | R4 | R5 | R6 | R7 | R8 | R9 |
|---|---|---|---|---|---|---|---|---|---|---|
| RAxML-NG | 21.8 | 23.3 | 25.7 | 23.2 | 29.8 | 23.6 | 23.8 | 23.9 | – | – |
| primary (B, rxfast) | 22.7 | 22.9 | 27.1 | 24.2 | 28.2 | 22.9 | 25.7 | 22.3 | 21.3 | 20.2 |
| constrained (A) | 20.9 | 23.9 | 25.2 | 25.1 | 31.2 | 22.6 | 24.5 | 24.9 | 22.7 | – |
| graft only | 30.6 | 31.5 | 35.6 | 36.1 | 37.7 | 28.5 | 33.9 | 33.9 | 30.1 | – |
| FastTree | 48.9 | 46.5 | 48.9 | 49.4 | 50.8 | 43.8 | 46.2 | 46.2 | 47.5 | 43.6 |
| backbone FastTree tree, 520 taxa | 12.4 | 12.5 | 12.9 | 11.7 | 16.9 | 12.1 | 14.2 | 15.5 | 15.3 | 13.0 |

Observations:
- **Robustness.** The primary is within ±2 points of RAxML-NG on every replicate. Its worst cases are R2 (+1.4) and R6 (+1.9).
- On R4, the replicate with the worst backbone (16.9% FN), the constrained arm suffers (31.2%) because the constraint freezes the backbone errors. The unconstrained polish recovers (28.2%, better than RAxML-NG's 29.8%). This is the one qualitative advantage of (B) over (A).
- The FP rate tracks FN within 0.4 points, since all trees are nearly binary.
- The R0–R4 RAxML-NG numbers equal the earlier pilot's (21.8/23.3/25.7/23.2/29.8), so the setup reproduces.

### 5.2 uDance (secondary; M1HF R0–R4)

uDance v1.6.5 ran single-threaded on the same backbone and backbone tree (`udance/NOTES.md`):
- It needed 3 local patches for a single-gene run.
- It finished R0, R1, R2 and R4. **R3 failed** in its final `stitch` rule; not debugged.
- It drops taxa: 874–959 of 1000 kept (unplaced queries, TreeShrink, and losses at stitching).

Every tree below is therefore pruned to uDance's leaf set before scoring (n = 4):

| | uDance | graft only | primary (B) | constrained (A) | RAxML-NG | IQ-TREE `--fast` | FastTree |
|---|---|---|---|---|---|---|---|
| mean FN on uDance's leaves | 27.2% | 29.2% | 22.1% | 21.9% | 21.7% | 33.8% | 44.1% |
| CPU min | 22.6 | 1.0 | 12.8 | 12.3 | 35.4 | 7.0 | 1.7 |

uDance is a little better than graft-only but 5 points worse than the ML pipelines, at about twice their CPU.
- This is not a fair test of uDance's design goal, which is multi-gene species trees with updatable backbones.
- On single-gene fragmentary data, its per-partition local re-estimation does not recover what a global ML polish does.

### 5.3 RNASim1K-HF (secondary; incomplete)

Only R0 finished, and only for three arms, before the stop:
- FastTree 58.9%;
- graft only 42.3%;
- **primary 34.7%** (16.1 CPU-min).

There is no paired RAxML-NG run from this pilot. The earlier pilot (`cs581/ml`, `results/test.jsonl`) generated the identical replicate: the same `make_frag.py` seed, and its FastTree FN is also 58.9%. Its RAxML-NG 2.0.3 FN on R0 is 34.0%, about 42 CPU-min. That suggests the same picture (+0.7 points at about 0.4× CPU), but n = 1 and the runs come from different pilots.

### 5.4 Does the EPA-ng fix matter? (B vs C)

| data | backbone tips | fragments with the same best edge, patched vs stock | LWR > 0.99, patched / stock | graft FN, patched / stock | + FastTree polish, patched / stock | FastTree on all |
|---|---|---|---|---|---|---|
| M1HF, R0–R8 | 517–526 | **100%** (all 9 reps) | 60–68% / identical | 33.1% / 33.1% | – | 47.2% |
| RNASim10K-HF R0 | 5,040 | 55.6% | 43.9% / 94.2% | 35.9% / 48.8% | 30.6% / 32.2% | 63.3% |
| RNASim10K-HF R1 | ~5,040 | 58.9% | 41.1% / 94.0% | 35.3% / 46.6% | 31.2% / 32.2% | 65.4% |

- The stock binary's best-placement log-likelihoods are inflated by about 1.5 × 10⁵ log units per fragment, which matches the epang pilot's diagnosis.
- Inside the pipeline, the bug costs **12 FN points** at the graft stage and **1.0–1.6 points** after a FastTree polish (2/2 replicates).
- In other words, polishing repairs most of the damage but not all.
- Cost at 10K, per replicate: about 10 CPU-min for graft (EPA-ng peak RSS 7.9 GB on the masked alignment), about 19 CPU-min with FastTree polish.
- A RAxML-NG polish or RAxML-NG reference at 10K was not run (hours per run), so we cannot say how this compares with full ML at that scale.

## 6. What was not run (stopped early)

Pre-registered but not finished:
- lnL re-evaluation under a fixed model (`code/lnl.py`; written, not run);
- the rest of the (A) grid: IQ-TREE `--fast` backbone, τ = 0.75;
- the other polishes at 1K: FastTree, IQ-TREE `--fast`, full RAxML-NG;
- the IQ-TREE `--fast` backbone for (B);
- IQ-TREE 3 default;
- RNASimHF R1–R4 and its RAxML-NG baseline;
- RAxML-NG on M1HF R8–R9.

The estimated-alignment setting (UPP) was not attempted. The earlier pilot found UPP with default settings gives about 91% FN trees on these data unless its backbone length filter is fixed.

## 7. Verdict for a 4-week CS581 project: **unclear (leaning not promising) as stated; promising if reframed to large trees**

**Why not promising as stated.**
- On 1000-taxon fragmentary data, the pre-registered pipeline is reliably RAxML-NG-level at about 1/3 of the CPU, but no better.
- The earlier pilot's constrained-ML heuristic already gave the same speed/accuracy point, so a 4-week project would add robustness evidence but no new capability.
- Place-then-optimise is known practice (UShER/matOptimize, PUmPER, uDance), which limits novelty to "benchmarked on the Warnow-lab fragmentary data". A plausible but modest course result.

**What is promising.** At 10K taxa:
- Placement plus a FastTree polish beats FastTree on everything by 33 points (31% vs 64% FN) for 1.6× its CPU.
- The EPA-ng bug measurably hurts such pipelines: 12 points at graft, about 1 point after polish.

A large-tree project, in which RAxML-NG is too slow to be the default, is where "backbone → place → polish" could actually change practice. It would cover:
- RNASim 10K/50K HF;
- polish with FastTree vs a RAxML-NG fast mode;
- references from the published Park et al. RNASim10K/50K RAxML-NG and GTM trees (IDB-7008049).

### Weeks 1–4 (if pursued, large-tree framing)

| week | work | done when |
|---|---|---|
| 1 | RNASim10K-HF, 5 reps (and 50K-HF, 2 reps), via `prep.py` + `make_frag.py`. Baselines: FastTree, IQ-TREE `--fast`, published Park et al. RAxML-NG and GTM trees on the non-fragmentary versions as context. Finish the 1K lnL and the remaining 1K arms. | baseline table |
| 2 | Pipelines at 10K: backbone FastTree vs IQ-TREE `--fast`; patched vs stock EPA-ng (and BSCAMPP for 50K); polish = FastTree, RAxML-NG fast mode with a time cap, or IQ-TREE `--fast`. | ≤ 3 frozen variants |
| 3 | Held-out reps; 16S.M HF (real data); one estimated-alignment arm (UPP with a fixed backbone length or WITCH). | all runs |
| 4 | Accuracy vs CPU curves; the size of the bug's effect vs tree size; write-up. | report and slides |

### Risks
- **No accuracy gain over RAxML-NG at 1K.** This is already observed. At 10K the reference is weaker (FastTree), so the gains are real but come against a weaker baseline.
- **Memory.** 5,000-tip EPA-ng needs about 8 GB even on a masked alignment. 50K needs BSCAMPP or `--memsave`.
- **Estimated alignments** can dominate the error (UPP backbone selection). Budget time to configure UPP/WITCH correctly.
- **RAxML-NG references at 10K+** take many CPU-hours. Rely on the published Park et al. trees, or on a 4-thread budget on the campus cluster.
- **uDance comparisons** need local patches for single-gene data and drop taxa. Any comparison must score on the common leaf set, as done here.
