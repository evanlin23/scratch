# Pilot: Is ASTRID-multi consistent under GDL? Does an "ASTRID-Pro" distance fix it?

*CS581 project pilot, 2026-10-09/10. Wall-clock about 4 h on 4 cores. Branch `claude/cs581-gdl`. Code is in `code/`, result files in `results/`.*

## 0. Verdict

**Promising, as a theory-led project. Not promising as a new method that beats ASTRAL-Pro or DISCO.**

1. **Theory (new, clean, verified numerically).**
   - Under pure GDL, with ideal rooting and tagging and the root counted, the expected *speciation-only ortholog* internode distance is always a tree metric.
   - It has the species-tree topology **when no branch is supercritical** (λ_e ≤ μ_e on every branch).
   - With one supercritical internal branch it converges to the **wrong** tree. Predicted vs observed quartet sums: 1.149 vs 1.1515 / 1.1508 (20,000 families).
   - A **survival-reweighted correction** (ASTRID-Pro-S) removes the bias. Survival probabilities estimated by reconciliation match theory (ŝ = 0.049 vs e⁻³ = 0.0498; 0.695 vs s_yx). The correction recovers the true tree in both failing configurations, given a correct first-pass species tree.
   - Iterating from its own wrong tree does not escape: the wrong tree is a fixed point.
2. **ASTRID-multi: no inconsistency found; the evidence leans consistent.**
   - The best candidate from the 4–5-taxon search was wrong at 20,000 families in 4 of 4 samples.
   - With **2 million** families (jackknife SE), all five quartet margins are positive: three clearly (z = 5–30) and two weakly (z ≈ 1.1).
   - On DISCO true gene trees (GDL + ILS), ASTRID-multi converges more slowly than the others: about 1 wrong branch at 1,000 genes. By 10,000 genes, 4 of 5 replicates are error-free.
   - So ASTRID-multi looks consistent but sample-inefficient under GDL. That is not a proof, and heavy-tailed copy numbers make Monte Carlo a poor tool for settling it.
3. **Accuracy on estimated gene trees.**
   - *FastMulRFS data* (252 held-out runs): ASTRID-Pro **ties** ASTRID-multi (p = 0.90) and is slightly behind ASTRID-DISCO. All ASTRID variants beat ASTRAL-Pro (p = 4e-7) and FastMulRFS.
   - *DISCO data, pre-registered held-out test* (33 runs, `results/PREREG_disco_heldout.md`): ASTRID-Pro **beats ASTRID-multi**, −0.0096 FN rate, 16/12/5, Holm-adjusted p = 0.021. The pre-registered criterion is met.
     - The gain comes from the high-duplication conditions. In gdl_1e-9_05 the FN rate drops from 0.094 to 0.051.
     - ASTRID-Pro **ties** ASTRID-DISCO (p = 0.18) and ASTRAL-Pro (p = 0.30).

**What a 4-week project should be.** "Distance-based species-tree estimation under GDL: when is it consistent?" Contents:
- the additivity theorem and the supercritical counterexample;
- the survival-reweighted correction;
- ASTRID-multi analysed with exact expectations;
- empirical confirmation that the orthology-restricted distance fixes ASTRID-multi's high-GDL degradation.

The method is a *better ASTRID-multi*, not a better ASTRAL-Pro. Section 9 has the plan.

## 1. Question

Source: lecture deck *Phylogenomics part 2*, slide 35, "(Some) Open Questions" (`../literature/slide_open_problems.md`, item 7):

- "Can we find a distance correction for GDL so that ASTRID and NJst are statistically consistent under GDL models?"
- "Unknown if distance-based species tree estimation (e.g., ASTRID-multi) is statistically consistent under GDL models."
- Related, on the same slide: "Which other existing species tree estimation methods are statistically consistent under GDL and DLCOAL models?" and "Is ASTRAL-Pro statistically consistent for GDL under some random model of error for rooting and tagging?"

## 2. Prior art

Full table with DOIs (checked against Crossref): `results/prior_art.md`.

- **Proved consistent under GDL:**
  - ASTRAL-one and ASTRAL-multi: Legried, Molloy, Warnow, Roch, *JCB* 2021, doi:10.1089/cmb.2020.0424.
  - Quartet methods under DLCoal: Markin & Eulenstein, *Bioinformatics* 2021, doi:10.1093/bioinformatics/btab414.
  - ASTRAL-Pro, given correct rooting and tagging: Zhang et al., *MBE* 2020, doi:10.1093/molbev/msaa139.
  - ASTRAL-DISCO, given correct rooting and tagging: Willson et al., *Syst Biol* 2022, doi:10.1093/sysbio/syab070.
  - FastMulRFS, when GDL is not adversarial: Molloy & Warnow, *Bioinformatics* 2020, doi:10.1093/bioinformatics/btaa444.
- **Distance methods.**
  - ASTRID-multi, NJst and STAG under GDL appear only in simulation studies: Legried et al. 2021; Willson, Roddur, Warnow, AlCoB 2021, doi:10.1007/978-3-030-74432-8_8; STAG, doi:10.1101/267914.
  - Closest analogue: NJst/ASTRID are *inconsistent* under random missing data (Rhodes, Nute, Warnow, arXiv 2001.07844). Gene loss also deletes taxa.
- **ASTRAL-Pro 2** (doi:10.1093/bioinformatics/btac620) is an algorithmic speed-up with no new theory.
- **Bottom line.** I found no published result, positive or negative, on consistency of ASTRID-multi or NJst under GDL. I also found no GDL-corrected internode distance with a proof. The search was not exhaustive: recent Warnow-lab and Mirarab-lab preprints should be checked once more before the proposal.

## 3. Data and baseline reproduction

| dataset | DOI | used |
|---|---|---|
| FastMulRFS data | 10.13012/B2IDB-5721322_V1 | 100 species. DL rate 1e-10, 2e-10, 5e-10 × Ne 1e7, 5e7, 10 replicates each. RAxML gene trees from 25 and 100 bp. First 25, 100 and 500 genes. Published species trees. |
| DISCO data | 10.13012/B2IDB-4050038_V1 | `trees.tar.gz` (5.7 GB, downloaded once). Conditions: `default`, `gdl_1e-9_1` and `gdl_1e-9_05` with estimated trees from 100 bp; `gtrees_10000_l1` with true gene trees (GDL + ILS). |
| GDL-comparison data | 10.13012/B2IDB-2418574_V1 | exp2 downloaded (18 MB). Not used further, because its conditions overlap with the DISCO data. |

**Baseline reproduction** (`code/rescore_published.py`, `results/rescore_fastmulrfs_published.csv`). I rescored the authors' published species trees against the true species trees with my own RF code and compared with their `data-error-and-timings-ntaxa-100.csv`:

| method | trees | exact match (\|Δ\| < 1e-5) |
|---|---|---|
| ASTRAL-multi | 1200 | **1200** |
| MulRF | 1200 | **1200** |
| FastMulRFS (single) | 1200 | **1200** |
| STAG | 370 | 279. The 91 mismatches are unexplained; the published STAG files are probably a different run from the CSV. |

The published ASTRID(-multi) trees are flagged as buggy in the README and were not used.

## 4. Methods implemented (`code/methods.py`, about 350 lines)

- **Root and tag.**
  - Rooting: the root minimises the number of duplication nodes, found by an O(n) rerooting DP. Ties are broken by a loss count at duplication nodes (ASTRAL-Pro/DISCO style).
  - Tagging: a node is a duplication iff its children's species sets overlap.
  - For simulated trees, the true tags and root can be used instead (`truetags`).
- **ASTRID-multi.** A re-implementation of ASTRID's `IndSpeciesMapping::average` (read from the ASTRID source):
  - per gene, average the internode distance over all copy pairs of species A and B;
  - then average over the genes that contain both species;
  - the degree-2 root is not counted.
- **ASTRID-Pro.** Same averaging, but only over orthologous pairs (LCA tagged speciation), counting only speciation nodes on the path. Variants:
  - `-min`: closest orthologous copy per gene, STAG-like;
  - `-w`: each counted node weighted by the support of its parent branch, on DISCO-data trees that carry support;
  - `count_root`: also count the gene-tree root (needed for the theory, Section 5);
  - ablations: orthologous pairs with all nodes counted; all pairs with only speciation nodes counted.
- **ASTRID-DISCO.** My re-implementation of DISCO's decomposition, because the DISCO GitHub repository was not reachable from this machine. At each duplication node (post-order) the largest child subtree is kept and the others are split off as single-copy trees; trees with at least 4 leaves are kept. ASTRID is then run on the result.
- **Tree building.** FastME 2 (`-m B -n B -s`: balanced minimum evolution with NNI and SPR). Missing matrix entries are filled with the maximum observed distance; this was never needed on the benchmark data.
- **ASTRAL-Pro.** The `astral-pro` binary from bioconda `aster` (v1.25.3.8), which identifies itself as *ASTRAL-Pro3*, the successor of ASTRAL-Pro 2. Default settings, 1 thread.
- **FastMulRFS, ASTRAL-multi, MulRF, STAG.** The authors' published trees on the same replicates (FastMulRFS data only).
- **Simulator** (`code/gdlsim.py`): pure GDL birth–death along a species tree with per-branch (λ, μ), no ILS, emitting true D/S tags. Also generates Yule species trees.

## 5. Theory (`results/theory.md`)

**Lemma.**
- Setting: pure GDL, true gene trees, true root and tags, and families not conditioned on size.
- For an orthologous pair (a, b) with species LCA L:
  E[d(a, b)] = 1 + Σ_{internal v on path_S(A,B), v ≠ L} s(off(v)).
  Here s(c) is the probability that a copy entering the branch above c survives.
- Why: the observed path keeps the speciation node at v iff the daughter copy that goes off the path survives. That event is independent of the existence of a and b.

**Four-point consequence.**
- For every quartet, the two larger sums are **equal**, so the limiting matrix is always a tree metric.
- Balanced quartets are always resolved correctly.
- An induced caterpillar (((A,B)x,C)y,D) is resolved correctly iff
  s_yx < 2P_xy + s_xA + s_xB + s_yC.
  - s_yx is the survival probability of a copy entering the branch from y toward x.
  - P_xy is the sum of off-path survivals strictly between x and y.
  - s_xA, s_xB and s_yC are the survivals of copies entering x's child branches and y's child branch toward C.

**Sufficient condition.** If λ_e ≤ μ_e on every branch, the expected copy number never grows. A union bound then gives the inequality, so ASTRID-Pro (ideal tags, root counted) followed by NJ or FastME is **statistically consistent**.

**Counterexample.**
- Species tree (((A:1,B:1)x:1,C:1)y:1,D:3). Branch x has (λ, μ) = (3, 0); branches A, B and C have (0, 3).
- Predicted limiting sums: AC|BD = AD|BC = 1 + 3e⁻³ = 1.149 vs AB|CD = 1 + s_yx.
- Simulated with 20,000 families (`results/counterexample*.jsonl`): AB|CD = 1.695 vs AC|BD = 1.1515 and AD|BC = 1.1508. With the root counted, all three sums shift by exactly +1.
- The tree is wrong at 2,000 and 20,000 families.
- On the same families these are correct: ASTRAL-Pro, ASTRID-multi, ASTRID-DISCO, and ASTRID-Pro with *inferred* tags (luckily: hidden paralogs get tagged as speciation).

**Implementation lesson found by the search.**
- With true tags, the root **must** be counted. ASTRID's convention of skipping the degree-2 root removes 1 from every pair that crosses the root split. That can make the root-split edge negative.
- Effect: 3 balanced 4-taxon configurations failed without root counting and were fixed by it (`results/recheck_root.jsonl`).
- After root counting, true-tag ASTRID-Pro failed in only 2 of the confirmed configurations. **Both contain a supercritical internal branch**, exactly as the theory predicts: x = (2, 0) in one and y = (0.3, 0) in the other.
- On estimated trees, counting the inferred root changes nothing: dev replicates, 108 runs, mean FN 9.04 vs 8.97, p = 0.23 (`results/dev_rootcount.jsonl`).

**Not covered:**
- conditioning on families with at least 4 species;
- estimated tags (hidden paralogy);
- the closest-copy variant;
- ILS;
- ASTRID-multi.

**A natural next step.** δ(A,B) − 1 is linear in the unknown survival probabilities s(·). These can be estimated from presence/absence and copy-number data, so a reweighted distance could restore the correct topology when a branch is supercritical. Untested.

## 6. Results

### 6.1 Random search for limiting failures (4–5 taxa, true gene trees; `code/search.py`, `confirm.py`)

- **Search.** 98 random configurations with per-branch λ ∈ {0, 0.3, 1, 2} and μ ∈ {0, 0.3, 1, 2, 4}, branch lengths 0.1–1.5, 2,000 families each. Number of configurations where each method had FN > 0:
  - ASTRID-multi 8;
  - true-tag ASTRID-Pro (root not counted) 10;
  - ASTRID-Pro with inferred tags 5;
  - closest-copy 6;
  - ASTRID-DISCO 1.
- **Confirmation.** The 14 flagged configurations were rerun at 5,000 and 20,000 families with ASTRAL-Pro (`results/confirm.jsonl`):
  - Many flags vanished: small-sample noise.
  - Two configurations had *every* method wrong, ASTRAL-Pro included. These are very lossy settings and are uninformative.
  - The true-tag ASTRID-Pro failures are explained by root exclusion or by supercritical branches (Section 5).
- **ASTRID-multi candidate.** Tree ((((A:1,B:2)x:0.3,C:2)y:1,D:1)z:0.1,E:2). Branch D has (4, 0); C, A and E have high loss. Results by sample size:

  | families | seeds | ASTRID-multi wrong | which quartets |
  |---|---|---|---|
  | 20,000 | 4 | 4 | ACDE → AD\|CE in all 4 |
  | 80,000 | 3 | 1 | different quartets |

  - ASTRAL-Pro and true-tag ASTRID-Pro are correct (`results/margin_candidate1*.jsonl`).
  - Copy numbers in D are geometric with mean about 55, so the per-gene averages are heavy-tailed and the matrix converges very slowly.
  - **High-precision follow-up** (`code/hiprec.py`, `results/hiprec_candidate1.md`): 2,000,000 families in 400 batches (2 seeds), limiting matrix with delete-one-batch jackknife SEs. No family exceeded the copy cap, so the model is not truncated.

    | quartet | margin (wrong − true split) | SE | z |
    |---|---|---|---|
    | ABCD | +0.355 | 0.012 | +30.5 |
    | ABCE | +0.498 | 0.093 | +5.4 |
    | ABDE | +0.616 | 0.050 | +12.3 |
    | ACDE | +0.119 | 0.100 | +1.2 |
    | BCDE | +0.106 | 0.098 | +1.1 |

  - **Conclusion.** The 20,000-family "failures" were sampling noise. All limiting margins are positive, so this configuration leans *consistent*, but two quartets are only at about 1 SE. Settling it, and ASTRID-multi in general, needs exact expectations (generating functions of the birth–death process), not sampling.

### 6.2 Error vs number of families on 30-taxon trees (true gene trees, pure GDL; `results/curve30.md`)

Mean FN out of 27 internal branches, 3–4 replicate species trees:

| setting (λ, μ) | families | ASTRID-multi | ASTRID-Pro (true tags) | ASTRID-Pro (inferred) | ASTRID-DISCO | ASTRAL-Pro |
|---|---|---|---|---|---|---|
| critical (2,2) | 25 | 1.67 | 1.00 | 0.33 | 0.67 | 0.33 |
| | 100 | 0.33 | 0 | 0 | 0 | 0.33 |
| | 500 | 0 | 0 | 0 | 0 | 0 |
| supercritical (2,1) | 25 | 0.50 | 0 | 0 | 0 | 0 |
| | 100 / 500 | 0 | 0 | 0 | 0 | 0 |
| adversarial: internal (3,0.5), terminal (0,3) | 25 | 3.50 | 1.75 | 1.00 | 1.00 | 2.00 |
| | 100 | 0.50 | 0 | 0 | 0 | 0 |
| | 500 | 0 | 0 | 0 | 0 | 0 |

- **No plateau** in generic settings. Failures need specific branch-level rate patterns.
- ASTRID-multi is consistently the slowest to converge.

### 6.3 FastMulRFS data, estimated gene trees (`results/fmrfs_comparison.md`, `results/fmrfs_runs.jsonl`)

**Protocol.**
- 6 conditions × 10 replicates × {25, 100} bp × {25, 100, 500} genes = 360 paired runs.
- Pre-declared split: replicates 01–03 are dev and pick the ASTRID-Pro variant; replicates 04–10 (252 runs) are held out.
- Dev picked the closest-copy variant (0.0916 vs 0.0925 for the mean variant, a negligible difference).
- Metric: species-tree FN rate. All trees are binary, so FN rate = RF rate.
- Tie means equal FN count.
- Test: two-sided Wilcoxon signed-rank, paired by condition, replicate, sequence length and number of genes.

**Held-out, mean FN rate:**

| bp | genes | ASTRID-multi | ASTRID-Pro (mean) | ASTRID-Pro (min) | ASTRID-DISCO | ASTRAL-Pro | FastMulRFS | ASTRAL-multi (pub) | MulRF (pub) |
|---|---|---|---|---|---|---|---|---|---|
| 25 | 25 | **0.229** | 0.244 | 0.252 | 0.254 | 0.284 | 0.326 | 0.357 | 0.359 |
| 25 | 100 | 0.105 | 0.106 | 0.109 | **0.103** | 0.130 | 0.153 | 0.227 | 0.163 |
| 25 | 500 | 0.050 | 0.047 | **0.046** | **0.046** | 0.058 | 0.073 | 0.137 | 0.070 |
| 100 | 25 | 0.101 | 0.094 | 0.092 | **0.089** | 0.100 | 0.121 | 0.143 | 0.116 |
| 100 | 100 | 0.045 | **0.041** | 0.043 | **0.041** | 0.047 | 0.057 | 0.076 | 0.051 |
| 100 | 500 | 0.022 | 0.019 | 0.020 | **0.017** | 0.019 | 0.024 | 0.036 | 0.019 |
| all | all | 0.0920 | 0.0919 | 0.0937 | **0.0916** | 0.1064 | 0.1256 | 0.1626 | 0.1298 |

**Held-out paired tests**, ASTRID-Pro (min, the pre-registered choice) vs each baseline. W = ASTRID-Pro better.

| baseline | mean diff | W/T/L | p |
|---|---|---|---|
| ASTRID-multi | +0.0018 | 90/72/90 | 0.90 |
| ASTRID-DISCO | +0.0021 | 67/85/100 | 0.015 |
| ASTRAL-Pro | −0.0127 | 131/45/76 | 3.6e-7 |
| FastMulRFS | −0.0319 | 174/33/45 | 1.9e-23 |
| ASTRAL-multi (published) | −0.0689 | 221/22/9 | 8.4e-38 |

**Strata.**
- With 25 bp gene trees (high gene-tree error), ASTRID-multi beats ASTRID-Pro: +0.0079, p = 0.024.
- With 100 bp, ASTRID-Pro beats ASTRID-multi: −0.0043, p = 0.004.

A plausible reading: orthology restriction throws away pairs, which hurts when gene trees are noisy, and the GDL correction helps when they are accurate.

### 6.4 DISCO data, estimated gene trees from 100 bp (`results/disco_comparison.md`)

Ten replicates per condition (`results/disco_comparison.md`, `results/disco_runs.jsonl`). Mean FN rate:

| condition | genes | ASTRID-multi | ASTRID-multi-w | ASTRID-Pro | ASTRID-Pro-min | ASTRID-Pro-w | ASTRID-DISCO | ASTRAL-Pro |
|---|---|---|---|---|---|---|---|---|
| default (dup 5e-10, loss/dup 1) | 1000 | 0.054 | 0.050 | 0.051 | 0.050 | 0.051 | **0.047** | 0.051 |
| gdl_1e-9_1 (dup 1e-9, loss/dup 1) | 1000 | 0.037 | 0.038 | 0.027 | **0.023** | 0.026 | 0.029 | 0.029 |
| gdl_1e-9_05 (dup 1e-9, loss/dup 0.5) | 100* | 0.089 | 0.086 | 0.046 | 0.063 | 0.044 | **0.042** | 0.051 |

\*Only the first 100 genes, because families average about 1000 leaves and the Python distance code is O(n²) per gene. gdl_1e-9_0 (about 3700 leaves per family) was skipped.

**Pre-registered held-out test** (`results/PREREG_disco_heldout.md`, `results/prereg_heldout.md`). The plan was written after seeing replicates 01–04 only:
- Set A: replicates 05–10 of the three conditions above.
- Set B: three conditions never run before (gdl_5e-10_05, ils_1e4, ils_2e8), replicates 01–05, 1000 genes.
- Primary method: ASTRID-Pro (mean). Holm correction over 3 comparisons.

| set | comparison | mean diff | W/T/L | p | Holm p |
|---|---|---|---|---|---|
| A ∪ B (n = 33) | ASTRID-Pro vs ASTRID-multi | −0.0096 | 16/12/5 | 0.007 | **0.021** |
| | vs ASTRID-DISCO | +0.0040 | 5/13/15 | 0.090 | 0.18 |
| | vs ASTRAL-Pro | +0.0028 | 8/11/14 | 0.30 | 0.30 |
| A only (n = 18) | vs ASTRID-multi | −0.0176 | 12/5/1 | 0.007 | 0.020 |
| B only (n = 15) | vs ASTRID-multi | 0.0000 | 4/7/4 | 0.57 | 0.66 |
| secondary, A ∪ B | ASTRID-Pro-w vs ASTRID-multi | −0.0118 | 16/15/2 | 0.0015 | – |

**Set B per-condition FN** (ASTRID-multi / ASTRID-Pro / ASTRID-DISCO / ASTRAL-Pro):

| condition | ASTRID-multi | ASTRID-Pro | ASTRID-DISCO | ASTRAL-Pro |
|---|---|---|---|---|
| gdl_5e-10_05 | 0.076 | 0.071 | 0.073 | 0.073 |
| ils_1e4 | 0.027 | 0.031 | 0.027 | 0.022 |
| ils_2e8 | 0.053 | 0.053 | 0.045 | 0.047 |

**Reading.** The orthology restriction helps where duplication is frequent and gene trees are accurate enough to tag. It does nothing at moderate GDL. ASTRID-DISCO, a different way to use orthology, is at least as good.

### 6.5 True gene trees with GDL + ILS, DISCO `gtrees_10000_l1` (up to 10,000 genes)

DISCO `gtrees_10000_l1`: 100 species, GDL + ILS, true gene trees, 5 replicates. Data: `results/disco_truetree_curve.jsonl` (full runs with ASTRAL-Pro) and `results/disco_truetree_stream.jsonl` (streaming runs, ASTRID methods only, for 10,000 genes).

Number of species-tree errors (FN out of 97) per replicate:

| genes | ASTRID-multi | ASTRID-Pro | ASTRID-DISCO | ASTRAL-Pro |
|---|---|---|---|---|
| 100 | 3, 2, 3, 6, 4 | 2, 0, 1, 4, 3 | 1, 1, 2, 3, 1 | 1, 0, 2, 5, 2 |
| 1000 | 0, 1, 0, 3, 1 | 0, 0, 0, 0, 0 | 0, 0, 0, 0, 0 | 0, 0, 0, 1, 0 |
| 3000 (streaming; reps 01, 03, 04, 05) | 0, 0, 2, 0 | 0, 0, 1, 0 | – | – |
| 10000 | 0, 1, 0, 0, 0 | 0, 0, 0, 0, 0 | 0 (rep 02 only) | 0 (rep 02 only) |

ASTRID-multi converges more slowly than the orthology-aware methods. At 10,000 genes it still has one error in replicate 02, but no plateau is evident across replicates.

### 6.6 Survival-reweighted correction (ASTRID-Pro-S; `code/correction.py`, `results/correction_small.jsonl`)

**How it works.**
- Map gene trees onto a rooted first-pass species tree.
- For each maximal gene subtree inside clade w, record whether its parent is a speciation node at parent(w). The fraction of "yes" estimates s(sibling of w), the survival probability of the other daughter.
- Count each visible speciation node on an orthologous path with weight 1/ŝ(off-path clade). The expected distance then becomes the species-tree node count, which is additive.

**Results on true gene trees, 20,000 families, true tags** (FN = species-tree errors):

| configuration | ASTRID-Pro (root counted) | ASTRID-multi | Pro-S, oracle first pass | Pro-S, first pass = ASTRID-multi | Pro-S iterated from its own tree |
|---|---|---|---|---|---|
| 4-taxon counterexample (two seeds) | 1, 1 | 0, 0 | **0, 0** | **0, 0** | 1 (stuck) |
| 5-taxon, supercritical x = (2, 0) | 1 | 0 | **0** | **0** | 1 (stuck) |

The estimated survival probabilities match theory: ŝ(C) = 0.049 vs 0.0498, and ŝ(x) = 0.68–0.70 vs the simulated s_yx = 0.695. Small upward bias for the cherry leaves (0.063 vs 0.050) comes from conditioning families on at least 2 species.

**Honest limit.** The correction needs a first pass that is right where the bias bites. Here ASTRID-multi was already right, so the correction has not yet been shown to add accuracy over its first pass. It was not tried on estimated gene trees.

## 7. Runtime

Single core, 100 species. Seconds per species tree, FastMulRFS data:

| step | 25 genes | 100 genes | 500 genes |
|---|---|---|---|
| root + tag (Python) | 0.22 | 0.91 | 4.45 |
| ASTRID-multi distances + FastME | 0.20 | 0.67 | 3.28 |
| ASTRID-Pro distances + FastME | 0.18 | 0.60 | 2.77 |
| ASTRID-DISCO (incl. root + tag) | 0.80 | 3.13 | 15.3 |
| ASTRAL-Pro3 | 3.44 | 13.4 | 87.7 |

The distance step is O(Σ n_g²) in pure Python/numpy. It became the bottleneck on the DISCO high-duplication conditions: about 1000 leaves per family for loss/dup 0.5, and 3728 for loss/dup 0, which I skipped. A C++ implementation, or ASTRID's own code with an ortholog mask, is needed for those.

## 8. Caveats

- My ASTRID-multi is a re-implementation. The official ASTRID-2 binary would need a bazel/JNI build, which I did not do, and the published FastMulRFS-data ASTRID trees are flagged buggy. So the ASTRID-multi numbers are not validated against the original program.
- DISCO is re-implemented (repository unreachable). The original decomposition rule may differ in tie handling.
- I ran ASTRAL-Pro3 (ASTER v1.25), not the 2022 ASTRAL-Pro 2 build. It is the same algorithm family and probably at least as accurate.
- On the DISCO data, the variant was chosen after seeing replicates 01–04. The held-out test (replicates 05–10 plus three new conditions) was pre-registered, but its runs for set A had already been computed (not inspected) when the plan was written. The ASTRID-Pro variant that wins differs between datasets: closest-copy was picked on the FastMulRFS dev split, the mean variant on DISCO data.
- The data are SimPhy conditions with modest GDL. On the FastMulRFS data, ASTRAL-type methods are not at their best with only 25–500 short-sequence genes.

## 9. A 4-week project, and risks

**Plan.**
- **Week 1: theory.**
  - Write the additivity theorem with full proof: root counting, NJ/FastME consistency, and the effect of conditioning on at least k species.
  - Make the supercritical counterexample analytic with closed-form birth–death survival.
  - State and prove the corrected distance, with consistency given a correct first pass. Possibly also a two-stage version whose first pass is any consistent method.
- **Week 2: ASTRID-multi.**
  - Exact expected per-gene-averaged distances on 4–5 taxa via generating functions, or 10⁸-family C++ Monte Carlo as a fallback.
  - Either a counterexample (answers the slide question negatively) or a conjecture with strong numerical support.
- **Week 3: implementation and benchmarks.**
  - Port the orthology-restricted distance into ASTRID's C++, an ortholog mask plus speciation counts, so the 1000- and 3700-leaf families become feasible.
  - Rerun the DISCO high-GDL conditions, including gdl_1e-9_0, and the GDL-comparison data. Validate against the real ASTRID-multi binary.
- **Week 4.** Write-up.

**Risks.**
- **Accuracy.** The gain over ASTRID-multi is confirmed only in high-GDL conditions, and ASTRID-Pro does not beat ASTRID-DISCO or ASTRAL-Pro. Frame it as "repairs ASTRID-multi", not "new best method".
- **ASTRID-multi may stay open.** Exact expectations of per-gene *averages* are ratio expectations, which may be hard. Fallback: analyse the pooled-pairs variant, which is linear, and report per-gene averaging as open.
- **"Too easy".** The instructor may find the theorem folklore-easy. Mitigation: the root-counting subtlety, the supercritical counterexample and the reconciliation-based correction are concrete, non-obvious results.
- **Estimated tags.** Tagging error and hidden paralogy are not covered by the theory. A sibling pilot (`claude/cs581-gdlcons`) studies ASTRAL-Pro under random rooting and tagging error and could share that model.
