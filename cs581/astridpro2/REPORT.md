# ASTRID-Pro, round 2: theorem, broader benchmark, scaling

*CS581 pilot, 2026-10-10. Branch `claude/cs581-astridpro2`. Code is in `code/`, raw results in `results/`, the
theorem in `theory.md`, and the pre-registration in `PREREG.md`. Machine: 4 cores and 15 GB RAM. All methods ran
single-threaded, but up to 4 jobs ran at once (see the runtime caveat).*

## 0. Verdict

**Promising for a 4-week CS581 project. Go, by the pre-registered criterion.** The criterion was "matches
ASTRAL-Pro3 and is at least 10× faster at 1000 genes × 100 taxa". Both parts are met, but the accuracy gain is over
ASTRAL-Pro3, not over the other fast distance methods.

1. **Theory (new, clean, checked).**
   - Under GDL with ideal rooting and tagging, ASTRID-Pro's limit is an **additive metric on the species tree**
     with explicit edge lengths β(c).
   - It is **consistent iff every interior β > 0**.
   - That holds whenever no interior branch is supercritical, and it holds on **every standard simulated GDL
     benchmark**: exact min β ≥ 0.46.
   - An exact supercritical counterexample (β = −0.267) converges to the wrong tree.
   - A survival-reweighted variant (ASTRID-Pro-S) is consistent under any rates, given a correct first pass.
2. **Accuracy on estimated gene trees.** Pooled simulated data from FastMulRFS (100 taxa) and DISCO (100 taxa,
   12 conditions). ASTRID-Pro FN rate minus each method (negative favours ASTRID-Pro):
   - vs **ASTRAL-Pro3**: **−0.0135 FN rate, 84 wins / 30 ties / 26 losses, Holm p = 6e-8** (n = 140: 120 FastMulRFS + 20 DISCO pairs). On DISCO alone: −0.0041, 8/9/3, p = 0.054, and ASTRAL-Pro3 also timed out on 4 of 24 DISCO runs. Our ASTRAL-Pro3 matches the authors' published ASTRAL-Pro trees on the same 120
     inputs (0.0890 vs 0.0904, p = 0.5), so the baseline is sound.
   - vs **ASTRID-multi**: −0.0034, 75/72/53, p = 0.023 (Holm 0.068, so *not* significant after correction). The gain is concentrated in the high-duplication DISCO conditions (dup 1e-9, or dup 5e-10 with
     no loss): −0.014, 20/4/3, p = 0.003.
   - **Ties** with ASTRID-DISCO (+0.0009, p = 0.32) and Asteroid (+0.0002, p = 0.98).
   - Beats DISCO+ASTRAL, FastMulRFS and DupLoss-2. Ties wQFM-GDL (the 2026 state of the art) on the FastMulRFS data: −0.004, p = 0.09, 54 pairs.
   - Empirical: 1KP C12 has 5 FN vs the 1KP ASTRAL reference, the same as ASTRAL-Pro3, against 10 for ASTRID-multi.
     Fungi16 and vertebrates188: all methods agree.
3. **Speed.**
   - 0.2 s (100 genes) and 1–40 s (1000 genes, up to 3,700-leaf families) on 100 taxa, against 2.5–28 min for
     ASTRAL-Pro3 on DISCO (1000 genes, 100 taxa, 1 thread, shared machine).
   - ASTRAL-Pro3 timed out (> 40 min) on 4 of 24 runs in the 3 heaviest DISCO conditions.
   - At 1000 taxa: 14 s. ASTRAL-Pro3 and FastMulRFS already exceed 20 min at 500 taxa.
   - 10,000 genes: 8 s. ASTRAL-Pro3 needs 15 min at 3,000 genes and exceeds 20 min at 10,000.
   - About **100–500× faster than ASTRAL-Pro3**, with memory independent of the number of genes (streaming).

**The honest framing.**
- ASTRID-Pro is a **very fast distance method, provably consistent under pure GDL with ideal rooting and
  tagging, that matches or beats ASTRAL-Pro3** in these benchmarks.
- It is *not* better than ASTRID-DISCO or Asteroid in accuracy. They are fast too, but neither has a GDL consistency
  proof.
- The project's novelty is the theorem (an answer to slide 35 for the orthology-restricted distance) plus a fast
  implementation. It does not answer the question for ASTRID-multi.


## 1. Theory (full statement and proofs in `theory.md`)

**Setting.** Pure GDL (any branch-specific, time-varying birth–death rates; no ILS), true gene trees, true root and
D/S tags. The ASTRID-Pro distance counts the speciation nodes on the path between orthologous copies, including the
LCA and the gene-tree root.

- **Thm 3 (exact limit).** The averaged matrix converges a.s. to a metric that is *additive on the species tree S*.
  Internal edge c has length

  β(c) = ½ [ s(c₁) + s(c₂) + s(sib c) − s(c) ]

  where s(·) is the survival probability of a copy entering a branch. These probabilities have closed forms for
  linear birth–death. Leaf edges are ½[1 + s(sib) − s(leaf)] ≥ 0.
- **Thm 4 (iff).**
  - If β(c) > 0 at every interior node, ASTRID-Pro + NJ/FastME is statistically consistent.
  - If some β(c) < 0, the limit violates the four-point condition. The true split around c becomes the *largest*
    quartet sum, so every distance method fails on that quartet.
  - This corrects the previous note, which said the limit is "always a tree metric".
- **Thm 5 (sufficient).** β(c) > 0 whenever the branch above c is not supercritical (λ ≤ μ). Leaf branches and root
  branches can have any rates. With no loss, β ≡ 1.
- **Counterexample.** Exact β = −0.267. Simulation agrees with the predicted matrix to 0.003 at 100,000 families,
  and ASTRID-Pro converges to the wrong tree. A second, random 6-taxon instance with β = −0.034 is also wrong at
  20,000 families.
- **Benchmarks are covered.** The exact minimum interior β on every DISCO and FastMulRFS simulation condition is
  ≥ 0.46, including the supercritical loss/dup = 0.5 conditions. Over 20,000 random-rate configurations on 4–8
  taxa, 1.15% fail, and all of those have a supercritical branch.
- **Prop 6 (correction).** Reweighting each counted node by 1/ŝ(off-path clade), with ŝ estimated by reconciling to
  a correct rooted first-pass tree, gives a limit with unit edge lengths. That is consistent under *any* GDL rates.
  It fixes both failing instances given a correct first pass.
- **Open.** Estimated tags (hidden paralogy), MinDup rooting, ILS, and ASTRID-multi itself. Counting the gene-tree
  root is required by the theorem with true roots. It is harmless or slightly harmful on estimated trees; see
  `astrid-pro-r0`.

## 2. Data

| dataset | source | used here |
|---|---|---|
| FastMulRFS data (= ASTRAL-Pro "S100") | Molloy & Warnow 2020, doi:10.13012/B2IDB-5721322_V1 | 100 species. DL rate 1e-10, 2e-10, 5e-10 × Ne 1e7, 5e7. RAxML gene trees from 25 and 100 bp. 100 genes, reps 01–10 (120 runs). The planned 500-gene runs were not reached. |
| DISCO data | Willson et al. 2022, doi:10.13012/B2IDB-4050038_V1 | 100 species, 1000 estimated gene trees (100 bp). 12 conditions: default; gdl_{1e-10, 5e-10, 1e-9}_{0, 0.5, 1} (families up to 3,700 leaves); ils_1e4; ils_2e8; missing_1000. Fast methods on reps 01–07 (80 runs); ASTRAL-Pro3 on reps 01–02. |
| scaling | DISCO species_1000 rep 01 (1000 species) induced on nested random subsets of 50–1000 species (~400 genes); DISCO gtrees_10000_l1 (100 species, 100 to 10,000 estimated genes, reps 01–02) | runtime and FN |
| empirical | fungi16 (7,180 trees; MrBayes reference) and 1KP C12 (9,237 trees; 83 species; reference = original 1KP ASTRAL tree, 80 shared species), both from the ASTRAL-Pro repository; vertebrates188 (31,612 trees; reference with polytomies) from the DISCO data | FN vs reference |
| **ASTRAL-Pro's own SimPhy gene trees** (Zhang et al. 2020, Dryad doi:10.6076/D11C7P) | **not used**: Dryad returns 403/401 from this machine | The DISCO grid spans the same axes (duplication rate × loss/dup ratio, and ILS). The FastMulRFS data *is* the ASTRAL-Pro S100 data. |

## 3. Methods (all single-threaded; `code/bench.py`)

- **ASTRID-Pro** (`code/apro.cpp`). The theorem's estimator.
  - Rooting: MinDup, with ties broken by fewest losses.
  - Tagging: species overlap.
  - Distance: orthologous pairs only, speciation nodes only, root counted.
  - Averaging: per-gene mean, then mean over genes.
  - Tree: FastME 2 (BalME + NNI + SPR).
  - Implementation: C++, O(Σ n_g·k_g) by aggregating at each LCA, streaming over genes.
  - Validation: identical (to 5e-7) to the previous pilot's Python implementation.
- **Variants.**
  - `-r0`: root not counted.
  - `-S`: survival-reweighted second pass. The first pass is ASTRID-Pro, rooted by fewest total duplications.
- **ASTRID-multi**: same code, all pairs and all nodes (ASTRID's averaging).
- **ASTRID-DISCO**: official DISCO v1.4.1 decomposition, then ASTRID.
- **DISCO+ASTRAL**: official DISCO, then ASTRAL-IV (ASTER v1.25).
- **ASTRAL-Pro3**: ASTER v1.25.3.8.
- **Asteroid**: built from GitHub, with its missing-data correction. It crashes on trees with < 4 leaves, which are
  filtered out for it, and on some high-duplication inputs.
- **FastMulRFS**: built from GitHub with FastRFS; the "single" tree.
- **wQFM-GDL**: v1.0.3 jar, tree mode.
- **DupLoss-2**: GitHub build.
- **SpeciesRax was not run**: it needs MSAs and GeneRax family files, and is too costly here.

Metric: species-tree FN rate. Paired by (dataset, condition, replicate, sequence length, #genes). Two-sided Wilcoxon
signed-rank and 95% bootstrap CI of the mean difference. Holm correction over the 4 pre-registered primary
comparisons (`PREREG.md`).

## 4. Accuracy results (full tables: `results/tables.md`)

**Check of the ASTRAL-Pro baseline** (`code/check_published_apro.py`). On the same 120 FastMulRFS inputs, the authors' published ASTRAL-Pro (v2) trees score 0.0904 and our ASTRAL-Pro3 scores 0.0890 (p = 0.5). ASTRID-Pro scores 0.0740 and beats the *published* trees by −0.016 (79/0/41, p = 5.5e-7).

### 4.1 Pooled simulated data (diff = ASTRID-Pro − method, FN rate; W = ASTRID-Pro better)

**Primary comparisons (Holm over 4)**

| ASTRID-Pro vs | n | mean diff (FN rate) | 95% CI | W/T/L | p (Wilcoxon) | Holm p |
|---|---|---|---|---|---|---|
| astral-pro3 | 140 | -0.0135 | [-0.0179, -0.0092] | 84/30/26 | 1.4e-08 | 5.6e-08 |
| astrid-multi | 200 | -0.0034 | [-0.0060, -0.0009] | 75/72/53 | 0.023 | 0.068 |
| astrid-disco | 199 | +0.0009 | [-0.0013, +0.0031] | 49/86/64 | 0.32 | 0.64 |
| asteroid | 188 | +0.0002 | [-0.0024, +0.0027] | 58/68/62 | 0.98 | 0.98 |

**Secondary comparisons (no correction)**

| ASTRID-Pro vs | n | mean diff (FN rate) | 95% CI | W/T/L | p (Wilcoxon) |
|---|---|---|---|---|---|
| astrid-pro-r0 | 200 | +0.0003 | [-0.0003, +0.0010] | 9/179/12 | 0.75 |
| astrid-pro-s | 200 | -0.0008 | [-0.0019, +0.0002] | 24/162/14 | 0.23 |
| disco-astral | 121 | -0.0235 | [-0.0302, -0.0172] | 79/22/20 | 7.5e-11 |
| fastmulrfs | 121 | -0.0240 | [-0.0299, -0.0184] | 87/18/16 | 7.8e-13 |
| wqfm-gdl | 54 | -0.0042 | [-0.0088, +0.0002] | 22/14/18 | 0.092 |
| duploss2 | 54 | -0.1075 | [-0.1275, -0.0888] | 53/1/0 | 2.3e-10 |

### 4.2 FastMulRFS data

**FastMulRFS data: ASTRID-Pro vs each method**

| ASTRID-Pro vs | n | mean diff (FN rate) | 95% CI | W/T/L | p (Wilcoxon) |
|---|---|---|---|---|---|
| astrid-pro-r0 | 120 | +0.0003 | [-0.0007, +0.0014] | 8/103/9 | 0.78 |
| astrid-pro-s | 120 | -0.0014 | [-0.0032, +0.0003] | 24/82/14 | 0.23 |
| astrid-multi | 120 | -0.0024 | [-0.0059, +0.0012] | 43/40/37 | 0.13 |
| astrid-disco | 120 | +0.0018 | [-0.0013, +0.0050] | 31/46/43 | 0.3 |
| asteroid | 120 | +0.0013 | [-0.0022, +0.0049] | 39/36/45 | 0.59 |
| astral-pro3 | 120 | -0.0150 | [-0.0201, -0.0101] | 76/21/23 | 5.2e-08 |
| disco-astral | 119 | -0.0237 | [-0.0307, -0.0172] | 77/22/20 | 7.9e-11 |
| fastmulrfs | 119 | -0.0246 | [-0.0306, -0.0189] | 87/17/15 | 4.9e-13 |
| wqfm-gdl | 54 | -0.0042 | [-0.0090, +0.0004] | 22/14/18 | 0.092 |
| duploss2 | 54 | -0.1075 | [-0.1279, -0.0878] | 53/1/0 | 2.3e-10 |

**FastMulRFS data: mean FN rate by sequence length and #genes**

| sqln | ngen | n | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl | duploss2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 100 | 60 | 0.042 | 0.042 | 0.042 | 0.046 | 0.041 | 0.040 | 0.048 | 0.050 (59) | 0.053 | **0.040 (27)** | 0.087 (27) |
| 25 | 100 | 60 | 0.106 | 0.105 | 0.109 | 0.107 | **0.104** | 0.105 | 0.130 | 0.145 | 0.144 (59) | 0.114 (27) | 0.273 (27) |

**FastMulRFS data: mean FN rate by condition (all sqln, ngen)**

| cond | n | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl | duploss2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dl1e-10_ps1e7 | 20 | 0.047 | 0.046 | 0.047 | 0.050 | 0.047 | 0.046 | 0.057 | 0.060 | 0.058 | **0.043 (10)** | 0.156 (10) |
| dl1e-10_ps5e7 | 20 | 0.069 | 0.069 | 0.067 | 0.068 | **0.066** | 0.070 | 0.070 | 0.074 | 0.081 | 0.068 (10) | 0.193 (10) |
| dl2e-10_ps1e7 | 20 | 0.065 | 0.064 | 0.069 | 0.065 | 0.064 | **0.063** | 0.076 | 0.082 | 0.089 | 0.069 (10) | 0.178 (10) |
| dl2e-10_ps5e7 | 20 | 0.070 | 0.070 | 0.071 | 0.069 | 0.069 | **0.066** | 0.084 | 0.087 (19) | 0.091 (19) | 0.079 (8) | 0.197 (8) |
| dl5e-10_ps1e7 | 20 | 0.090 | 0.092 | 0.094 | 0.099 | 0.093 | **0.086** | 0.121 | 0.136 | 0.128 | 0.110 (8) | 0.164 (8) |
| dl5e-10_ps5e7 | 20 | 0.103 | 0.102 | 0.105 | 0.108 | **0.094** | 0.105 | 0.126 | 0.148 | 0.142 | 0.104 (8) | 0.196 (8) |

### 4.3 DISCO data

**DISCO data: ASTRID-Pro vs each method**

| ASTRID-Pro vs | n | mean diff (FN rate) | 95% CI | W/T/L | p (Wilcoxon) |
|---|---|---|---|---|---|
| astrid-pro-r0 | 80 | +0.0003 | [-0.0003, +0.0008] | 1/76/3 | 0.75 |
| astrid-pro-s | 80 | +0.0000 | [+0.0000, +0.0000] | 0/80/0 | 1 |
| astrid-multi | 80 | -0.0048 | [-0.0085, -0.0015] | 32/32/16 | 0.016 |
| astrid-disco | 79 | -0.0005 | [-0.0036, +0.0022] | 18/40/21 | 0.83 |
| asteroid | 68 | -0.0017 | [-0.0045, +0.0011] | 19/32/17 | 0.27 |
| astral-pro3 | 20 | -0.0041 | [-0.0087, -0.0000] | 8/9/3 | 0.054 |
| disco-astral | 2 | -0.0102 | [-0.0102, -0.0102] | 2/0/0 | 0.5 |
| fastmulrfs | 2 | +0.0102 | [+0.0000, +0.0204] | 0/1/1 | 1 |

Per-condition means are below. Columns with fewer replicates (count in parentheses) are not paired with the others; use the paired table above.

**DISCO data: mean FN rate by condition (100 bp, 1000 genes)**

| cond | n | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs |
|---|---|---|---|---|---|---|---|---|---|---|
| default | 7 | 0.042 | 0.042 | 0.042 | 0.045 | 0.036 | 0.042 | **0.036 (2)** | 0.071 (1) | 0.061 (1) |
| gdl_1e-10_0 | 7 | 0.041 | 0.041 | 0.041 | 0.042 | 0.042 | **0.038** | 0.082 (2) | – | – |
| gdl_1e-10_05 | 7 | 0.044 | 0.042 | 0.044 | 0.044 | 0.044 | **0.041** | 0.046 (2) | – | – |
| gdl_1e-10_1 | 7 | 0.028 | 0.028 | 0.028 | **0.025** | 0.028 | 0.026 | 0.051 (2) | 0.071 (1) | 0.041 (1) |
| gdl_1e-9_0 | 6 | **0.058** | **0.058** | **0.058** | 0.092 | 0.077 | – | – | – | – |
| gdl_1e-9_05 | 7 | 0.031 | 0.029 | 0.031 | 0.039 | 0.029 (6) | 0.037 (3) | **0.000 (1)** | – | – |
| gdl_1e-9_1 | 7 | 0.034 | 0.034 | 0.034 | 0.045 | **0.031** | 0.047 | 0.051 (2) | – | – |
| gdl_5e-10_0 | 7 | 0.022 | 0.022 | 0.022 | 0.028 | 0.026 | 0.024 (5) | **0.020 (1)** | – | – |
| gdl_5e-10_05 | 7 | **0.054** | **0.054** | **0.054** | 0.060 | **0.054** | 0.063 | 0.056 (2) | – | – |
| ils_1e4 | 6 | 0.037 | 0.037 | 0.037 | **0.032** | 0.032 | 0.034 | 0.036 (2) | – | – |
| ils_2e8 | 6 | 0.053 | 0.053 | 0.053 | 0.054 | **0.048** | 0.048 | 0.051 (2) | – | – |
| missing_1000 | 6 | 0.039 | 0.039 | 0.039 | **0.034** | 0.043 | **0.034** | 0.041 (2) | – | – |

Reading: ASTRID-Pro improves on ASTRID-multi exactly where the theory says orthology restriction matters, the high-duplication conditions (gdl_1e-9_*: 0.023–0.061 vs 0.036–0.092). It is no better under pure ILS (ils_*) or missing data. ASTRID-Pro-S never differs from ASTRID-Pro on DISCO, which fits β > 0 on all of these conditions.

### 4.4 Empirical

**FN rate vs reference tree**

| cond | n | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1kp_c12 | 1 | **0.072** | **0.072** | 0.087 | 0.145 | 0.087 | – | **0.072** | **0.072** | – | – |
| fungi16 | 1 | **0.077** | **0.077** | **0.077** | **0.077** | **0.077** | **0.077** | **0.077** | **0.077** | **0.077** | **0.077** |
| vertebrates188 | 1 | **0.038** | **0.038** | **0.038** | **0.038** | **0.038** | – | – | – | – | – |

## 5. Runtime and scaling (seconds; 1 thread per method; shared 4-core machine)

**FastMulRFS data: mean runtime (s, 1 thread)**

| ngen | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl | duploss2 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 100 | 0.15 | 0.16 | 0.53 | 0.17 | 1.01 | 1.86 | 10.02 | 14.82 | 20.20 | 73.70 | 11.60 |

**DISCO data: mean runtime (s, 1 thread)**

| cond | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs |
|---|---|---|---|---|---|---|---|---|---|
| default | 0.82 | 0.83 | 3.07 | 0.89 | 8.07 | 11.01 | 384.71 | 891.78 | 404.98 |
| gdl_1e-10_0 | 0.83 | 0.85 | 3.11 | 1.03 | 6.91 | 1.03 | 556.47 | – | – |
| gdl_1e-10_05 | 0.86 | 0.88 | 2.89 | 0.88 | 6.18 | 8.93 | 458.58 | – | – |
| gdl_1e-10_1 | 0.69 | 0.76 | 2.75 | 0.90 | 5.67 | 11.91 | 308.91 | 636.33 | 726.43 |
| gdl_1e-9_0 | 41.12 | 39.98 | 168.97 | 162.71 | 294.54 | – | – | – | – |
| gdl_1e-9_05 | 8.72 | 8.91 | 40.42 | 66.07 | 78.14 | 25.02 | 586.02 | – | – |
| gdl_1e-9_1 | 1.36 | 1.39 | 5.05 | 3.43 | 11.93 | 13.84 | 565.49 | – | – |
| gdl_5e-10_0 | 3.61 | 3.60 | 13.10 | 9.53 | 31.54 | 12.39 | 783.47 | – | – |
| gdl_5e-10_05 | 1.49 | 1.54 | 5.61 | 2.84 | 14.54 | 16.64 | 1028.41 | – | – |
| ils_1e4 | 0.94 | 0.93 | 3.57 | 0.91 | 8.49 | 10.19 | 330.15 | – | – |
| ils_2e8 | 0.88 | 0.92 | 3.51 | 1.12 | 8.81 | 11.23 | 355.29 | – | – |
| missing_1000 | 0.50 | 0.53 | 1.86 | 0.67 | 4.37 | 3.75 | 129.99 | – | – |

**Runtime (s) vs #species (species_1000 rep 01, ~400-435 genes)**

| cond | astrid-pro | astrid-pro-s | astrid-multi | astrid-disco | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl |
|---|---|---|---|---|---|---|---|---|
| taxa_100 | 0.59 | 1.74 | 1.35 | 5.02 | 84.60 | 81.10 | 36.27 | 349.69 |
| taxa_1000 | 14.32 | 90.56 | 135.73 | 670.66 | – | – | – | – |
| taxa_200 | 1.49 | 6.19 | 2.61 | 9.17 | 541.56 | 485.31 | 175.77 | – |
| taxa_50 | 0.23 | 0.58 | 0.29 | 1.74 | 14.08 | 16.23 | 11.52 | 66.64 |
| taxa_500 | 3.63 | 26.62 | 42.11 | 48.26 | – | – | – | – |

**FN rate vs #species**

| cond | n | astrid-pro | astrid-pro-s | astrid-multi | astrid-disco | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl |
|---|---|---|---|---|---|---|---|---|---|
| taxa_100 | 1 | 0.082 | 0.082 | 0.093 | 0.082 | 0.082 | 0.093 | **0.072** | **0.072** |
| taxa_1000 | 1 | **0.066** | **0.066** | 0.085 | 0.077 | – | – | – | – |
| taxa_200 | 1 | **0.061** | **0.061** | 0.086 | 0.076 | 0.071 | 0.081 | 0.091 | – |
| taxa_50 | 1 | 0.085 | 0.085 | 0.106 | 0.064 | **0.043** | **0.043** | **0.043** | 0.064 |
| taxa_500 | 1 | 0.064 | 0.064 | 0.068 | **0.060** | – | – | – | – |

**Runtime (s) vs #genes (gtrees_10000_l1, 100 species)**

| ngen | astrid-pro | astrid-pro-s | astrid-multi | astrid-disco | astral-pro3 | fastmulrfs |
|---|---|---|---|---|---|---|
| 100 | 0.17 | 0.82 | 0.23 | 1.11 | 12.48 | 9.27 |
| 300 | 0.37 | 1.54 | 0.64 | 3.55 | 53.59 | 29.05 |
| 1000 | 1.18 | 4.47 | 1.27 | 9.51 | 340.44 | 136.85 |
| 3000 | 2.82 | 11.33 | 3.14 | 28.37 | 904.49 | – |
| 10000 | 8.00 | 36.16 | 9.18 | 95.29 | – | – |

**FN rate vs #genes**

| ngen | n | astrid-pro | astrid-pro-s | astrid-multi | astrid-disco | astral-pro3 | fastmulrfs |
|---|---|---|---|---|---|---|---|
| 100 | 2 | **0.051** | **0.051** | 0.056 | 0.071 | 0.066 | 0.168 |
| 1000 | 2 | 0.036 | 0.036 | 0.041 | **0.031** | 0.041 | 0.046 |
| 10000 | 2 | **0.031** | **0.031** | 0.036 | **0.031** | – | – |
| 300 | 2 | 0.041 | 0.041 | **0.036** | 0.051 | 0.061 | 0.087 |
| 3000 | 2 | 0.036 | 0.036 | 0.031 | 0.036 | **0.020 (1)** | – |

**Empirical runtime (s)**

| cond | astrid-pro | astrid-pro-r0 | astrid-pro-s | astrid-multi | astrid-disco | asteroid | astral-pro3 | disco-astral | fastmulrfs | wqfm-gdl |
|---|---|---|---|---|---|---|---|---|---|---|
| 1kp_c12 | 6.89 | 6.69 | 26.28 | 6.83 | 46.92 | – | 3347.17 | 2962.00 | – | – |
| fungi16 | 0.53 | 0.53 | 1.22 | 0.69 | 2.88 | 6.02 | 11.55 | 23.05 | 8.94 | 17.07 |
| vertebrates188 | 21.51 | 22.54 | 79.32 | 22.98 | 168.53 | – | – | – | – | – |

**Failed runs (timeouts and crashes; excluded from all paired tables)**

| file | method | condition | error | n |
|---|---|---|---|---|
| disco_runs.jsonl | asteroid | gdl_1e-9_0 | CalledProcessError | 6 |
| disco_runs.jsonl | asteroid | gdl_1e-9_05 | CalledProcessError | 3 |
| disco_runs.jsonl | asteroid | gdl_5e-10_0 | CalledProcessError | 1 |
| disco_runs.jsonl | astral-pro3 | gdl_1e-9_0 | TimeoutExpired | 2 |
| disco_runs.jsonl | astral-pro3 | gdl_1e-9_05 | TimeoutExpired | 1 |
| disco_runs.jsonl | astral-pro3 | gdl_5e-10_0 | TimeoutExpired | 1 |
| disco_runs.jsonl | wqfm-gdl | gdl_1e-10_1 | CalledProcessError | 1 |
| empirical_runs.jsonl | asteroid | 1kp_c12 | CalledProcessError | 1 |
| empirical_runs.jsonl | asteroid | vertebrates188 | CalledProcessError | 1 |
| empirical_runs.jsonl | fastmulrfs | 1kp_c12 | CalledProcessError | 1 |
| empirical_runs.jsonl | wqfm-gdl | 1kp_c12 | TimeoutExpired | 1 |
| fmrfs_runs.jsonl | disco-astral | dl2e-10_ps5e7 | CalledProcessError | 1 |
| fmrfs_runs.jsonl | fastmulrfs | dl2e-10_ps5e7 | CalledProcessError | 1 |
| scaling_runs.jsonl | asteroid | taxa_200 | CalledProcessError | 1 |
| scaling_runs.jsonl | astral-pro3 | gtrees_10000_l1 | TimeoutExpired | 2 |
| scaling_runs.jsonl | astral-pro3 | taxa_1000 | TimeoutExpired | 1 |
| scaling_runs.jsonl | astral-pro3 | taxa_500 | TimeoutExpired | 1 |
| scaling_runs.jsonl | fastmulrfs | gtrees_10000_l1 | TimeoutExpired | 2 |
| scaling_runs.jsonl | fastmulrfs | taxa_1000 | TimeoutExpired | 1 |
| scaling_runs.jsonl | fastmulrfs | taxa_500 | TimeoutExpired | 1 |



## 6. Caveats

- **Runtimes were measured on a shared machine.** Up to 4–5 single-threaded jobs ran at once on 4 cores (wQFM-GDL
  is multi-threaded Java). Absolute seconds are inflated, perhaps up to 2×, for every method alike. The ratios
  (100–500×) are far larger than this noise.
- **ASTRAL-Pro3 comparisons on DISCO cover reps 01–02 only** (20 finished pairs). Three high-duplication conditions
  are missing because ASTRAL-Pro3 timed out at 40 min. The pooled ASTRAL-Pro3 comparison is therefore driven by the
  FastMulRFS data (120 pairs).
- **Not fully held out.** The previous pilot ran an ASTRID-Pro variant on the FastMulRFS data and on 6 of the 12
  DISCO conditions. This round tuned nothing, because the variant was fixed by the theorem before any run
  (`PREREG.md`). The queue was trimmed twice, for throughput only, before the affected comparisons were looked at.
  Notes are in `code/mkjobs_round*.py`.
- DupLoss-2 and wQFM-GDL ran on the FastMulRFS data at 100 genes (reps 01–04) only. FastMulRFS and DISCO+ASTRAL ran
  on DISCO only for 1–2 replicates.
- Asteroid failed on the largest families (gdl_1e-9_0, gdl_1e-9_05) and on two empirical sets.
- **The vertebrates188 reference has polytomies**, so only FN is meaningful there.
- **Theory gaps.** The theory assumes true tags and roots. Estimated tags (hidden paralogs) and ILS are open. Prop 6
  needs a correct first pass.

## 7. A 4-week project, and risks

- **Week 1: theory write-up.**
  - Thm 3–5 and Prop 6 with full proofs (done in draft).
  - Extend Lemma 1 to MinDup-rooted, overlap-tagged *true* gene trees: bound the effect of hidden duplications.
  - Prove or refute: no ILS + no hidden duplications ⇒ MinDup gives the true root.
  - State and test the counting-the-root rule under ILS.
- **Week 2: ILS.**
  - Under DLCoal, characterise when an orthologous pair's expected S-node count stays additive.
  - At minimum, prove the two limits: GDL only (done) and MSC only (= NJst with root not counted).
  - Run the exact-β calculator with coalescent depth added. Simulate DLCoal counterexamples on 4–5 taxa with
    SimPhy/gdlsim.
- **Week 3: experiments.**
  - Finish the ASTRAL-Pro3 runs on all DISCO replicates, using a multi-threaded ASTRAL-Pro3 to get them done.
  - Clean single-job timing on an idle machine.
  - Add SpeciesRax / AleRax where alignments exist.
  - Add an adversarial benchmark with rate-heterogeneous branches (β < 0), where ASTRID-Pro should fail, ASTRID-Pro-S
    should recover, and ASTRAL-Pro should be fine.
- **Week 4: write-up.**

**Risks.**
1. **The accuracy story is "as good as ASTRAL-Pro, much faster", not "best".** ASTRID-DISCO and Asteroid are equally
   accurate and also fast. The defensible claim is the consistency proof plus speed and scalability.
2. **The theorem may be judged too easy** once written: it is the branching property plus Buneman. Mitigation: the
   iff condition with explicit β, the certification of the standard benchmarks, the counterexample, and the
   corrected estimator are concrete and non-obvious.
3. **ILS extension may not close.** Fallback: present GDL-only and MSC-only as proved, and GDL + ILS as empirical.
4. **ASTRID-multi stays open.** Say so explicitly. Do not claim to answer the slide question for ASTRID-multi.
5. **Compute.** ASTRAL-Pro3 at 1000 genes takes 5–40 min per run here. Plan multi-threaded runs, or a campus cluster,
   for the baselines.

