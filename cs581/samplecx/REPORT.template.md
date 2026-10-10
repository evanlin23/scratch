# Pilot: sample complexity of ASTRID/NJst under the MSC, and a missing-data-corrected ASTRID

This pilots one candidate project for CS581 (Fall 2026). The branch is `claude/cs581-samplecx`.

- **Code:** `code/`
- **Small result files and plots:** `results/`
- **Literature notes with DOIs:** `prior_art.md`
- **Raw data and per-run logs:** `/opt/data/nute` and `/opt/runs` (not committed)

All numbers below come from the scripts named in each section, so each can be regenerated.

## 0. Verdict up front

| Sub-idea | Verdict | Why |
|---|---|---|
| **A. Sample complexity of ASTRID/NJst vs ASTRAL** (slide question) | **Promising, as an evaluation and theory-guided study** | **What the pilot shows:** the f-dependence is about the same for all three methods. The n-dependence is not. On caterpillar trees ASTRID/NJst need more and more genes than ASTRAL as n grows: at n = 64 and f = 0.1 CU, k₉₅ ≈ 1190 vs ≈ 490 (30 replicates). In the n = 128 test at f = 0.2 CU (20 replicates), k₉₅ is about 350–400 for ASTRID vs about 150–285 for ASTRAL. On balanced trees the gap is smaller and shows up only at n ≥ 64 (n = 64, f = 0.1: 771 vs 467). On the 25-taxon Nute et al. data the three methods are tied. **Why it is worth pursuing:** this matches the open conjecture in Roch (2018), that internode-distance methods need ≥ linear-in-n genes while ASTRAL needs O(log n), and nobody has measured it beyond n = 8. **Caveat:** the evidence is suggestive, not yet conclusive (20–30 replicates per cell; several bootstrap CIs overlap). |
| **B. Missing-data-corrected ASTRID** (new method) | **Not promising** | **The derivation works:** I derived the unique unbiased per-node correction, 1 − (−p/(1−p))^k, which makes ASTRID statistically consistent under i.i.d. deletion. **But it does not help empirically.** On held-out replicates, no bounded variant beats plain ASTRID (Δ FN = +0.0001, p = 0.71). The exact correction blows up for p ≥ 0.5 (+0.33 FN). The ASTRID bias it removes is small: a perfect correction gains at most about 0.005–0.009 FN in the infinite-gene limit. **The niche is already filled:** Asteroid (Morel et al. 2023), a consistent missing-data-aware internode method, is already published, and it is the only method here that beats ASTRID (−0.008 FN, p = 0.0035, RAxML gene trees). |

**Recommendation.** If this topic is chosen, do **A** and frame it as "an empirical test of Roch's n-scaling conjecture for ASTRID". **B** can be one short section (the derivation, plus "the bias is real but small; Asteroid already addresses it").

## 1. The questions and where they come from

- **Slide question.** The slide is *Phylogenomics part 2, slide 35, "(Some) Open Questions"*: "What is the sample complexity for ASTRID and NJst under the MSC?" (listed as item 8 in `cs581/literature/slide_open_problems.md`).
- **Missing-data question.** Rhodes, Nute & Warnow (arXiv:2001.07844) show that ASTRID/NJst are statistically **inconsistent** under i.i.d. taxon deletion; the published record is the 2020 correction to Nute et al. 2018, doi:10.1186/s12864-020-6540-1. Can the internode distance be corrected? This is entry F of `wide_scan.md` §2.3.

## 2. Prior art (details and DOIs in `prior_art.md`)

**Sample complexity.**

| Method | What is proven | Source |
|---|---|---|
| ASTRAL (exact ASTRAL*) | **Θ(f⁻² log n)**: matching upper bound m > 20 log(n/ε)/f² and worst-case lower bound | Shekhar, Roch & Mirarab, TCBB 2018, doi:10.1109/TCBB.2017.2757930 |
| NJst / ASTRID | **No proven sample-complexity bound.** Consistency only (Allman, Degnan & Rhodes, TCBB 2018, doi:10.1109/TCBB.2016.2604812). Roch 2018 (RECOMB-CG, doi:10.1007/978-3-030-00834-5_11) proves a *variance lower bound* for the internode distance: Var ≥ C·d/m, and ≥ C·n/m for the largest entry on some trees. He conjectures that NJ on all pairwise distances needs m ≥ linear in n, vs log n for ASTRAL, and lists the tight bound as open. | as given |
| Any method | m = Ω(1/f) | Dasarathy, Nowak & Roch, TCBB 2015, doi:10.1109/TCBB.2014.2361685 |
| GLASS (uses branch lengths) | O(1/f) | Mossel & Roch, TCBB 2010, doi:10.1109/TCBB.2008.66 |
| METAL | O(f⁻² log n) | Dasarathy, Nowak & Roch 2015 (as above) |
| STEAC | ∝ f⁻² | per Shekhar et al. 2018 |

**Empirical comparison.** The only direct comparison of genes needed is Shekhar et al. 2018.

- **Setup:** n = 8 only, target P(correct) ≈ 0.9, f from 0.005 to 0.1.
- **Scaling:** ASTRAL-II and NJst both scale as about 1/f².
- **ASTRAL vs NJst:**
  - Caterpillar and balanced trees: ASTRAL needs 10–20% fewer genes.
  - Double-quartet tree: NJst needs fewer.

**Gap.** Nothing published on n > 8, on per-branch k₉₅(f), or on the conjectured n-scaling.

**Missing data.**
- **Prior work:**
  - **Asteroid** (Morel, Williams & Stamatakis, *Bioinformatics* 2023, doi:10.1093/bioinformatics/btac832) is a missing-data-aware internode-distance method. Each gene is scored against the induced species tree. The authors claim consistency under any per-gene random deletion. **This is the main prior art for sub-idea B.**
  - Weighted ASTRID (Liu & Warnow, AMB 2023, doi:10.1186/s13015-023-00230-6) does not address missing data.
  - OCTAL (doi:10.1186/s13015-018-0124-5) completes gene trees.
- **Not found:** a 2020–2026 paper on an inverse-probability (Horvitz–Thompson) correction of ASTRID.

## 3. Data, tools and baseline reproduction

**Data and simulator.**
- **Nute et al. 2018 data** (doi:10.13012/B2IDB-7735354_V1):
  - Taxa and conditions: 25 ingroup taxa + outgroup; tree heights 10M / 2M / 500K generations (low / high / very high ILS); 2 speciation rates; 20 replicates.
  - Contents: true and RAxML gene trees; complete data plus i.i.d. missing data at 30% and 60%; published ASTRAL (4.10.5), ASTRID (v1), MP-EST and SVDquartets trees.
  - I estimated the effective population size from rooted-triplet concordance (`code/estimate_ne.py`): Ne = 2.02×10⁵ (2.00, 2.03 and 2.01×10⁵ per height). So branch length in CU = generations / 2×10⁵, as in the ASTRAL-II SimPhy setup.
- **Own MSC simulator** (`code/msc.py`): one lineage per species. Checked against theory: P(gene quartet = species quartet) = 1 − ⅔e^{−f} is reproduced to within Monte Carlo error at f = 0.05, 0.2 and 1.

**Tools.**
- **ASTRAL:** ASTER `astral` v1.25 (bioconda `aster`). I used `-r 1 -s 0` in the large sweeps, because ASTRAL time is superlinear in k (n = 64, k = 5000 took 387 s per run).
- **FastME:** 2.1.6 (bioconda).
- **Asteroid:** built from source.
- **ASTRID-2** (github pranjalv123/ASTRID) **failed to build** (bazel dependency `@bazel_skylib` not resolvable). I re-implemented ASTRID instead, which is mean internode distance + FastME BME with NNI+SPR (`code/stree.py`), and NJst (same distance + NJ).

**Baseline reproduction** (`code/baseline.py` → `results/baseline_*.csv`):
1. **Table S1 of Nute et al. reproduced** (average gene-tree/species-tree distance AD, gene-tree estimation error GTEE, total discord TD; 1000 genes; 20 replicates). All 18 means agree to 2 decimals, e.g. deep/500K: AD 0.75, GTEE 0.44, TD 0.81. Standard deviations agree to ±0.01.
2. **Published ASTRID distance matrices reproduced exactly** (max |Δ| = 2e-15). This works only if I count **edges in the rooted gene tree**, i.e. ASTRID v1 counts the root node of rooted inputs, which is a small deviation from the NJst definition.
3. **Rerunning ASTRID on the same genes gives the same species-tree topology as the published ASTRID tree:**
   - 96% of complete-data / true-gene-tree datasets.
   - 85% with RAxML gene trees.
   - 65–88% with missing data.

   The difference is in the tree search (FastME version and options), not the matrix. Mean FN of my rerun vs the published trees (1000 genes, full data) agrees within 0.002 in every condition.
4. **Rescoring the published species trees** (`results/baseline_published.csv`, 1000 genes, complete data, mean FN):

| condition | ASTRAL true | ASTRID true | MP-EST true | SVDq true | ASTRAL RAxML | ASTRID RAxML |
|---|---|---|---|---|---|---|
| deep, 10M | 0.004 | 0.004 | 0.007 | 0.098 | 0.063 | 0.057 |
| deep, 2M | 0.024 | 0.022 | 0.026 | 0.109 | 0.030 | 0.030 |
| deep, 500K | 0.061 | 0.054 | 0.089 | 0.202 | 0.100 | 0.111 |
| recent, 10M | 0.000 | 0.000 | 0.000 | 0.039 | 0.011 | 0.011 |
| recent, 2M | 0.013 | 0.015 | 0.020 | 0.072 | 0.028 | 0.022 |
| recent, 500K | 0.061 | 0.059 | 0.085 | 0.211 | 0.161 | 0.150 |

## 4. Sample-complexity curves

### 4a. Nute et al. data, true vs estimated gene trees, 3 ILS levels

**Script:** `code/nute_curves.py` and `code/analyze_nute.py`.
**Plots:** `results/nute_fn_vs_k.png` and `results/nute_k95.png`.
**Setup:** 6 conditions × 20 replicates; k = 10–1000 (first k genes); per-branch recovery recorded with the branch length in CU.

- **The three methods are tied at every ILS level and k.** For example, very high ILS, true gene trees, k = 1000: FN is 0.062 (ASTRAL), 0.059 (ASTRID) and 0.058 (NJst). With RAxML trees at the same setting: 0.132 / 0.128 / 0.126.
- **Per-branch recovery at k = 1000** (`results/nute_branch_bins.csv`):

| f bin (CU) | true: ASTRAL | true: ASTRID | RAxML: ASTRAL | RAxML: ASTRID |
|---|---|---|---|---|
| [0, 0.01) | 0.48 | 0.47 | 0.42 | 0.39 |
| [0.01, 0.03) | 0.79 | 0.82 | 0.62 | 0.65 |
| [0.03, 0.1) | 0.99 | 0.99 | 0.90 | 0.89 |
| [0.1, 0.3) | 1.00 | 1.00 | 0.93 | 0.94 |

- **Logistic fit** of P(recover) on (log k, log f) for branches with f < 1:
  - k₉₅(f = 0.1): about 220–250 genes with true gene trees, about 1300 with RAxML trees, for all three methods.
  - **Gene-tree estimation error costs about 5× more genes.**
  - The pooled fit gives a log–log slope of about −1.2 to −1.35, flatter than f⁻². But this pools branches in different tree contexts, and f < 0.03 is extrapolated beyond k = 1000. Treat the slope as descriptive only.

### 4b. Controlled MSC simulation: every internal branch = f

**Script:** `code/sc_sim.py` and `code/analyze_sc.py`.
**Plots:** `results/sc_curves_tree.png`, `results/sc_k95_vs_f.png` and `results/sc_k95_vs_n.png`.

**Setup.**
- **Trees:** caterpillar and balanced, n = 8, 16, 32, 64; f = 0.05, 0.1, 0.2, 0.5 CU.
- **Genes:** true gene trees, nested k grid 10 … 5000.
- **Replicates:** __SCREPS__ per cell.
- **k₉₅:** the smallest k at which P(species tree exactly correct) ≥ 0.95, log-interpolated on the grid.
- **ASTRAL cap:** k ≤ 2000 at n = 32 and k ≤ 1000 at n = 64. A ">" entry means the target was not reached within the run grid.

**k₉₅ (whole tree) — caterpillar:**

__SCTABLE_CAT__

**k₉₅ (whole tree) — balanced:**

__SCTABLE_BAL__

**Findings.**
- **Dependence on f.** Within each (shape, n), k₉₅ scales as about f^−1.4 to f^−2.2 (`results/sc_sim_slopes.csv`; vs 1−e^{−f}: −1.5 to −2.25). This is consistent with the f⁻² theory for ASTRAL, with a coarse grid. **There is no evidence that ASTRID's f-exponent differs from ASTRAL's.**
- **Dependence on n — the interesting part.** The constant c = k₉₅·f² (fit over f):

| | n = 8 | 16 | 32 | 64 |
|---|---|---|---|---|
| caterpillar, ASTRAL | __C_cat_astral__ |
| caterpillar, ASTRID | __C_cat_astrid__ |
| balanced, ASTRAL | __C_bal_astral__ |
| balanced, ASTRID | __C_bal_astrid__ |

- On caterpillars, ASTRID's constant grows clearly faster than ASTRAL's (about 4.9 → 15.3 vs 4.0 → 7.6 from n = 8 to 64). On balanced trees the two grow alike up to n = 64 (3.2 → 8.6 vs 3.4 → 7.6). Caterpillars have the longest leaf-to-leaf paths, so internode distances have the largest variance; Roch 2018's lower bound is driven by exactly this.

### 4c. Focused n-scaling test (f = 0.2, n up to 128)

**Script:** `code/nscale.py` and `code/analyze_nscale.py`.
**Plot:** `results/nscale.png`.
**Setup:** fine k grid; __NSREPS__ replicates; 90% bootstrap CIs over replicates in `results/nscale_k95.csv`.

__NSTABLE__

**Interpretation.**
- **ASTRAL vs ASTRID in n.** ASTRAL's k₉₅ grows slowly in n, consistent with log n. ASTRID's grows faster.
- **Conjectured scaling not yet tested.** Over n = 8–128 the fitted exponents in n are __NSEXP__. Neither method looks linear in n at f = 0.2. With 20 replicates per cell, k₉₅ near 0.95 is noisy: the bootstrap CIs overlap at several n. The clearest separations are balanced n = 128 (ASTRAL 150 [135–245] vs ASTRID 346 [271–373]) and caterpillar n = 32 (96 [93–98] vs 200 [186–363]).
- **What a project needs.** 100+ replicates per cell, n up to 256–512, and smaller f. Those runs are cheap for ASTRID (seconds). ASTRAL is the bottleneck.

## 5. New method: missing-data-corrected ASTRID

**Derivation** (`code/stree.py`, `gene_contrib`).

- **Setup.** Under i.i.d. deletion with rate p (q = 1−p), take an internal node v of the full gene tree on the i–j path. It stays a node of the restricted tree iff at least one of the K_v leaves on its off-path side survives. We observe k_v ~ Bin(K_v, q), and the node only if k_v ≥ 1.
- **Unique unbiased per-node weight.** The weight that depends only on k is

  **w(k) = 1 − (−p/q)^k**, with w(0) = 0,

  because Σ_k C(K,k) q^k p^{K−k}(1 − (−p/q)^k) = (q+p)^K − (p−p)^K = 1 for every K ≥ 1.
- **Consistency.** Summing w over the observed path nodes gives an unbiased estimate of the complete-data internode distance, conditional on i and j being present. Averaged over genes, it converges to the complete-data NJst matrix. That matrix is additive on the species tree, so ASTRID/NJst with this weight is **statistically consistent under i.i.d. deletion for every p < 1**. p is estimated as 1 − (mean leaves per gene)/n.
- **Catch: variance.** |p/q| ≥ 1 when p ≥ 0.5, so the weights alternate in sign and grow geometrically; at p = 0.6, w(10) = −57.
- **Variants tested.**
  - **astrid-ht:** exact weights.
  - **astrid-plugin:** bounded, biased, w = 1/(1 − p^{k/q}).
  - **ht1:** first-order, w(1) = 1/q, otherwise 1.
  - **astrid-cmp:** the completion alternative. Each gene tree is completed by greedy OCTAL-like insertion against a first-pass ASTRID tree rooted at the outgroup, then ASTRID is rerun.

**How big is the bias being corrected?** (`code/asym_bias.py` → `results/asym_bias.csv`)

- **Method.** The expected observed internode count under deletion is exactly Σ_v(1 − p^{K_v}). So the infinite-gene limit of plain ASTRID can be computed from full gene trees (4000 MSC genes, paired with the complete-data matrix from the same genes).
- **Equal-branch caterpillar and balanced trees** (n = 8–32, f = 0.05–1): the limit tree **never** changes, even at p = 0.9.
- **Random Yule trees** (n = 10–30, mean branch length 0.05–1 CU, 90 trees): the limit tree differs from the complete-data tree in 17% / 18% / 24% / 31% of trees at p = 0.3 / 0.5 / 0.7 / 0.9. It is worse in 10–19% and better in 4–7%. The mean FN cost is only +0.005 (p = 0.3) to +0.009 (p = 0.9), on a base of 0.054.
- **So this is the ceiling for any correction:** less than 0.01 FN, and only with very many genes.

**Experiment** (`code/md_exp.py` → `results/md_wilcoxon.csv`, `md_meanFN.csv`, `md_runtime.csv`).

- **Data:**
  - Nute et al. rand-30 and rand-60 (published i.i.d. deletion; k = 50, 200, 1000).
  - My own i.i.d. deletion at p = 0.5 and 0.7 on the complete replicates (k = 1000).
  - All 6 conditions; true and RAxML gene trees; 1920 datasets.
- **Split:** odd replicates = DEV (used to choose the variant), even replicates = held-out TEST.
- **Metric:** species-tree FN rate.
- **Statistics:**
  - Paired by dataset.
  - Tie band: |Δ| < 1e-9. FN moves in steps of 1/23, so a tie means equal FN.
  - Two-sided Wilcoxon signed-rank test, zeros dropped.

DEV screen (also `/opt/runs/dev_variants.jsonl`, rand-60 and sim-70):
- **Exact HT is unusable at p ≥ 0.5:** +0.56 to +0.80 FN.
- **At p = 0.3 exact HT is slightly better than ASTRID:** −0.003 FN, not significant (p = 0.13–0.65).
- **ht1:** +0.0045 FN, W/T/L 110/465/145.
- **plugin:** +0.0002 FN, W/T/L 55/606/59.
- **Completion:** slightly worse.

The DEV screen picks **plugin** as the variant to carry forward.

**Held-out TEST, all missing-data settings pooled** (480 paired datasets per row; Δ = method − comparison, negative = better):

| gene trees | method | vs | mean FN | Δ FN | W/T/L | p |
|---|---|---|---|---|---|---|
| true | astrid-plugin | ASTRID | 0.1253 | +0.0001 | 29/418/33 | 0.71 |
| true | astrid-ht | ASTRID | 0.4567 | +0.3315 | 25/193/262 | <1e-10 |
| true | astrid-cmp | ASTRID | 0.1312 | +0.0060 | 8/416/56 | <1e-6 |
| true | Asteroid | ASTRID | 0.1207 | −0.0044 | 135/237/108 | 0.13 |
| true | ASTRAL | ASTRID | 0.1264 | +0.0013 | 125/231/124 | 0.53 |
| true | astrid-plugin | ASTRAL | 0.1253 | −0.0012 | 126/229/125 | 0.71 |
| RAxML | astrid-plugin | ASTRID | 0.1940 | −0.0005 | 48/400/32 | 0.57 |
| RAxML | astrid-ht | ASTRID | 0.4947 | +0.3002 | 46/165/269 | <1e-10 |
| RAxML | astrid-cmp | ASTRID | 0.1999 | +0.0053 | 16/398/66 | <1e-6 |
| RAxML | Asteroid | ASTRID | 0.1865 | −0.0081 | 162/197/121 | 0.0035 |
| RAxML | ASTRAL | ASTRID | 0.1955 | +0.0009 | 151/172/157 | 0.86 |
| RAxML | astrid-plugin | Asteroid | 0.1940 | +0.0075 | 119/205/156 | 0.010 |

**Per setting (TEST):**
- **Asteroid is the only consistent winner, and only at high deletion.**
  - rand-60, RAxML: −0.0135, p = 0.04.
  - sim-70, RAxML: −0.0174, p = 0.005.
- **astrid-plugin is never significantly different from ASTRID** (all p > 0.13).
- **My ASTRID vs the published ASTRID v1:** mine is −0.005 FN (p = 0.045 true, 0.019 RAxML). The difference is the unrooted vs rooted distance and the FastME search, and it is not a method contribution.

**Runtime** (seconds per species tree, 26 taxa, mean over all runs, 1 thread, machine shared with other jobs):

| k | ASTRID (mine, Python) | astrid-ht / plugin | astrid-cmp | Asteroid (C++) | ASTRAL (ASTER) |
|---|---|---|---|---|---|
| 50 | 0.08 | 0.08 | 0.26 | 0.04 | 0.23 |
| 200 | 0.27 | 0.26 | 0.95 | 0.11 | 0.82 |
| 1000 | 1.08 | 1.08 | 4.46 | 0.46 | 4.26 |

The correction itself costs nothing extra. My Python distance code dominates the time, and a C++ ASTRID would be about 10× faster.

## 6. Verdict, a 4-week plan, and risks

**Verdict: overall promising, but only for part A (sample-complexity study). Part B is not promising.**

**4-week plan for A**, "Does ASTRID need more genes than ASTRAL as n grows? An empirical test of Roch's conjecture":
1. **Week 1.**
   - Port the internode-distance accumulation to C++ or numba.
   - Switch to ASTRAL's default search and use an exact-ASTRAL check for n ≤ 16.
   - Add Asteroid and wASTRID as extra internode methods.
   - Fix the protocol: per-branch and whole-tree k₉₅ with 200 replicates per cell, using bisection on k instead of a grid, as in Shekhar et al.
2. **Week 2.**
   - **Trees:** caterpillar, balanced, and random Yule trees.
   - **Grid:** n = 8 … 512, f = 0.02 … 0.5, all with true gene trees.
   - **Fits:** k₉₅ vs n against the candidate forms a·log n, a·n^b and a·n; k₉₅ vs f against the f⁻² prediction.
   - **Key output:** a plot of k₉₅(n) for ASTRAL vs ASTRID on caterpillars.
3. **Week 3.**
   - **Estimated gene trees:** re-run with them (seq-gen + FastTree at 2–3 sequence lengths), and on Nute et al. / ASTRAL-II data.
   - **Mechanism:** measure the variance of the internode distance for the pairs that drive errors, and check whether "short-path" distance methods remove the n-penalty, as Roch suggests.
4. **Week 4.** Write-up.
   - **Optional theory:** a lower-bound construction, a caterpillar where NJ on internode distances needs Ω(n·f⁻²) genes.
   - **Missing-data section:** the unbiased weight and why it fails (variance), as a short negative result.

**Risks.**
- **(i) The n-effect may be weaker than it looks.** The current evidence uses 20–30 replicates, and many k₉₅ CIs overlap. It may turn out to be a constant factor, not a different rate. That is still a publishable answer for a class project, but less striking.
- **(ii) Overlap with Shekhar et al. 2018 at n = 8.** The novelty is n > 8 and the n-scaling.
- **(iii) ASTRAL cost limits large-n, large-k runs.** ASTRAL at n = 64 and k = 5000 took 387 s per run on this 4-core machine; the whole pilot's simulation used about 3.5 h on 4 shared cores. Mitigation: bisection on k and fewer, targeted cells.
- **(iv) This is an evaluation study**, not a new method. The course may favor method projects.
- **(v) Search heuristics: checked, minor.** I compared ASTER's reduced search (`-r 1 -s 0`, used in the sweeps) with the default on 60 paired runs (`code/astral_search_check.py` → `results/astral_search_check.json`; n = 32–64, f = 0.1–0.2). It returned the identical tree in 59/60 runs. In the one other run it was slightly worse, so ASTRAL's k₉₅ here is, if anything, very slightly pessimistic.
