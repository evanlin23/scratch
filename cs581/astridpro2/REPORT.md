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
   - vs **ASTRAL-Pro3**: __APRO__. Our ASTRAL-Pro3 matches the authors' published ASTRAL-Pro trees on the same
     inputs (0.0855 vs 0.0863, p = 0.86), so the baseline is sound.
   - vs **ASTRID-multi**: __MULTI__. The gain is concentrated in the high-duplication DISCO conditions: −0.016,
     14/3/1, p = 0.003.
   - **Ties** with ASTRID-DISCO (__DISCO__) and Asteroid (__AST__).
   - Beats DISCO+ASTRAL, FastMulRFS and DupLoss-2. Ties wQFM-GDL (the 2026 state of the art) on the FastMulRFS data.
   - Empirical: 1KP C12 has 5 FN vs the 1KP ASTRAL reference, the same as ASTRAL-Pro3, against 10 for ASTRID-multi.
     Fungi16 and vertebrates188: all methods agree.
3. **Speed.**
   - 0.2 s (100 genes) and 1–40 s (1000 genes, up to 3,700-leaf families) on 100 taxa, against 6–29 min for
     ASTRAL-Pro3 on DISCO (1000 genes, 100 taxa, 1 thread, shared machine).
   - ASTRAL-Pro3 timed out (> 40 min) on the 3 heaviest DISCO conditions.
   - At 1000 taxa: 14 s, where ASTRAL-Pro3 and FastMulRFS exceed 20 min at 500 taxa.
   - 10,000 genes: 8 s.
   - About **100–500× faster than ASTRAL-Pro3**, with memory independent of the number of genes (streaming).

**The honest framing.**
- ASTRID-Pro is a **provably consistent, very fast distance method that matches or beats ASTRAL-Pro3** in these
  benchmarks.
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
| FastMulRFS data (= ASTRAL-Pro "S100") | Molloy & Warnow 2020, doi:10.13012/B2IDB-5721322_V1 | 100 species. DL rate 1e-10, 2e-10, 5e-10 × Ne 1e7, 5e7. RAxML gene trees from 25 and 100 bp. Reps 01–10 at 100 genes; reps 01–05 at 500 genes. |
| DISCO data | Willson et al. 2022, doi:10.13012/B2IDB-4050038_V1 | 100 species, 1000 estimated gene trees (100 bp). 12 conditions: default; gdl_{1e-10, 5e-10, 1e-9}_{0, 0.5, 1} (families up to 3,700 leaves); ils_1e4; ils_2e8; missing_1000. Fast methods on reps 01–05; ASTRAL-Pro3 on reps 01–02. |
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

__TABLES__

## 6. Caveats

- **Runtimes were measured on a shared machine.** Up to 4–5 single-threaded jobs ran at once on 4 cores (wQFM-GDL
  is multi-threaded Java). Absolute seconds are inflated, perhaps up to 2×, for every method alike. The ratios
  (100–500×) are far larger than this noise.
- **ASTRAL-Pro3 comparisons on DISCO cover reps 01–02 only** (14 finished pairs). Three high-duplication conditions
  are missing because ASTRAL-Pro3 timed out at 40 min. The pooled ASTRAL-Pro3 comparison is therefore driven by the
  FastMulRFS data (57 pairs).
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

