# Pilot: is ASTRAL-Pro consistent under rooting/tagging error? Which methods are consistent under GDL / DLCOAL?

*CS581 project pilot, 2026-10-10, about 2.5 h wall-clock on 4 cores. Branch `claude/cs581-gdlcons`. Code is in `code/` (see `code/README.md`); results are in `results/`. The theory note is `results/theory.md`; the prior-art table with DOIs is `results/prior_art.md`.*

## 0. Verdict: **promising** (for question 1); unclear for question 2

**Question 1 (ASTRAL-Pro under rooting/tagging error): in general, no.** I found three kinds of 4-taxon counterexamples under pure GDL, all with *true* gene trees and branch-specific rates. In each one, ASTRAL-Pro with the true root and true tags converges to the right tree:

| failure | rooting | tagging | example | ASTRAL-Pro wrong in |
|---|---|---|---|---|
| (a) hidden paralogy | correct | species-overlap rule (ASTRAL-Pro's own, given the root) | cand2 | 4/4 datasets of 10,000 families; stock ASTRAL-Pro3 also 4/4, and 4/4 in an 8-taxon embedding |
| (b) **random rooting error** | re-root at a uniformly random edge w.p. p = 0.3 | overlap rule | cand3 (correct root is consistent here) | 4/4 at 10,000 |
| (c) **ASTRAL-Pro's own min-duplication rooting** | ASTRAL-Pro3 heuristic | ASTRAL-Pro3 | cand4 (moderate rates: about 20 copies; correct-root overlap tags are consistent) | 20/20 at 1,000 and 4/4 at 5,000; 4/4 in an 8-taxon embedding |

In (c), the stock software puts the root inside the high-turnover C lineage in 43% of families and at the true position in only 15%.

**Exact theory.** For (a) and for random hidden-paralog mislabelling `d2sv(q)` I derived the limiting ASTRAL-Pro scores in closed form, as an integral over duplication times with exact birth–death generating functions (`results/theory.md` Section 3).
- Consistency holds iff O + H_AB > max(H_AC, H_BC), where O counts orthologous classes and the H terms count hidden-paralog classes.
- Under random mislabelling there is a sharp threshold q* = O / (H_wrong − H_AB).
- The formula matches ASTRAL-Pro3's own scores within about 1.5 SE:
  - correct − AC|BD = −0.060 ± 0.003 observed vs −0.0555 predicted (cand2, overlap tags);
  - predicted q* = 0.079, observed sign change between 0.10 and 0.15; the margins at all six q values tested lie within 1.5 SE of the prediction.

**How common, and caveats.**
- Under the formula, only **12 of 15,434** random rate configurations are inconsistent with correct-root overlap tags. All have very high duplication or turnover on one branch, and one of them is *critical* (λ = μ), so "no supercritical branch" is not sufficient.
- On all standard-looking settings (16- and 50-taxon Yule trees, λ = μ, the sibling's adversarial rates), every error model converged to zero error by 1,000 families (most by 200); errors cost samples, not the limit. On generic small-tree grids, their margins stayed ≥ 0.51.
- So consistency fails in **corners** of parameter space, not in typical simulations.
- "Naive" tag flips that put S on nodes whose children share species break ASTRAL-Pro at q ≈ 0.035 even where overlap tags are fine (cand1). ASTRAL-Pro's tagger can never output such labels, so the slide question needs an error model that respects its tagging rule (Def. 1).

**Question 2 (which other methods are consistent).**
- In the generic atlas (true gene trees, 4 GDL and 4 DLCOAL settings with SimPhy, up to 5,000 families), every method converges, or nearly so. No clear plateau appears. STAG is slowest; ASTRID-multi has small residual errors under the highest-duplication DLCOAL setting.
- On the adversarial 4-taxon pools, methods split by mechanism, and **no method is robust on all of them** (Section 5.2):
  - ASTRAL-Pro, ASTRAL-DISCO and ASTRID-DISCO share ASTRAL-Pro's rooting and tagging. They fail on the hidden-paralog pool (cand2) and the own-rooting pool (cand4): 100% wrong.
  - STAG is right on cand2 but wrong on cand4.
  - FastMulRFS and ASTRID-multi are right on cand4 but wrong on cand1 (4/4 at 10,000 families), a pool where ASTRAL-Pro is fine.
  - ASTRAL-multi, the only multi-copy method with a GDL proof, sits at near ties on cand1 and cand2.
- The literature has proofs only for ASTRAL-one/multi (GDL, DLCOAL), ASTRAL-Pro and ASTRAL-DISCO (GDL, *correct tags*), and FastMulRFS (non-adversarial GDL). DLCOAL proofs for tag-based methods exist only as recent, possibly shaky preprints or papers (Section 2).

**Why "promising".**
- The project has a clean, new, checkable core: an exact formula, explicit counterexamples in the shipped software (both from tagging and from rooting), and a threshold law.
- There is a natural theorem to aim for: consistency of ASTRAL-Pro with overlap tags when duplication and turnover are bounded. There is also a natural method question: robust rooting.
- It reuses the sibling session's simulator and fits in 4 weeks (Section 7).

**The main risk.** The counterexamples use unusual branch-specific rates (cand2 is extreme, cand4 moderate), so a reviewer may call them pathological. The closest prior work is Parsons et al. 2026 (DLCOAL tagging correctness). I could not read its theorems, and it must be checked before committing to the project.

## 1. Questions (verbatim)

Source: lecture deck *Phylogenomics part 2*, slide 35, "(Some) Open Questions" (`../notes/lectures/CS581-phylogenomics-2026-part2.txt`, lines 208–216; `../literature/slide_open_problems.md`):

1. "Is ASTRAL-Pro statistically consistent for GDL under some random model of error for rooting and tagging?"
2. "Which other existing species tree estimation methods are statistically consistent under GDL and DLCOAL models?"

Context: slide 34 says quartet-based estimation is consistent under ILS (Allman et al.), GDL (Legried et al.) and DLCOAL (Markin & Eulenstein).

## 2. Prior art

Full table with 26 rows and verified DOIs: `results/prior_art.md`. [A] means read in the paper text or abstract; [M] means memory or inference.

| method | GDL | DLCOAL | rooting/tagging assumption | source |
|---|---|---|---|---|
| ASTRAL-one / ASTRAL-multi | **proved** | **proved** | none (untagged MUL-trees) | Legried, Molloy, Warnow, Roch, JCB 2021, doi:10.1089/cmb.2020.0424; Markin & Eulenstein, Bioinformatics 2021, doi:10.1093/bioinformatics/btab414. Sample complexity: Hill, Legried, Roch, AAP 2022, doi:10.1214/22-aap1799 |
| ASTRAL-Pro (MLQST) | **proved, given correct root + tags** (Thm 2); errors only "suspected" fine [A] | open; a "major open question" [A: talk abstract] | correct tags | Zhang, Scornavacca, Molloy, Mirarab, MBE 2020, doi:10.1093/molbev/msaa139 |
| ASTRAL-Pro 2 / 3 | no new theory | — | — | doi:10.1093/bioinformatics/btac620; ASTER README |
| ASTRAL-Pro, exclusion-only objective | — | claimed in a preprint (v1). In v2 the consistency wording is gone from title and abstract [A: abstracts only] | new definition of "correct duplication" under ILS | Parsons, Liu, Dua, Markin, Molloy, bioRxiv 2026, doi:10.64898/2026.01.20.700722 |
| ASTRAL-DISCO | **proved, given correct root + tags** | claimed (wQFM-DISCO paper); "correct tagging" under ILS is undefined there [A] | correct tags | Willson et al., Syst Biol 2022, doi:10.1093/sysbio/syab070; Hakim, Ratul, Bayzid, Bioinf Adv 2024, doi:10.1093/bioadv/vbae189 |
| ASTRID-DISCO, ASTRID-multi, NJst | unproven. NJst/ASTRID are *inconsistent* under random missing data | unproven | — | Rhodes, Nute, Warnow, arXiv 2001.07844; sibling session `claude/cs581-gdl` |
| FastMulRFS | **proved if no "adversarial GDL"** (Thm 6) | — | none | Molloy & Warnow, Bioinformatics 2020, doi:10.1093/bioinformatics/btaa444 |
| STAG, SpeciesRax, DupLoss-2 | empirical only | empirical only | — | doi:10.1101/267914; doi:10.1093/molbev/msab365; doi:10.1093/sysbio/syaf073 |
| DupTree / gene-tree parsimony | **open** (listed as open by Molloy & Warnow 2020 and Willson et al. 2021) | GTP is *inconsistent* under the MSC alone | rooted | Wehe et al. 2008, doi:10.1093/bioinformatics/btn230; Sapoval & Nakhleh, RECOMB-CG 2026, doi:10.1007/978-3-032-26891-4_9 |

**Not found.**
- Any paper on ASTRAL-Pro or DISCO under *rooting/tagging error*.
- Any result showing that the species-overlap (LCA) tagging rule can make ASTRAL-Pro inconsistent.

Caveat: I could not read Parsons et al.'s theorems (bioRxiv rate-limited); their "accuracy of A-Pro tagging" section may overlap with Section 5.2 below.

## 3. Methods and tools

**Simulators.**
- *GDL.* The sibling session's pure-GDL simulator (`gdlsim.py`): branch-specific (λ, μ), true root, true D/S labels, families with at least 4 species kept.
- *DLCOAL.* SimPhy 1.0.2: 16 or 50 taxa, height 2·10⁶ generations, Ne 2·10⁵ or 10⁶, dup/loss 3·10⁻⁷ to 1.5·10⁻⁶ per generation, 20,000 families per replicate.

**Error models** (`core.apply_error`; definitions in `results/theory.md` Section 2):
- `true`: true root and tags.
- `ovl`: true root, species-overlap tags.
- `rovl(p)`: re-root at a random edge with probability p, then overlap tags.
- `rtrue(p)`.
- `rsp-X`: root on a leaf of species X.
- `d2sv(q)`, `flipv(q)`: errors that respect ASTRAL-Pro's Def. 1.
- `d2s(q)`, `flip(q)`, `s2d(q)`: naive flips.
- `own`: ASTRAL-Pro3's own rooting/tagging.

**Feeding root and tags to ASTRAL-Pro.**
- The CLI cannot accept them, so I patched ASTER (29 lines; `code/astral-pro-fixed.patch`).
- With `APRO_FIXED=1` the binary keeps the input root and reads the tags from branch lengths. Otherwise it is the stock ASTRAL-Pro3 (v1.25.3.8 source).
- Quartet supports come from `-C -u 3` (`freqQuad.csv`) on the true species tree.
- Pitfall: the `t2`/`t3` labels in that file are not stable across runs. Always parse the split string; one early analysis had to be redone because of this.

**Atlas methods.**

| method | implementation |
|---|---|
| ASTRAL-Pro3 | bioconda `aster` |
| ASTRID-multi | the sibling's re-implementation + FastME |
| ASTRID-DISCO | DISCO v1.4.1 (`disco.py`) + ASTRID re-implementation |
| ASTRAL-DISCO | DISCO + ASTER astral4 |
| FastMulRFS | built from source: preprocess v3 + FastRFS `.single` |
| STAG | original Python 2 code; uses only families containing all species |
| DupTree | **not installable**: GitHub needs authentication, no conda package |
| ASTRAL-multi | ASTER astral4 with a gene→species map; used as a control on the 4-taxon pools |

**Design.** The datasets are nested: the first K families of each 20,000-family file form the K-family dataset.
- K ∈ {50, 200, 1,000, 5,000}, plus 20,000 for the DLCOAL settings.
- 3 replicate species trees per setting.
- Metric: FN rate = missed true bipartitions / (n − 3).

## 4. Results: ASTRAL-Pro under rooting/tagging error

### 4.1 Generic settings: no plateau for any error model

**Curves** (`results/curves_summary.md`, `results/apro_error_curves.png`).
- **Settings.** 16 taxa: critical λ = μ = 1; high 3/3; the sibling's adversarial (3, 0.5)/(0, 3); and a root duplication burst. 50 taxa: lossy, and adversarial.
- **Error models.** true, ovl, rovl 0.1/0.3/1, flip 0.05/0.15/0.3, d2s 0.5/1, s2d 0.5.
- **Result.** Every model reaches FN = 0 by 1,000 families, most by 200. Naive all-S tagging (d2s q = 1) is slowest: it needs 1,000 families on gdl-adv and 5,000 on gdl50-lossy.
- **Errors cost samples, not consistency, in these settings.** On gdl50-lossy at 50 families (2 replicates), the FN rate is 0.032 with true tags, 0.043 with ovl, 0.074 with rovl(0.3), 0.149 with d2s(0.5) and 0.160 with flip(0.3).

**Pure GDL with true gene trees has no gene-tree noise.** So the informative quantity is the limiting quartet *margin*, (correct − best wrong) / total ASTRAL-Pro support.

| scan | configurations | minimum margin |
|---|---|---|
| 4–5-taxon grid (`qgrid`) | 54 | ovl 0.72, rovl(1) 0.57, own 0.55 |
| systematic misrooting (`rootbias`) | 18 | 0.52 |
| random search (`search4`) | 13 informative | no negative margin with real signal |

Naive flips shrink the margin steadily as q grows: d2s(q = 1) has median 0.48 and minimum 0.017. They never go negative in this scan.

### 4.2 Exact theory and a counterexample for ASTRAL-Pro's own tagging

Details and numbers: `results/theory.md` Sections 3–4.

**Formula.** For (((A,B)x,C)y,D) with duplications above the ABC clade, the limiting per-family scores are:
- O + H_AB for the true topology;
- H_AC for AC|BD;
- H_BC for AD|BC.

O counts orthologous classes, and H_τ counts hidden-paralog duplication nodes of each pattern. Both are computed exactly (`code/predict_margin.py`).

**Counterexample cand2.** y-branch (λ, μ, t) = (8, 0.5, 0.5); x with no events; A (1, 4, 1); B (4, 8, 1); C (0.5, 1, 4); D a single copy.

Block standard errors, 40,000 families:

| tags | correct − AC\|BD, per family | predicted |
|---|---|---|
| ovl | **−0.060 ± 0.003** | −0.0555 |
| own (stock ASTRAL-Pro3) | **−0.136 ± 0.004** | — |
| d2sv(0.5) | −0.024 ± 0.002 | −0.0254 |
| true | +0.0049 ± 0.0003 | +0.0048 |

Species-tree error (fraction of disjoint datasets wrong; `results/curve4_cand2.jsonl`):

| tags | 100 families | 500 | 2,000 | 10,000 |
|---|---|---|---|---|
| true (sparse signal: about 0.005 orthologous classes per family) | 24/40 | 4/40 | 0/20 | 0/4 |
| ovl | 35/40 | 39/40 | 20/20 | 4/4 |
| **own** | 38/40 | **40/40** | **20/20** | **4/4** |

- **Embedding.** In an 8-taxon tree containing cand2, stock ASTRAL-Pro3 is wrong in 4/4 datasets of 5,000 families; true tags are right in 4/4.
- **Threshold.** For `d2sv(q)` the predicted q* is 0.079; the observed sign change lies between q = 0.10 and 0.15 (table in theory.md).
- **Prevalence.** 12/15,434 random configurations are inconsistent with overlap tags (11 have a supercritical y-branch; the 12th is critical, λ = μ = 8 over length 2). A finite q* exists in 353 configurations, and in only 12 of these is q* < 1.

### 4.2b Rooting error: random re-rooting (cand3) and ASTRAL-Pro's own rooting (cand4)

Details: `results/theory.md` Section 4b.

**cand3.**
- Rates: y (8, 0, 0.2); x (2, 0.5, 0.05); A (0.5, 4, 1); B (8, 8, 4); C (0, 8, 0.2).
- With the correct root, overlap tags give correct − AC|BD = +0.010 ± 0.003, and own gives +0.016 ± 0.004. Both are consistent: 0/4 datasets wrong at 10,000 families.
- Re-rooting each family at a uniformly random edge with probability **p = 0.3** gives −0.008 ± 0.003, and the species tree is wrong in **4/4** datasets at 10,000 families.
  - With p = 1: −0.055 ± 0.004, wrong 20/20 at 2,000.
  - Rooting on an A leaf (`rsp-A`): −0.086 ± 0.004.

**cand4.**
- Rates: y (4, 1, 1); x (1, 0.5, 0.5); A, B (0, 8, 0.5); C (8, 8, 0.2).

| rooting/tags | 1,000 families | 5,000 families |
|---|---|---|
| true | 0/20 | 0/4 |
| correct root + overlap tags | 3/20 | 0/4 |
| rovl(0.3) | 5/20 | 0/4 |
| **stock ASTRAL-Pro3** | **20/20** | **4/4** |

- Where ASTRAL-Pro3 roots 2,000 families: 868 inside C, 487 on the A or B lineage, and only 304 at the true root split (ABC | D).
- So min-duplication rooting is pulled toward the high-turnover lineage.
- In an 8-taxon embedding of cand4, stock ASTRAL-Pro3 is wrong in 4/4 datasets of 5,000 families; true tags and correct-root overlap tags are right in 4/4.
- **How common?** `code/rootscan.py` drew 56 random rate configurations with feasible copy numbers, all with consistent correct-root overlap tags, at 3,000 families each.
  - ASTRAL-Pro3's own rooting lowered the margin in 41/56 (median 0.97 vs 1.00) but never made it negative.
  - The minimum was +0.003, again with a high-turnover C branch (λ = μ = 8).
  - Side check: the formula's predicted overlap-tag margin matched simulation, with median |difference| < 0.001 and maximum 0.029.

### 4.3 Naive random flips (outside Def. 1)

**Configuration cand1:** y (4, 0, 1), x (0.5, 2, 0.01), A (0.25, 4), B (0.25, 1), C (0.25, 8). Here overlap tags, own and `d2sv` are all consistent (own: smallest margin +0.07 per family).

| naive model | correct − AC\|BD | species tree at 10,000 families |
|---|---|---|
| d2s(0.02) | +0.024 ± 0.002 | right (wrong 0/4) |
| d2s(0.035) | +0.001 ± 0.003 | — |
| d2s(0.05) | −0.022 ± 0.003 | — |
| d2s(0.1) | −0.109 ± 0.006 | wrong 4/4 |
| flip(0.1) | — | wrong 4/4 |

Why: an S label on a node whose children share species makes ASTRAL-Pro count copy tuples with multiplicity. The total support grows 35× at q = 0.1, and a small bias in that mass wins.

### 4.4 How good is ASTRAL-Pro3's own tagging on true trees?

Source: `results/tagging_accuracy.jsonl`, 1,000 families per setting.

| setting | true duplication fraction | own tag errors (D→S / S→D) | overlap rule on true root: error rate |
|---|---|---|---|
| gdl-crit | 0.24 | 1.3% (122 / 129) | 0.6% |
| gdl-high | 0.47 | 2.8% (702 / 174) | 2.2% |
| gdl-adv | 0.31 | 1.3% (1,111 / 462) | 0.9% |
| gdl-root | 0.13 | 2.3% (743 / 0) | 2.0% |

- Most errors are hidden paralogs (D→S), as Observation 1 predicts.
- The exact root position matches the truth in only 12–21% of families. That is expected: Claim 1 of Zhang et al. says any root along a speciation path gives the same score. So this metric overstates rooting error.

## 5. Results: method atlas (question 2)

### 5.1 Generic GDL and DLCOAL settings

Tables: `results/curves_summary.md`. Figure: `results/atlas_curves.png`.

Each cell below is the mean FN rate over 3 replicate species trees. Selected columns; everything is in `results/curves_summary.md`.

| setting | families | ASTRAL-Pro3 | ASTRID-multi | ASTRID-DISCO | ASTRAL-DISCO | FastMulRFS | STAG |
|---|---|---|---|---|---|---|---|
| gdl-crit / -high / -adv / -root (16 taxa) | 200 | 0 | 0 | 0 | 0 | 0 | 0 |
| gdl50-lossy (50 taxa) | 1,000 | 0 | 0.021 | 0 | 0 | 0 | (fails: no family has all species) |
| dlc-mod | 20,000 | 0 | 0 | 0 | 0 | 0 | 0 |
| dlc-ils (Ne = 10⁶) | 5,000 / 20,000 | 0.026 / 0.026 | 0.026 / 0 | 0 / 0 | 0.026 / 0.026 | 0.026 / 0 | 0.077 / 0 |
| dlc-high (dup = loss = 1.5·10⁻⁶) | 5,000 / 20,000* | 0 / 0 | 0.051 / 0.077 | 0.026 / 0 | 0 / 0 | 0 / 0 | 0.128 / 0.077 |
| dlc-asym (dup 1.5·10⁻⁶ > loss 3·10⁻⁷) | 20,000 | 0 | 0 | 0 | 0 | 0 | 0 |
| dlc50-ils (50 taxa, 1 rep) | 5,000 | 0 | 0 | 0 | 0 | — | — |

\*Only one replicate at 20,000: SimPhy produced fewer than 20,000 families with ≥ 4 species in replicates 0 and 2 (18,909 and 16,785).

- **Pure GDL with true gene trees is easy for every method.** All 16-taxon settings reach FN = 0 by 200 families, including the sibling's adversarial rates, so there is no plateau. At 50 families the only errors are FastMulRFS on gdl-high (0.31).
- **DLCOAL.** Every remaining error at 20,000 families is a single bipartition on one 3,407-generation internal branch of the replicate-1 species tree (about 0.002 coalescent units in dlc-ils). That is not evidence of inconsistency:
  - methods take turns failing it (ASTRAL-type in dlc-ils, ASTRID-multi and STAG in dlc-high);
  - ASTRAL-Pro's quartet margin there at 20,000 families is −1.4%, about the size of the sampling noise.
- **STAG** is the slowest. It cannot run at all when no family contains every species, which happened in the 50-taxon and high-duplication settings.
- **Bottom line for question 2 on generic data:** no method shows a plateau. The adversarial pools (5.2) are where they differ.


### 5.2 Adversarial 4-taxon pools

Fraction of disjoint datasets on which each method returns a wrong 4-taxon tree. True gene trees in every case; ASTRAL-multi is the ASTER control.

Full tables: `results/counterexample_summary.md`. Figure: `results/counterexample_curves.png` (top row: ASTRAL-Pro variants; bottom row: atlas).

| method | cand2: hidden paralogs (10,000 families) | cand4: own rooting fails (5,000) | cand1: ASTRAL-Pro fine (10,000) |
|---|---|---|---|
| ASTRAL-Pro3 | **4/4** | **4/4** | 0/4 |
| ASTRAL-DISCO | **4/4** | **4/4** | 0/4 |
| ASTRID-DISCO | **4/4** | **4/4** | 0/4 |
| ASTRAL-multi (astral4) | 2/4 (block margin −0.003 ± 0.002: tie) | 0/4 (margin +0.021 ± 0.003) | tie (±0.002) |
| ASTRID-multi | 3/4 (14/20 at 2,000) | 0/4 | **4/4** (12/20 at 2,000) |
| FastMulRFS | 3/4 (12/20 at 2,000) | 0/4 | **4/4** (20/20 at 2,000) |
| STAG | **0/4** | **4/4** | 0/4 |

- **No method is robust on all three pools.**
- The DISCO pipelines inherit ASTRAL-Pro's rooting and tagging, so they fail exactly when ASTRAL-Pro does.
- STAG survives hidden paralogs (closest-copy distances) but fails when a high-turnover lineage dominates (cand4).
- FastMulRFS fails on cand1. That is consistent with its theorem: duplications above the ABC clade followed by loss create "adversarial" bipartitions there.
- ASTRID-multi is wrong in 4/4 datasets of 10,000 families on cand1. But its limiting four-point sums are a near tie: AC|BD is shorter than AB|CD by 0.024 ± 0.016, out of sums of about 19.7, which is 1.5 SE (`code/astrid4.py`, `results/astrid4.log`). So it is a *candidate* for the sibling session's open question (is ASTRID-multi consistent under GDL?), not a demonstration.
- ASTRAL-multi sits at near ties on cand1 and cand2, so I make no claim either way.
- All of these use 4-taxon trees with extreme branch-specific rates. They show *how* each method can fail, not that the failures are typical.


## 6. Theory derived

`results/theory.md`:
1. **The error-as-selection view.** Rooting and tagging only select which quartets count. Their unrooted topologies are fixed.
2. **The one-sidedness of overlap-tag errors** on true roots under GDL (D→S on hidden paralogs only).
3. **The exact 4-taxon limiting-score formula** with its consistency condition O + H_AB > max(H_AC, H_BC), and the threshold q* = O / (max H_wrong − H_AB) for random hidden-paralog mislabelling.
4. **Rooting counterexamples** (Section 4b of theory.md). Random re-rooting with p = 0.3 (cand3) and ASTRAL-Pro's own min-duplication rooting (cand4) each break consistency, even where correct-root overlap tagging is consistent. No formula yet.
5. **A refutation of the easy proof route.** The per-node association inequalities (I1, I2) that would give a general proof fail on 4.6% (I1) and 7% (I2) of random survivable rate configurations.

All of these are checked numerically; none are written up as formal proofs yet.

## 7. What a 4-week project would look like

| week | work |
|---|---|
| 1 | Write the 4-taxon proposition rigorously: class decomposition, the integral formula, conditioning on the species filter, and the version with a root branch. Read Parsons et al. 2026 in full; get the PDF from a library or the authors. |
| 2 | **Positive result.** Try to prove ASTRAL-Pro with overlap tags is consistent under GDL when duplication is bounded, for example λ_e < μ_e with bounded turnover λ_e·t_e. "No supercritical branch" alone is *not* enough: the one non-supercritical failure in the scan is critical with huge turnover (λ = μ = 8 on a branch of length 2; formula only, not simulated). Use the formula to map the inconsistency region in (λ_y·T, loss asymmetry). Extend the formula to random rooting (`rovl(p)`). |
| 3 | **n-taxon and empirical side.** Show inconsistency for n > 4 (both 8-taxon embeddings already work empirically). Characterize when min-duplication rooting is attracted to high-turnover lineages (cand4). Rerun the atlas on the adversarial region, including estimated gene trees from sequences, to see whether the effect survives gene-tree error. Run under DLCOAL with SimPhy, where tag errors go both ways. |
| 4 | Write-up. Optional: a "tag-robust" ASTRAL-Pro variant that downweights S-tagged nodes whose children have very unbalanced copy numbers, or uses closest-copy information as STAG does. Test whether it fixes cand2 without hurting generic accuracy. |

**Deliverables either way:** an exact formula, explicit counterexamples in shipped software (tagging and rooting), a threshold law, and an atlas.

## 8. Risks

- **Pathological regime.** The tagging counterexample (cand2) needs λ ≈ 8 on one branch and heavy loss; only 0.08% of random configurations fail with correct-root overlap tags. The rooting counterexample (cand4) uses moderate growth but a high-turnover tip branch (λ = μ = 8 over 0.2). The honest framing: correct rooting and tagging are *necessary* for the theorem, and overlap tagging is safe in all but extreme regimes. Proving the "safe" part is the harder half. In a random scan, min-duplication rooting never flipped a margin (0/56); the failure needs a high-turnover lineage.
- **Overlap with Parsons et al. 2026.** That preprint studies tagging correctness and ASTRAL-Pro consistency under DLCOAL. If it already contains a hidden-paralogy counterexample under GDL, the novelty shrinks to the exact formula and threshold law. Check first.
- **ASTRAL-multi near-ties.** On both pools ASTRAL-multi's limiting margin is within ±0.003 of zero (not significant). A reviewer may ask whether Legried et al.'s theorem covers branch-specific rates [M: I believe their model uses uniform rates]. That needs checking, since a failure there would be a separate, larger claim.
- **Approximations.** The formula covers duplications on one branch only. The general case needs nested patterns (root branch, duplications on several internal branches).
- **Tools.** The ASTRID-multi and DISCO+ASTRID results use re-implementations (the sibling's code); DISCO itself is the original. DupTree is missing.

## 9. Reproduce

Set up the tools as described in `code/README.md` (micromamba env, tool builds, ASTER patch). Then:

```
python code/predict_margin.py 'A=1,4,1;B=4,8,1;C=0.5,1,4;x=0,0,4;y=8,0.5,0.5' 0.1 0.5 1
python code/curve4.py 'A=1,4,1;B=4,8,1;C=0.5,1,4;x=0,0,4;y=8,0.5,0.5;D=0,0,1' 40000 out.jsonl own:0,ovl:0,true:0 100,500,2000,10000
python code/run_curve.py gdl-crit 0 out.jsonl all      # after python code/gen_data.py gdl-crit 0
```
