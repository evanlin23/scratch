# Pilot: is ASTRAL-Pro consistent under rooting/tagging error? Which methods are consistent under GDL / DLCOAL?

*CS581 project pilot, 2026-10-10, about 2.5 h wall-clock on 4 cores. Branch `claude/cs581-gdlcons`. Code is in `code/` (see `code/README.md`); results are in `results/`. The theory note is `results/theory.md`; the prior-art table with DOIs is `results/prior_art.md`.*

## 0. Verdict: **promising** (for question 1); unclear for question 2

**Question 1 (ASTRAL-Pro under rooting/tagging error).** The answer depends on the error model, and that dependence can be made exact. On a 4-taxon caterpillar I derived the limiting ASTRAL-Pro quartet scores in closed form, as an integral over duplication times with exact birth–death generating functions. The formula matches simulation to within standard errors. It gives three results:

1. **The unmodified ASTRAL-Pro3 binary is statistically inconsistent under pure GDL, even with true gene trees.**
   - Configuration: one very supercritical short branch above a three-taxon clade, followed by heavy, unequal loss.
   - ASTRAL-Pro3 returns the wrong tree in 4/4 datasets of 10,000 families, 20/20 of 2,000 and 4/4 in an 8-taxon embedding.
   - Given the *true* tags it converges to the right tree (wrong in 0/24 datasets of ≥ 2,000 families, plus 0/4 in the embedding). Below 500 families the orthologous signal is too sparse.
   - Cause: hidden paralogs. These are duplications whose two sides become species-disjoint after complementary losses. Any species-overlap rule, including ASTRAL-Pro's own, must label them speciations, and here they systematically favour a wrong quartet.
   - This does not contradict Zhang et al.'s theorem, which assumes correct tags. It shows the assumption is not just technical.
2. **Under random errors that respect ASTRAL-Pro's tagging rule** (mislabel each hidden paralog with probability q), there is a sharp threshold q* = O / (H_wrong − H_correct).
   - O is the expected number of orthologous quartet classes; the H terms are the expected numbers of hidden-paralog classes of each topology.
   - Predicted q* = 0.079; observed sign change between 0.10 and 0.15, and the predicted margins match within 1 SE at five q values.
3. **How common?** The exact formula on 20,000 random rate configurations finds only **12 of 15,434** inconsistent with overlap tags, almost all with λ ≫ μ on one branch.
   - On all standard-looking settings (16- and 50-taxon Yule trees, λ = μ, the sibling's adversarial rates) every rooting/tagging error model that respects Def. 1 converged to zero error by 50–200 families. On the generic small-tree grids, their quartet margins stayed ≥ 0.51.
   - "Naive" random tag flips that create S labels on nodes whose children share species are a different story. They break ASTRAL-Pro at q ≈ 0.035 even where overlap tagging is fine. That model is outside what ASTRAL-Pro's tagger can output, so the slide question needs an error model that respects Def. 1.

**Question 2 (which other methods are consistent).**
- In the generic atlas (true gene trees, 4 GDL and 4 DLCOAL settings with SimPhy, up to 5,000 families), every method converges, or nearly so. No clear plateau appears. STAG is slowest; ASTRID-multi has small residual errors under the highest-duplication DLCOAL setting.
- On the adversarial 4-taxon pool, methods split by mechanism:
  - ASTRAL-Pro, ASTRAL-DISCO and ASTRID-DISCO use ASTRAL-Pro-style tagging and are all fooled (100% wrong at 2,000 families).
  - ASTRID-multi and FastMulRFS look near-random.
  - **STAG is always right**, because its closest-copy distances ignore hidden paralogs.
- The literature has proofs only for ASTRAL-one/multi (GDL, DLCOAL), ASTRAL-Pro and ASTRAL-DISCO (GDL, *correct tags*), and FastMulRFS (non-adversarial GDL). DLCOAL proofs for tag-based methods exist only as recent, possibly shaky preprints or papers (Section 2).

**Why "promising".**
- The project has a clean, new, checkable core: an exact formula, an explicit counterexample in the shipped software, and a threshold law.
- There is a natural theorem to aim for: consistency of ASTRAL-Pro with overlap tags when no branch is supercritical.
- It reuses the sibling session's simulator and fits in 4 weeks (Section 7).

**The main risk.** The counterexample needs extreme rates, so a reviewer may call it pathological. The closest prior work is Parsons et al. 2026 (DLCOAL tagging correctness). I could not read its theorems, and it must be checked before committing to the project.

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
- On 16-taxon settings (critical λ = μ = 1, high 3/3, the sibling's adversarial (3, 0.5)/(0, 3), and a root duplication burst) and 50-taxon settings, every ASTRAL-Pro error model reaches FN = 0 by 50–200 families.
- Models covered: true, ovl, rovl 0.1/0.3/1, flip 0.05/0.15/0.3, d2s 0.5, s2d 0.5.
- The only non-zero entry is naive all-S tagging (d2s q = 1) on gdl-adv, at 50–200 families. It fixes itself by 1,000.

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

__ATLAS__

### 5.2 Adversarial 4-taxon pools

Fraction of disjoint datasets on which each method returns a wrong 4-taxon tree. True gene trees in every case; ASTRAL-multi is the ASTER control.

__ATLAS4__

## 6. Theory derived

`results/theory.md`:
1. **The error-as-selection view.** Rooting and tagging only select which quartets count. Their unrooted topologies are fixed.
2. **The one-sidedness of overlap-tag errors** on true roots under GDL (D→S on hidden paralogs only).
3. **The exact 4-taxon limiting-score formula** with its consistency condition O + H_AB > max(H_AC, H_BC), and the threshold q* = O / (max H_wrong − H_AB) for random hidden-paralog mislabelling.
4. **A refutation of the easy proof route.** The per-node association inequalities (I1, I2) that would give a general proof fail on 4.6% (I1) and 7% (I2) of random survivable rate configurations.

All of these are checked numerically; none are written up as formal proofs yet.

## 7. What a 4-week project would look like

| week | work |
|---|---|
| 1 | Write the 4-taxon proposition rigorously: class decomposition, the integral formula, conditioning on the species filter, and the version with a root branch. Read Parsons et al. 2026 in full; get the PDF from a library or the authors. |
| 2 | **Positive result.** Try to prove ASTRAL-Pro with overlap tags is consistent under GDL when duplication is bounded, for example λ_e < μ_e with bounded turnover λ_e·t_e. "No supercritical branch" alone is *not* enough: the one non-supercritical failure in the scan is critical with huge turnover (λ = μ = 8 on a branch of length 2; formula only, not simulated). Use the formula to map the inconsistency region in (λ_y·T, loss asymmetry). Extend the formula to random rooting (`rovl(p)`). |
| 3 | **n-taxon and empirical side.** Show inconsistency for n > 4 (the 8-taxon embedding already works empirically). Rerun the atlas on the adversarial region, including estimated gene trees from sequences, to see whether the effect survives gene-tree error. Run under DLCOAL with SimPhy, where tag errors go both ways. |
| 4 | Write-up. Optional: a "tag-robust" ASTRAL-Pro variant that downweights S-tagged nodes whose children have very unbalanced copy numbers, or uses closest-copy information as STAG does. Test whether it fixes cand2 without hurting generic accuracy. |

**Deliverables either way:** an exact formula, an explicit counterexample in shipped software, a threshold law, and an atlas.

## 8. Risks

- **Pathological regime.** The counterexample needs λ ≈ 8 on one branch and heavy loss (0.08% of random configurations). The honest framing: correct tagging is *necessary* for the theorem, and overlap tagging is safe in all but extreme supercritical regimes. Proving the "safe" part is the harder half.
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
