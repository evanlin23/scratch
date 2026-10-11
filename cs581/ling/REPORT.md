# Pilot: stochastic models of linguistic character evolution and testing tree methods under them

CS581 Fall 2026 project-idea pilot (candidate #12 in `../literature/slide_open_problems.md`).
Branch `claude/cs581-ling`. Time spent: ~4.5 h wall-clock. Nothing here is a decision; this is
evidence for or against spending 4 weeks on the idea.

## 1. Question and where it comes from

The linguistics lecture deck (`../notes/lectures/MIT2016-warnow-histling.txt`) ends with:

* **Future research:** "more investigation of statistical methods based on good stochastic
  models" and "realistic parametric models of linguistic evolution and method development
  under these parametric models".
* **Modelling issues:** units, polymorphism, homoplasy, non-treelike evolution (borrowing),
  heterotachy, deviation from the lexical clock. The deck notes that the WERN model has homoplasy
  and borrowing but no polymorphism; the Nicholls–Gray stochastic Dollo model has polymorphism but no
  borrowing; and Gray–Atkinson binary encoding has both polymorphism and homoplasy but no borrowing.
* **Main points:** linguistic data evolve differently from molecular data, and all methods "need
  to be carefully tested". The deck also criticises binary (presence/absence) encoding of
  multistate characters.

Pilot question: can a small, honest simulator with all six features, plus paired method
comparisons, give a project with a new and defensible result?

## 2. Prior art (details, DOIs and verification status in `results/prior_art_notes.md`)

| work | DOI | simulation? | gist |
|---|---|---|---|
| Ringe, Warnow & Taylor 2002, TPS | 10.1111/1467-968X.00091 | no | Screened IE data (24 languages, 294 chars); weighted maximum compatibility tree, ~95% of characters compatible |
| Nakhleh, Ringe & Warnow 2005, *Language* | 10.1353/lan.2005.0078 | no | Perfect phylogenetic networks: 3 contact edges explain the residual incompatibilities |
| Nakhleh, Warnow, Ringe & Evans 2005, TPS | 10.1111/j.1467-968X.2005.00149.x | no | Methods × datasets on IE; all methods except UPGMA find the major subgroups, Anatolian+Tocharian and Greco-Armenian |
| Warnow, Evans, Ringe & Nakhleh 2006 (CUP chapter) | none found | model | WERN model: homoplasy plus borrowing on networks, no polymorphism |
| Nicholls & Gray 2008, JRSS-B | 10.1111/j.1467-9868.2007.00648.x | — | Stochastic Dollo (TraitLab) |
| Ryder & Nicholls 2011, JRSS-C | 10.1111/j.1467-9876.2010.00743.x | — | SD with missing data |
| Kelly & Nicholls 2017, AoAS | 10.1214/17-AOAS1040 | yes | SD with lateral transfer (borrowing-aware Bayesian) |
| Greenhill, Currie & Gray 2009, Proc B | 10.1098/rspb.2008.1944 | yes | SD-simulated borrowing; topology robust to local borrowing, not to global borrowing (same model used to generate and analyse) |
| Bouckaert et al. 2012, *Science* | 10.1126/science.1219669 | — | Bayesian phylogeography; Anatolian origin |
| Barbançon et al. 2013, *Diachronica* | 10.1075/dia.30.2.01bar | **yes**, 3,584 datasets | WERN model, 30 taxa, 0–3 contact edges; UPGMA < NJ < Gray–Atkinson < MP/MC; weighting helps only on screened data. **No polymorphism, no ML.** |
| Chang, Cathcart, Hall & Garrett 2015, *Language* | 10.1353/lan.2015.0005 | — | Ancestry constraints; steppe dates |
| Neureiter et al. 2022, HSSC (contacTrees) | 10.1057/s41599-022-01211-7 | yes | Bayesian model with contact edges |
| Heggarty et al. 2023, *Science* (IE-CoR) | 10.1126/science.abg0818 | — | New IE cognate database |
| **Canby, Evans, Ringe & Warnow 2024, TPS** | **10.1111/1467-968X.12289** | **yes**, 2,560 datasets | **Adds lexical polymorphism to WERN.** Polymorphism-aware MP ("MP4") is best in almost all conditions; Gray–Atkinson binary (MrBayes) is second; SD, NJ and UPGMA are worse. Public Java simulator (LingPhyloSimulator) and a corrected 370-character IE dataset that keeps polymorphism. |
| King 2026, PLOS CB | 10.1371/journal.pcbi.1014312 | yes | Calibration and adequacy problems in the standard BEAST2 linguistic setup |

**Novelty:** Canby et al. 2024 already covers most of what this candidate proposed (polymorphism plus
borrowing plus homoplasy, simulated with paired method comparisons, from Warnow's own group). The gaps
that remain, as listed in the notes and confirmed by this pilot:
* ML on multistate vs binary encodings;
* how polymorphism is *coded* before analysis, holding the method fixed;
* data-driven borrowing filters;
* compatibility on polymorphic data;
* missing data;
* network methods scored against known contact edges;
* model adequacy of the generators against real IE summary statistics.

## 3. Data (reproduction)

* The CPHL page (http://tandy.cs.illinois.edu/histling.html) links the Ringe–Taylor files at Rice:
  `IEDATA_112603` (screened, 294 chars = 259 L + 13 M + 22 P) and `IEDATA_050704` (unscreened,
  375 chars, 53 marked `*` = removed by screening). Copied to `data/`.
* Lines marked `!` (22 P + 8 M) are the characters that RWT required to be compatible. Canby's CSV
  gives exactly these 30 characters weight 10000.
* Also fetched: Canby et al.'s corrected dataset (`data/IE_canby2024.csv`, polymorphic cells `a/b`).

**Reproduced properties** (`results/ie_repro.md`, `code/ie_repro.py`):
* **Compatibility on the reference tree.** On an RWT/NRW-style reference tree transcribed from the
  deck's figure (the resolution inside subgroups is my own guess), 277/294 = **94.2%** of the screened
  characters are compatible, and all 30 `!` characters are. The deck says "95%".
  * The best weighted-compatibility tree found (the `!` characters forced compatible) has 278/294
    compatible and is RF 4 from the reference.
  * Unweighted MC/MP find 280/294, but only by making 3 `!` characters incompatible. This is
    exactly why RWT forced those characters to be compatible.
  * The unscreened data: 307/375 compatible on the reference tree.
* **Subgroups.** MP, MC, WMC, capped parsimony, ML-Mk and ML-binary all recover the seven
  multi-language major subgroups (Anatolian, Tocharian, Indo-Iranian, Italic, Celtic, Germanic,
  Balto-Slavic), plus Greco-Armenian and (except ML-binary on the screened data) Anatolian+Tocharian.
  This matches the deck's claim. NJ misses Anatolian+Tocharian.
  * The Satem core is recovered only by the weighted methods and ML. This matches "the Satem core is
    not always reconstructed" and "the choice of method matters".
* **Polymorphic data.** Weighted MP on Canby's polymorphic dataset finds **two distinct equally
  optimal trees** (score 58, `results/canby_mp.md`), consistent with Canby et al.'s report of two
  optimal MP trees. I did not compare them topology-by-topology with the paper's figures.
  * Note: my Fitch treats a polymorphic cell as an ambiguity, not as Canby's MP4.

## 4. Simulator (`code/sim.py`)

Time-sliced simulation on a dated tree with contact edges.

**Tree and rates**
* Yule topology, 24 leaves (as in RWT). About 30% "ancient" leaves, made by cutting their terminal
  edges short.
* Lognormal edge-rate multipliers (dlc = 0.3) model deviation from the lexical clock.
* Gamma(α = 1) rates across characters.
* Lognormal per-(character, edge) multipliers (het = 0.3) model heterotachy.

**Lexical characters**
* A change creates a new cognate class, so there is no back-mutation.
* In 10% of characters a change goes, with probability ½, to a previously used class (parallel
  semantic shift).
* Polymorphism: on a monomorphic slot, a change adds a synonym with probability `p_poly`. On a
  polymorphic slot, a change is a loss with probability 0.6, otherwise a replacement.

**Morphological and phonological characters**
* 20% of morphological characters revert to an older state with probability ½.
* 30% of phonological characters are binary 0↔1 (natural sound changes: parallel and back mutation).
* All other morphological and phonological characters are infinite-state. Neither type is
  polymorphic.

**Borrowing**
* Each contact edge joins two coexisting, non-sister lineages over a window of length 0.15.
* Each character is borrowed (with a fixed direction per contact edge) with probability
  b_L = 0.10, b_P = 0.03 or b_M = 0.
* A borrowed word replaces the recipient's word, or is added as a synonym.

**Coding**
* Multistate methods see a *randomly resolved* word per polymorphic cell, as a lexicographer would
  record it.
* The binary encoding keeps all words present.

**Calibration against IE** (`results/sim_summary.md`)

| statistic | simulated | real IE |
|---|---|---|
| lexical states per character | 11.5–12.6 | 14.3 (RWT), 16.4 (Canby) |
| morphological states per character | 8.7–9.8 | 8.8–9.9 |
| phonological states per character | 2.2–2.8 | 2.3 |
| lexical characters with ≥1 polymorphic cell (moderate) | 0.81 | 0.67 (Canby) |
| polymorphic cells (moderate) | 0.18 | 0.086 |

**Polymorphism is about 2× too high per cell**, and higher still in "poly-high". This is the main
realism gap.

**Conditions** (`code/run_sim.py`): 6 conditions × 30 replicates.

| condition | settings |
|---|---|
| clean | no homoplasy, borrowing or polymorphism |
| moderate | 1 contact edge |
| borrow3 | 3 contact edges, b_L = 0.15 |
| poly-high | p_poly = 0.5 |
| homoplasy | homoplastic fractions 25% L / 50% M / 60% P, so the M/P weights are misspecified |
| lexonly | 250 lexical characters only, 2 contact edges |

Character counts: 250 L + 15 M + 25 P. Replicates 0–9 are **training** (used only to choose the cap
below); replicates 10–29 are **test**, giving 120 test datasets.

## 5. Methods

**Parsimony and compatibility (one search engine for all).** A numba Fitch search
(`code/fitch.py`) does SPR hill-climbing from an NJ start plus 8 random starts. It minimises
Σ_c w_c·g(extra_c), where extra_c is the number of extra steps for character c:

| method | g(e) | weights |
|---|---|---|
| MP | g(e) = e | — |
| MC | g(e) = [e > 0] | — |
| WMC | g(e) = [e > 0] | M and P characters ×5 |
| **Cap-k** (the pilot's new idea) | g(e) = min(e, k) | — |

* Ties are broken by parsimony (ε = 1e-4).
* Cap-k is "compatibility tolerating bounded homoplasy", a soft data-driven screening. It
  interpolates between MC (k = 1) and MP (k = ∞).
* **Search check:** the search never scored worse than the true tree under its own objective
  (48/48 checks, `results/search_check.txt`).
* TNT downloaded but was not run: the sandbox refused to execute the unvetted binary.

**Other methods**
* **NJ** on the −ln(1−p) corrected proportion of characters with different states.
* **ML-Mk**: IQ-TREE 2.0.7, `-st MORPH -m MK+G4`, multistate.
* **ML-bin**: IQ-TREE `-st BIN -m GTR2+FO+G4` on the Gray–Atkinson binary encoding. This is an ML
  stand-in for their Bayesian analysis.
* **MP-poly**: Fitch with polymorphic cells as state sets.

**Metric and tests**
* Tree FN rate. Trees are binary, so FN = FP = RF/(2(n−3)); FN moves in steps of 1/21 = 0.048.
* Two-sided paired Wilcoxon signed-rank tests on the test replicates. The tie band is exact ties
  (|ΔFN| < 1e-9), which are dropped; win/tie/loss counts are reported alongside.

## 6. Results (full tables: `results/sim_summary.md`, `results/supp_summary.md`)

**Mean FN rate on the test replicates** (20 per condition; lower is better)

| method | clean | moderate | poly-high | homoplasy | borrow3 | lexonly | all | sec/rep |
|---|---|---|---|---|---|---|---|---|
| NJ | 0.098 | 0.083 | 0.131 | 0.143 | 0.133 | 0.129 | 0.119 | 0.02 |
| MP | 0.057 | 0.100 | 0.140 | 0.143 | 0.152 | 0.133 | 0.121 | 0.6 |
| MP-poly (sets) | 0.055 | 0.129 | 0.190 | 0.212 | 0.171 | 0.212 | 0.162 | 0.5 |
| **MC** | 0.057 | 0.079 | 0.136 | 0.129 | 0.131 | 0.124 | **0.109** | 0.6 |
| WMC (M,P ×5) | 0.057 | 0.079 | 0.143 | 0.140 | 0.140 | 0.124 | 0.114 | 0.6 |
| Cap2 | 0.057 | 0.095 | 0.148 | 0.133 | 0.136 | 0.133 | 0.117 | 0.6 |
| Cap3 | 0.055 | 0.095 | 0.143 | 0.152 | 0.143 | 0.140 | 0.121 | 0.6 |
| ML-Mk | 0.086 | 0.117 | 0.157 | 0.157 | 0.138 | 0.148 | 0.134 | 64 |
| **ML-bin** | 0.083 | 0.076 | 0.086 | 0.086 | 0.095 | 0.093 | **0.087** | 5.0 |
| ML-bin, resolved data (supp.) | 0.083 | 0.093 | 0.143 | 0.117 | 0.112 | 0.129 | 0.113 | ~5 |
| MP on binary encoding (supp.) | 0.129 | 0.098 | 0.114 | 0.138 | 0.131 | 0.133 | 0.124 | ~1 |

**Key paired tests** (all 120 test replicates; ΔFN = A − B; W/T/L counts are A better / tie / A worse)

| A vs B | mean ΔFN | W/T/L | p |
|---|---|---|---|
| Cap2 vs MC (**new idea**) | +0.008 | 12/77/31 | 0.039 (Cap2 worse) |
| Cap2 vs MP | −0.004 | 30/69/21 | 0.30 |
| MC vs MP | −0.012 | 39/64/17 | 0.007 |
| WMC vs MC | +0.005 | 9/100/11 | 0.12 |
| MC vs ML-Mk | −0.025 | 58/36/26 | 0.001 |
| MC vs ML-bin | +0.023 | 36/29/55 | 0.005 (ML-bin better) |
| ML-bin vs ML-bin on resolved data | −0.026 | 57/50/13 | 4e-6 |
| ML-bin on resolved data vs MC | +0.004 | 36/38/46 | 0.80 |
| ML-bin vs MP on binary encoding | −0.037 | 63/37/20 | 4e-7 |
| MP-poly vs MP (randomly resolved) | +0.041 | 20/39/61 | 3e-6 |

**What the numbers say**

1. **The pilot's new idea (capped parsimony) is a negative result.**
   * On the training replicates the best cap was k = 1, which is plain MC (MC 0.094, Cap2 0.104,
     Cap3 0.102, MP 0.108).
   * On the test replicates Cap2 is slightly *worse* than MC (p = 0.04) and no different from MP.
   * Weighted caps were worse still.
   * Tolerating 1–2 extra steps per character lets borrowed and homoplastic characters back in.
     Under this model, throwing them out (MC) works better.
2. **MC ≥ MP**, significant overall (p = 0.007), which agrees with Barbançon et al. **WMC does not
   help** (100/120 ties); in the "homoplasy" condition the ×5 weights are wrong and WMC is slightly
   worse.
3. **ML-Mk is the worst character-based method** (0.134), even in the clean condition. Its cost is
   about 64 s per replicate, against 0.6 s for the parsimony searches. An equal-rates Mk model with
   up to 24 states is badly misspecified for "always a new cognate" characters.
4. **Surprise: ML on the binary encoding is best overall** (0.087). This reverses the
   Barbançon/Canby ranking, where Gray–Atkinson binary was behind MP. The supplement locates the
   cause:
   * Most of the gap is **polymorphism information**. Binary ML on the resolved data loses 0.026 FN
     (p = 4e-6) and then ties MC (p = 0.8).
   * Random resolution, as I coded it for the multistate methods, throws information away.
   * Fitch with polymorphic cells as state sets (MP-poly) is *worse* still (+0.041 over MP), because
     treating polymorphism as ambiguity is the wrong model.
   * The binary *model* also adds something over binary *parsimony* (p = 4e-7).
   * Caveat: my polymorphism rate is about 2× the real IE rate per cell, which inflates this
     effect. In the clean condition, with no polymorphism, MC/MP (0.057) beat ML-bin (0.083).
   * I did **not** implement Canby's MP4, the right polymorphism-aware competitor. So this is not
     evidence against their result; it does show that the coding step decides the ranking.
5. **NJ ≈ MP** overall (p = 0.92). NJ is clearly worse only in the clean condition. Barbançon et al.
   found NJ clearly worse. My trees may be more clock-like and less heterotachous than theirs.

## 7. Verdict: **unclear, leaning not promising as originally framed**

* **Against:**
  * The core idea (simulate polymorphism + homoplasy + borrowing and test tree methods) was done in
    2024 by Warnow's own group (Canby et al.), with a public simulator, 2,560 datasets and the
    corrected IE data. A 4-week replication-plus-perturbation risks being seen as derivative.
  * The pilot's new method (capped parsimony) did not beat MC.
* **For:** the pilot found one concrete, cheap and still-open question with an effect large enough
  to measure: **how polymorphism is coded and modelled before tree estimation.** Binary ML kept
  polymorphism and won; multistate methods given resolved data lost about 0.02–0.05 FN. The notes
  list this as untested gaps #1 and #2. Neither Barbançon nor Canby ran ML, and neither isolated the
  coding step.

**What a 4-week project would look like** (if chosen)
1. **Week 1.** Install LingPhyloSimulator (Canby) and generate their model conditions (30 taxa,
   0–3 contact edges, 5 polymorphism levels). Re-calibrate my generator to real IE polymorphism
   (about 8.6% of cells) as a second, independent generator.
2. **Week 2.** Hold the estimator fixed and vary the coding: drop polymorphic characters; random
   resolution; split coding; MP4 (Canby); binary encoding. Methods: MP/MC, ML-binary (IQ-TREE /
   RAxML-NG with +ASC), ML-multistate.
3. **Week 3.** One new method that is polymorphism-aware *and* multistate. Options: a compatibility
   criterion over word sets (a character is compatible if some choice of one word per cell is
   compatible; Bonet et al. 1999 give the theory), or a binary-ML-guided screen.
4. **Week 4.** Paired tests across both generators, re-analysis of Canby's IE data, write-up.

**Risks**
* Novelty against Canby 2024. Mitigation: frame the project as "coding vs method" and "ML vs
  parsimony", which they did not test.
* How realistic the simulation is (this pilot's polymorphism is too high; all generators are
  calibrated only to IE). Mitigation: two generators plus summary-statistic matching.
* Small tree sizes (24–30 taxa), so FN differences are in steps of about 0.05 and need ≥100 paired
  replicates.
* IQ-TREE Mk runtime (about 1 min per replicate). Binary ML is about 5 s.
* Data access is **not** a risk: every dataset was public and fetched in minutes.

## 8. Files

* `code/`
  * `trees.py`: trees, Newick, splits, NJ.
  * `fitch.py`: numba Fitch and SPR search with the generic objective.
  * `methods.py`: method wrappers.
  * `sim.py`: the simulator.
  * `run_sim.py`: the experiment (restartable).
  * `analyze.py`: summary tables and tests.
  * `supp_encoding.py`: the encoding supplement.
  * `ie_repro.py`, `canby_repro.py`: the IE reproductions.
  * `search_check.py`: the search-adequacy check.
* `results/`
  * `prior_art_notes.md`: the literature notes.
  * `ie_repro.md`, `ie_trees.nwk`, `canby_mp.md`: IE reproduction results.
  * `sim_runs.jsonl` (180 runs), `sim_summary.md`: main simulation results.
  * `supp_runs.jsonl`, `supp_summary.md`: the supplement.
  * `search_check.txt`: the search check.
* `data/`: the RWT screened and unscreened data (CPHL) and Canby 2024's `ie_dataset.csv`
  (LingPhyloSimulator repo).
* Reproduce: `cd code`, then run `python3 ie_repro.py`, `python3 run_sim.py 30 4`,
  `python3 analyze.py` and `python3 supp_encoding.py 4`. Needs numpy, scipy, numba and `iqtree2`.
