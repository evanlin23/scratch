# Do alignment criteria predict ML tree accuracy? Do less-compressed (soft-MAGUS) alignments give better trees?

Pilot for a CS581 (Fall 2026) project idea (`literature/slide_open_problems.md`, item 2).
Branch `claude/cs581-alncrit`; code in `code/`, result files and plots in `results/`.

## 1. Questions and where they come from

From the deck *MSA results from data* (`notes/lectures/MSA-results-from-data-2026.txt`):

* slide "Alignment criteria": *"Which alignment criteria are predictive of tree accuracy?"* and *"How should we design MSA methods to produce best accuracy?"*
* slide "Observations – II" (Löytynoja & Goldman 2008): *"Most alignment methods 'over-align' (produce compressed alignments)… Compression results in over-estimations of branch lengths, under-estimation of insertions"*
* slide "Results so far" / "Observations" (Nelesen et al. 2008): *"Tree accuracy is not that well correlated with alignment accuracy."*

The concrete question for this repo: the soft-constraint MAGUS variant `slow-soft-m3` has lower SP error than MAGUS (about −0.6 to −0.75 points of (SPFN+SPFP)/2) and produces alignments that are about 20–30% longer, with lower SPFP. Does that give better ML trees?

## 2. Prior art (details and DOIs: [`prior_art.md`](prior_art.md))

* **Nelesen et al. 2008, PSB** (doi:10.1142/9789812776136_0004): better guide trees barely changed SP error (SPFN) but improved RAxML trees. FTA alignments with *worse* SPFN gave better trees. "Not all errors are of equal importance."
* **Liu, Linder & Warnow 2010, PLoS Currents** (doi:10.1371/currents.RRN1198): "alignment error, measured using SP-FN, is not particularly predictive of tree error". Prank+GT had high SPFN but good trees; Opal was the reverse.
* **Liu et al. 2009 Science (SATé)** (doi:10.1126/science.1171243) and **SATé-II 2012** (doi:10.1093/sysbio/syr095): better alignments led to better trees. Likelihood with gaps treated as missing data is not a usable alignment criterion.
* **Liu & Warnow 2012 PLoS ONE** (doi:10.1371/journal.pone.0033104): treelength-optimal alignments are worse both as alignments and for trees.
* **Mirarab et al. 2015 PASTA** (doi:10.1089/cmb.2014.0156) and **Nguyen et al. 2015 UPP** (doi:10.1186/s13059-015-0688-z): UPP has lower SP error but lower TC than PASTA, and PASTA gives better trees. TC ranked the methods the right way for trees; SP did not.
* **Nute & Warnow 2016** (doi:10.1186/s12864-016-3101-8): TC gains tracked tree gains better than SP gains.
* **Nute, Saleh & Warnow 2019 Syst Biol** (doi:10.1093/sysbio/syy068): introduces the expansion ratio (estimated length / true length). PRANK's long alignments gave *poor* trees on hard protein simulations.
* **Smirnov & Warnow 2021 MAGUS** (doi:10.1093/bioinformatics/btaa992): reports alignment error only, no trees.
* **Löytynoja & Goldman 2008 Science** (doi:10.1126/science.1158395): aligners over-align. The paper evaluated indel and branch-length inference, not topology.
* **Wang et al. 2011 TCBB** (doi:10.1109/TCBB.2009.68) and **Ogden & Rosenberg 2006 Syst Biol** (doi:10.1080/10635150500541730): alignment accuracy and tree accuracy are positively but loosely correlated. The link is strongest on hard data and depends on tree shape.
* **Dessimoz & Gil 2010** (doi:10.1186/gb-2010-11-4-r37) and **Tan et al. 2015** (doi:10.1093/sysbio/syv033): gaps carry signal, and removing columns usually hurts trees.
* **Chowdhury & Garai 2017 Genomics** (doi:10.1016/j.ygeno.2017.06.007): a catalogue of alignment objective functions; no tree evaluation.

**Gap.** We found no study that separates SPFN from SPFP (over- vs under-alignment) as predictors of ML tree error, or that tests whether less-compressed alignments give better ML trees. The existing evidence is anecdotal and points both ways (Prank+GT in 2010 vs PRANK in 2019).

## 3. Data and methods

* **Published alignments**: MAGUS paper Results.zip (Illinois Data Bank, doi:10.13012/B2IDB-2643961_V1), streamed with the HttpFile range reader (`code/fetch_published.py`). Per replicate it contains `true_align.txt`, `true_tree.tre`, `gcm.txt` (MAGUS), `gcm_slow.txt` (MAGUS(Slow)), `pasta_align.txt` (PASTA, 3 iterations), `pasta_1_align.txt`, `pasta_4_align.txt`, `pasta_3_gcm_align.txt` (PASTA(3)+GCM), `pasta*.tre` and timing/log files.
  * Coverage: 241 replicate directories over ROSE 1000L1–3, 1000M1–4, 1000S1–3 (20 each; 1000M1 has 19 in the scores file), RNASim 1K (20), RNASim 10K (10), 16S.3/16S.T/16S.M/16S.B.ALL (1 each) and BAliBASE (8). `gcm_slow` is missing for RNASim10K and 16S.B.ALL.
  * 1000M1 additionally has about 70 MAGUS parameter-sweep alignments (`gcm_<K>_<…>[_hmm].txt`, 19 reps). These are unused here but make a natural within-condition design for a full project.
  * The full listing is `results/results_zip_listing.txt`.
* **True trees**: `true_tree.tre` in Results.zip, identical (RF = 0, checked on 1000M2/R0) to `rose.tt` in Datasets.zip (`/opt/data/Datasets/ROSE/<cond>/R<i>/rose.tt`). It is the model tree with zero-length edges contracted (992 internal edges on 1000M2/R0). `random.tree` is the binary model tree. RNASim: `/opt/data/Datasets/RNASim/1000/R<i>/true_tree.tre`.
* **Trees**: FastTree 2.1 `-nt -gtr -gamma -nosupport` (`code/runtrees.py`). `-nosupport` and `-fastest` gave the same topology as the default on a test alignment. IQ-TREE 3.1.4 `--fast -m GTR+G4` was run on the soft-MAGUS set as time allowed.
* **Tree error**: FN rate = missing true bipartitions / true internal bipartitions (`code/treeerr.py`). It matches the ML session's dendropy `treeerr.py` exactly (e.g. 0.1069 for MAGUS on 1000M2/R0, identical to their baseline).
* **Alignment criteria**: FastSP SPFN, SPFP, TC and compression (= estimated length / true length; reused from `experiments/validation/published_scores.jsonl` for published files, and computed by `gcmx.experiment` for reruns). Also indel events (internal gap runs) relative to the true alignment (`code/alnstats.py`).
* **Soft-MAGUS reruns**: 21 ROSE/RNASim replicates whose MAGUS (paper settings) subalignments and backbones were cached by the fanout workers (`inputs.tar.xz` on `claude/cs581-worker-*`) (`code/run_soft.sh`).
  * Three merges on identical inputs: `default` (MAGUS merge; reproduces the full MAGUS run per validation layer 3), `slow` (HMM-extended backbones = MAGUS(Slow) evidence) and `slow-soft-m3`.
  * These are the same development replicates on which soft-m3 was first measured, not held-out ones.
* **Analysis**:
  * Within-replicate (replicate fixed effects) correlations, partial correlations and OLS with cluster-robust SEs (`code/analyze_pub.py`).
  * Paired Wilcoxon tests (`code/analyze_soft.py`).
  * A union analysis with up to 9 alignments per replicate (`code/analyze_union.py`).
  * Controlled perturbations of the true alignment (`code/perturb.py`, `code/analyze_perturb.py`).

## 4. Criteria vs tree error on the published alignments

* **Coverage**: 41 replicates × 6 alignments (TRUE, MAGUS, MAGUS(Slow), PASTA(3), PASTA(1), PASTA(3)+GCM), 3–4 replicates per condition over all 10 ROSE conditions and RNASim 1K. That is 246 FastTree trees; full tables are in [`results/criteria.md`](results/criteria.md) and [`results/pub_table.csv`](results/pub_table.csv).
* **Run order and coverage limit**: trees were run replicate by replicate across conditions, so coverage stays balanced. 1224 trees were queued; the time budget allowed 251.

**Means over all 41 replicates**

| alignment | SPFN % | SPFP % | compression | FastTree FN % | FN − MAGUS (paired) | W/T/L (tie 0.1) | Wilcoxon p |
|---|---|---|---|---|---|---|---|
| TRUE | 0 | 0 | 1.000 | 9.45 | −1.05 | 34/1/6 | 1.4e-06 |
| MAGUS | 7.21 | 6.53 | 1.001 | 10.50 | – | – | – |
| MAGUS(Slow) | 6.68 | 6.32 | 1.007 | 10.56 | +0.06 | 20/1/20 | 0.91 |
| PASTA(3) | 9.36 | 9.15 | 0.869 | 10.78 | +0.28 | 18/3/20 | 0.14 |
| PASTA(1) | 9.58 | 9.35 | 0.872 | 11.08 | +0.58 | 17/1/23 | 0.050 |
| PASTA(3)+GCM | 7.97 | 7.50 | 0.894 | 10.73 | +0.23 | 21/0/20 | 0.59 |

Per-condition values are in `results/criteria.md`. W/T/L counts are wins/ties/losses relative to MAGUS.

**Observations**

* PASTA's alignments are strongly compressed (about 13% shorter than true; 29% on 1000L3) and have about 2.4 points more SP error than MAGUS. Yet their trees are only 0.2–0.6 FN points worse.
* Every estimated alignment costs about 1–1.6 FN points relative to the true alignment. That gap is the total headroom any aligner could win back for FastTree.

**Which criterion predicts tree error?** Fixed effects per replicate (so only differences between alignments of the same data count), with cluster-robust SEs:

| | SPFN | SPFP | (SPFN+SPFP)/2 | TC error | log compression |
|---|---|---|---|---|---|
| within-replicate r, incl. TRUE (n=246) | 0.647 | 0.646 | 0.647 | 0.514 | −0.288 |
| within-replicate r, estimated only (n=205) | 0.475 | 0.475 | 0.477 | 0.377 | −0.068 (n.s.) |
| within R², estimated only | 0.226 | 0.226 | 0.228 | 0.142 | 0.005 |
| pooled Spearman ρ, no fixed effects (estimated) | 0.552 | 0.557 | 0.555 | **0.665** | 0.010 |

* **SPFN vs SPFP cannot be separated with these data.** Within a replicate SPFN and SPFP correlate at r = 0.99, because real aligners make both kinds of error together.
  * Partial correlations: r(FN, SPFN | SPFP) = 0.05 and r(FN, SPFP | SPFN) = 0.05, both n.s.
  * In two-predictor models the coefficients are unstable and their signs flip between subsets. Example: soft set estimated-only gives SPFN −0.62, SPFP +1.01; the union set gives 0.06 / 0.20, n.s.
  * Any claim that "SPFP matters more" from such a regression would be an artefact of collinearity.
* **SP error is the best within-replicate predictor; TC is the best across conditions.**
  * Within a replicate, the SP measures (either one or their mean) explain about 23% of the variance in tree-FN differences between estimated alignments, or 42% with TRUE included. TC explains less.
  * Without replicate fixed effects (pooled), TC ranks highest (ρ = 0.67 vs 0.56), consistent with Nguyen 2015 and Nute 2016. TC mainly tracks how hard the condition is.
* **Compression on its own does not predict tree error.**
  * At *fixed* SP error, a longer alignment is *slightly worse*: log-compression coefficient +3.6 FN points per unit of log length ratio (p = 0.007); partial r = +0.28 (p = 4e-5, estimated alignments).
  * This is the opposite of "less compressed is better".
* **The association is weak within conditions.** Per-condition within-replicate r(FN, SP error) ranges from −0.29 (1000M3) to +0.83 (1000L1), with 3–4 replicates per condition (`results/criteria.md`).
  * This reproduces the slide's "tree accuracy is not that well correlated with alignment accuracy".
  * FastTree's sensitivity to small alignment changes adds about ±1 FN point of noise per replicate.

Plots:
* [`results/within_rep_scatter.png`](results/within_rep_scatter.png): FN vs SPFN, SPFP and log compression, demeaned within replicate.
* [`results/fn_by_method.png`](results/fn_by_method.png): mean FN per condition and method.

**Union analysis** ([`results/union.md`](results/union.md), [`results/union_scatter.png`](results/union_scatter.png)): 21 replicates in both sets, 6 published alignments plus 3 reruns, n = 160 estimated alignments.
* Same picture: SP r ≈ 0.52 (R² 0.27), TC R² 0.20, log compression R² 0.05 (n.s.).
* SPFN–SPFP correlation is still 0.98 even with the longer soft-m3 alignments included.

## 5. Controlled perturbations: is all alignment error equally bad for trees?

Because observational data cannot separate the two error types, I perturbed the TRUE alignment of the same 21 replicates (`code/perturb.py`) and ran FastTree on the results:

* `split_r`: with probability r, a column is split into two halves. This is FN-only under-alignment and makes the alignment longer.
* `shift_0.1`: a residue next to a gap swaps into that gap. In ROSE this lands mostly in sparse columns, so it is FN-dominated.
* A third perturbation, merging compatible adjacent columns (FP-only over-alignment), produces negligible SPFP (< 0.01%) on ROSE: the only compatible columns are sparse insertion columns. Pure-FP error is not realistically constructible here, which is part of why SPFN and SPFP co-vary.

Results ([`results/perturb.md`](results/perturb.md), [`results/perturb_scatter.png`](results/perturb_scatter.png)), 21 replicates:

| alignment | SPFN % | SPFP % | compression | FN excess over TRUE tree (pts) | p (≠ 0) | excess per SPFN point |
|---|---|---|---|---|---|---|
| split_0.1 | 4.89 | 0 | 1.054 | +0.06 | 0.47 | 0.014 |
| shift_0.1 | 8.39 | 0.19 | 0.997 | +0.35 | 0.14 | 0.066 |
| split_0.2 | 9.75 | 0 | 1.107 | +0.49 | 0.042 | 0.051 |
| MAGUS (rerun) | 6.12 | 5.55 | 0.997 | +0.63 | 0.014 | 0.119 |
| MAGUS(Slow) | 5.61 | 5.35 | 0.997 | +0.71 | 0.005 | 0.144 |
| slow-soft-m3 | 5.25 | 4.91 | 1.185 | +0.62 | 0.023 | 0.130 |

* **At equal SPFN, real aligner errors cost trees 2–9× more than random FN-only errors.** split_0.1 has about the same SPFN as MAGUS but gives better trees: −0.57 FN points, W/T/L 13/3/5, p = 0.020. MAGUS at 6% SPFN costs as much as random splitting at 10% SPFN.
* **Two explanations are confounded.** The extra cost of real errors may come from their FP component (over-alignment) or from their being *systematic* (similar sequences mis-aligned the same way, which creates false shared states). This pilot cannot tell them apart.
* **The extra length itself is nearly free.** split_0.1 makes the alignment 5% longer and costs +0.06.

## 6. Soft-MAGUS trees

21 replicates (ROSE L1–3, M2–4, S1–3, RNASim; cached MAGUS inputs) with merges on identical inputs. Full table: [`results/soft_trees.md`](results/soft_trees.md).

| | MAGUS | MAGUS(Slow) | slow-soft-m3 |
|---|---|---|---|
| SPFN / SPFP % | 6.12 / 5.55 | 5.61 / 5.35 | 5.25 / 4.91 |
| (SPFN+SPFP)/2 % | 5.84 | 5.48 | 5.08 |
| compression (est/true) | 0.997 | 0.997 | 1.185 |
| FastTree FN % (21 reps) | 9.53 | 9.62 | 9.53 |
| IQ-TREE --fast FN % (13 reps) | 8.95 | 9.10 | 8.98 |

**Alignment error.** slow-soft-m3 − MAGUS = −0.75 points, W/T/L 19/0/2 (tie 0.05), p = 2.4e-5. The improvement is confirmed.

**Tree error.** Two-sided Wilcoxon; tie band |ΔFN| ≤ 0.1 points:

| comparison | tree method | n | mean ΔFN (pts) | W/T/L | p |
|---|---|---|---|---|---|
| slow-soft-m3 − MAGUS | FastTree | 21 | **−0.00** | 9/2/10 | 0.90 |
| slow-soft-m3 − MAGUS(Slow) | FastTree | 21 | −0.09 | 12/2/7 | 0.37 |
| MAGUS(Slow) − MAGUS | FastTree | 21 | +0.09 | 7/3/11 | 0.42 |
| MAGUS − TRUE | FastTree | 21 | +0.63 | 14 worse / 2 / 5 better | 0.014 |
| slow-soft-m3 − MAGUS | IQ-TREE --fast | 13 | +0.02 | 5/4/4 | 1.0 |
| slow-soft-m3 − MAGUS(Slow) | IQ-TREE --fast | 13 | −0.12 | 4/5/4 | 0.64 |
| MAGUS − TRUE | IQ-TREE --fast | 13 | +0.21 | 7 worse / 0 / 6 better | 0.79 |

* **No tree gain from soft-m3.** Its alignments are 19% longer than MAGUS's and have 0.75 points less SP error. Yet the trees are *exactly* as accurate as MAGUS trees with FastTree (mean difference 0.00), and within 0.02 points with IQ-TREE --fast.
* **The headroom is small.** With FastTree, MAGUS trees are only 0.63 points worse than TRUE-alignment trees on these replicates, and 0.21 points with IQ-TREE --fast (n.s.). A 13% cut in SP error would be worth about 0.08 FN points at best if tree error scaled with SP error. The perturbation results suggest the error soft-m3 removes is the cheap kind.
* **MAGUS(Slow) does not help trees either** (+0.09, n.s.), although its SP error is 0.36 points lower.

## 7. Runtime (4 cores; tree jobs ran 4–6 at a time, so per-job times are inflated by about 1.5×)

* **FastTree** GTR+G, single thread: mean 99 s per tree on ROSE 1000 taxa (median 91 s) and 206 s on RNASim 1K. 382 trees total.
* **IQ-TREE 3 --fast** GTR+G4, 1 thread: mean 183 s. 52 runs.
* **Soft-MAGUS per replicate on cached inputs:**
  * HMM extension of the 10 backbones: 118 s (4 jobs).
  * Merges, 3 concurrent at 1 thread each: MAGUS 43 s, MAGUS(Slow) 64 s, slow-soft-m3 200 s.
  * Total: about 5–8 min per replicate.
  * A full MAGUS run (decomposition + MAFFT) on the worker machines took 1352 s per replicate; it was skipped here by reusing cached inputs.
* **Wall clock**: about 3.6 h of compute (22:50–02:25 UTC) plus analysis. All published-alignment files were streamed from Results.zip in about 10 min.

## 8. Verdict

**The soft-MAGUS → better-trees hypothesis is not promising.** The alignment gain is real (−0.75 SP points, p = 2e-5), but it gives zero tree gain. The result holds with both FastTree (n = 21) and IQ-TREE --fast (n = 13), on the same replicates where the alignment gain was measured. The headroom against true-alignment trees is only about 0.2–0.6 FN points. A 4-week project built on "soft-MAGUS gives better trees" would very likely end negative.

**The broader criteria question is unclear, but it has a promising angle.** The pilot reproduced the textbook statements:
* Within replicates, SP error explains about 23% of tree-error differences between estimated alignments.
* TC is better across conditions.
* Compression is uninformative, or slightly the wrong way.
* SPFN and SPFP are inseparable in real alignments (r = 0.99).

The new result is the perturbation contrast: at equal SPFN, real aligner errors cost 2–9× more tree accuracy than random FN-only errors. This suggests a project about *which* alignment errors matter for trees, rather than about lowering SP error.

**What a 4-week project would look like**
1. **Week 1: more data.** Finish the published-alignment tree set: all 20 replicates × 6 methods × 11 conditions (about 1300 FastTree runs; about 36 core-hours, done in a day on a campus cluster). Add the 1000M1 MAGUS parameter sweep (about 70 alignment variants × 19 replicates). It gives many alignments of the same data whose SPFN, SPFP and length vary more independently, which is the best observational lever on the SPFN-vs-SPFP question.
2. **Week 2: designed error types**, separating random vs systematic and FN vs FP error:
   * clade-structured mis-alignment: shift the same residues in all members of a clade;
   * random shifts;
   * FP created by merging *incompatible* columns with resolution of displaced residues;
   * errors concentrated in long vs short branches.
   Measure the tree cost per unit SPFN and SPFP for each type.
3. **Week 3: tree-aware criteria.** Propose and test criteria that weight homology errors by their phylogenetic impact. Candidates: SP error counted only for pairs from different sides of true-tree edges, or the error rate among parsimony-informative columns. Test whether they predict tree error out of sample better than SP or TC, held out by condition.
4. **Week 4: write-up and one design implication.** For example: do the cheap error types justify spending less effort on the parts of MAGUS that remove them?

**Risks**
* FastTree tree noise (about ±1 FN point per replicate) forces many replicates.
* The true-alignment headroom is only about 1 FN point, so effects are small.
* Results may be specific to ROSE indel models; RNASim behaves differently (compression 1.43 for MAGUS).
* A "tree-aware" criterion needs the true tree, so it is an analysis tool, not an objective an aligner can optimize directly.
* IQ-TREE results are on 13 replicates only, and IQ-TREE shrinks the alignment effect (MAGUS − TRUE was 0.21 points, n.s.).

## 9. Files

* `code/fetch_published.py`: streams selected files from Results.zip.
* `code/runtrees.py`, `code/treeerr.py`: trees and FN/FP.
* `code/run_soft.sh`, `code/soft_final.sh`: soft-MAGUS alignments and their trees.
* `code/perturb.py`, `code/run_perturb.sh`: controlled perturbations.
* `code/analyze_pub.py`, `code/analyze_soft.py`, `code/analyze_union.py`, `code/analyze_perturb.py`: the analyses.
* `code/alnstats.py`: indel counts.
* `results/*.jsonl`: raw per-tree and per-alignment results.
* `results/*.md`: full tables.
* `results/*.png`: plots.
* `prior_art.md`: literature notes with verified DOIs.
* Large files (alignments and trees) are under /opt/data/published, /opt/runs/{alncrit,soft,perturb} and are not committed.
