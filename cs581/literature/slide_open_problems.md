# Open questions in the CS581 Fall 2026 slides, turned into 4-week projects

Prepared 2026-10-09. Sources are the slide decks linked from the lectures page
(text in `notes/lectures/*.txt`; the linguistics deck `MIT2016-warnow-histling` was
added today). Each item quotes the slide, then gives a concrete project, the
data, what this repo already has, and an honest guess at the chance of a clear
result. Ideas already being piloted are marked **[running]**.

## Where the questions come from

| deck | slide | question (verbatim or near-verbatim) |
|---|---|---|
| First day | "Open problems (and possible course projects)" | MSA: merging two alignments; consensus alignments. Supertrees at 100,000+ species; approximation algorithms. Species trees/networks: heterogeneity from multiple causes incl. HGT; theoretical guarantees; scalability. New models of tree shape and sequence evolution; identifiability/consistency; evaluating methods under new models. Deep learning for large-scale tree and network estimation. Applications (linguistics, microbiome, protein function, tumors). |
| Phylogenomics part 2 | slide 35 "(Some) Open Questions" | Which species-tree methods are consistent under GDL / DLCOAL? Is ASTRAL-Pro consistent under random rooting/tagging error? Sample complexity of ASTRID and NJst under the MSC? A distance correction for GDL so ASTRID/NJst are consistent? Can we improve DISCO? Better quartet amalgamation? Scale concatenation? What about networks? Also: "Unknown if distance-based species tree estimation (e.g., ASTRID-multi) is statistically consistent under GDL models." |
| Divide-and-conquer trees | summary | "GTM does NOT allow blending, it is unlikely GTM is the best that can be done. Open problem: Develop a better DTM approach that allows blending." |
| Maximum likelihood | take-home | "Current heuristics either have computational challenges on large datasets (RAxML, IQ-TREE) or are not as accurate (FastTree). We need better heuristics." |
| MSA results from data | criteria / observations | "Which alignment criteria are predictive of tree accuracy?" "How should we design MSA methods to produce best accuracy?" "Most alignment methods over-align (compressed alignments) -> over-estimated branch lengths, under-estimated insertions." "Tree accuracy is not that well correlated with alignment accuracy." "BAli-Phy is best on simulated data, but not great on biological data. Why?" |
| Networks (CAMUS) | future work | Scale beyond 200 species; "develop other techniques to estimate the underlying tree T"; more complex than level-1. |
| Historical linguistics | future research | Realistic parametric models of linguistic evolution (polymorphism, homoplasy, borrowing, heterotachy) and methods developed and tested under them. |

## Candidate projects

### MSA (closest to what is already built)

1. **Learned evidence weights for GCM** (first-day "deep learning" x "merging alignments").
   MAGUS weights an alignment-graph edge by how many backbones put two columns
   together. Our oracle runs show that evidence quality is the main error source
   (true backbones: 8.23% -> 4.86% on 1000M2). Train a small classifier (gradient
   boosting or an MLP) to predict whether an edge is a true homology from cheap
   features: backbone vote count, mean HMMER posterior, column gappiness, subset
   disagreement, and the p-distance between the two subsets. Use the predicted
   probability (or its log-odds) as the edge weight. Training labels are free
   (true alignments of the ROSE replicates); test on held-out conditions,
   BAliBASE and 16S.
   *Have:* the harness, cached inputs for 59 replicates, PP-weighted graph code (`fastgraph`).
   *Cost:* no extra alignment runtime; graph size unchanged.
   *Chance:* ~35% for a clear gain over count weights (simple reweighting `weight-frac` was worse; PP weights helped a little).

2. **Do less-compressed alignments give better ML trees?** (MSA-results deck:
   "which alignment criteria predict tree accuracy", "over-alignment").
   Soft-constraint MAGUS produces alignments ~30% longer than MAGUS, with lower SPFP.
   Estimate FastTree/IQ-TREE trees on MAGUS, MAGUS(Slow), soft-MAGUS and PASTA
   alignments of the same replicates; correlate tree FN with SPFN, SPFP, TC,
   length ratio and compression. This connects the MSA work with the partner's ML
   interest. **[ML session has a tree harness on these datasets and a column-masking pilot running]**
   *Chance:* the evaluation always yields a result; ~40% that the less-compressed alignments give measurably better trees.

3. **Consensus alignments** (first day). **[running: 228 replicates]**
   Interim result: 0.78 points better than the input it picks as primary, but a
   tie with MAGUS(Slow) and 0.11 worse than the best input in hindsight. A
   variant worth one more try: use the collection to *annotate column
   reliability* (projects page) and mask unreliable columns before tree estimation.

4. **Merging two alignments** (first day). Exact two-alignment DP under the MWT
   objective is already implemented (`progressive.py`); it was not better than
   GCM for many-way merges. As a *pairwise* merger inside PASTA (instead of
   transitivity merging) it is untested. *Chance:* ~20%; PASTA is no longer the best pipeline.

### Trees: divide-and-conquer and ML

5. **DTM with blending** (divide-and-conquer deck). **[running: GTM session]**
   Interim: on the published RNASim1000 / Cox1-HET data, even the best possible
   merger of GTM's subset trees is only 0.1-0.5 FN points better than GTM, so the
   headroom there is small; a simulated caterpillar or unbalanced condition is
   where blending can matter.

6. **Better ML heuristics** (ML deck). E.g. FastTree -> limited IQ-TREE SPR
   polishing restricted to low-support regions; or GTM starting trees for
   IQ-TREE at 1000-10,000 taxa. **[ML session exploring]**

### Species trees (phylogenomics part 2 slide)

7. **Is ASTRID-multi consistent under GDL? + ASTRID-Pro** (two listed questions).
   Empirical consistency test: error of ASTRID-multi vs number of gene families,
   with true gene trees, at high duplication/loss rates. If error plateaus above
   zero, that is evidence of inconsistency; then the orthology-restricted,
   speciation-node-only distance ("ASTRID-Pro") as the fix. Data: FastMulRFS and
   DISCO data (published species trees). ~300 lines of Python; runs in seconds.
   *Chance:* ~50% to beat ASTRID-multi; a clean negative or theoretical finding is also reportable.

8. **Sample complexity of ASTRID/NJst under the MSC** (listed). Empirical curves:
   species-tree error vs number of genes, with true and estimated gene trees,
   ASTRID vs ASTRAL, across ILS levels; fit the rate. Uses published simulated
   data. An evaluation study (less "new method"), very low risk.

9. **Improve DISCO** (listed): re-root and re-tag gene-family trees by
   reconciliation against a first-pass species tree, then decompose
   ("DISCO-R"; see `wide_scan.md` D). *Chance:* ~35%.

10. **Better quartet amalgamation** (listed): crowded field (TREE-QMC, wQFM-TREE, wASTRAL). *Chance:* ~20%.

### Networks, models, applications

11. **CAMUS base tree** (networks deck: "other techniques to estimate the
    underlying tree T"): replace the ASTRAL base tree by a quartet-filtered or
    network-aware estimate. Data with true networks exist; no CAMUS outputs published. *Chance:* ~30%.

12. **Linguistic evolution models** (linguistics deck): simulate characters with
    polymorphism + homoplasy + borrowing and test tree methods. Rich but
    open-ended; data for Indo-European are on the histling page (not checked).

## Recommendation (not a commitment)

- If we stay in MSA: **#1 (learned evidence weights)** is the most "new method"
  and costs no runtime, which is what soft-MAGUS lacks. **#2** is a cheap
  companion that also uses the ML work.
- If we leave MSA: **#7 (ASTRID under GDL)** answers two questions printed on the
  same slide, is very cheap, and has published baselines.
