---
title: "CS581 Project Proposal: Soft-Constraint MAGUS — Merging All Sub-alignments at Once"
author: "[Your name] · [NetID] · October 2026 (draft)"
---

**Problem.** Divide-and-conquer aligners (SATé, PASTA, MAGUS) split the input sequences into
subsets, align each subset with MAFFT-L-INS-i, and merge the *subset alignments*, which are
kept as **hard constraints**. MAGUS [1] improved on PASTA [2] with one change: instead of
merging pairs of subset alignments and completing the merge by transitivity, its Graph
Clustering Merger (GCM) merges **all subset alignments at once**. GCM builds a graph whose
nodes are subset-alignment columns, weighted by how often 10 "backbone" alignments put their
letters together, clusters it with MCL and orders the clusters with A\*. MAGUS is still among
the most accurate methods for large, full-length nucleotide data (it ties UPP2 on 16S [3] and
TWILIGHT on RNASim-10K, and beats TWILIGHT on long-branch simulations [4]). It is not
state-of-the-art for proteins (learnMSA2, Muscle5, FAMSA2), for millions of sequences
(TWILIGHT, FAMSA2) or for fragmentary data (UPP2, WITCH). Warnow lists "the best way to
merge disjoint alignments" as an open problem [5].

**Pilot observation: the hard constraints, not the merge search, limit accuracy.** I replaced
MAGUS's inputs with the true alignment ("oracle" runs; identical subsets everywhere):
PILOT_ORACLE_SENTENCE
Errors *inside* one subset alignment propagate to every cross-subset pair, and no merge that
keeps the constraints can undo them. Under GCM's own objective (Maximum Weight Trace
alignment merging, MWT-AM, NP-hard [6]), the oracle merge scores *lower* than MAGUS's merge.
So searching the same objective harder does not help: GCM with the FM + optimizer search
[6] and my own exact progressive pairwise merger (dynamic programming per pair, as in
WITCH-NG [7]) change accuracy by at most a few tenths of a point.

**Proposed method: soft constraints by splitting.** Split each subset alignment into *m*
groups of similar sequences (average-linkage clustering on p-distance), induce one alignment
per group, and let GCM merge **all k·m groups at once**, using the same backbones. With
*m* = 1 this is exactly MAGUS; with one sequence per group it is MAGUS's unconstrained mode.
In between, backbone evidence can re-decide homology *between* groups of the same subset,
while the reliable within-group alignments stay fixed. This applies the PASTA → MAGUS idea
("merge everything at once") one level further down. It needs no changes to MAGUS: about
100 lines of Python prepare the inputs (already prototyped). I will study:

1. the number of groups *m* ∈ {1, 2, 3, 4};
2. similarity-based vs random grouping;
3. whether to add the original subset alignments as down-weighted extra evidence;
4. a simple adaptive rule that splits only divergent subsets (high average p-distance).

**Evaluation: apples-to-apples with the MAGUS paper.**

- *Datasets:* exactly those of MAGUS's main comparison (Illinois Data Bank
  doi:10.13012/B2IDB-2643961_V1):
  - ROSE 1000L1–L3, 1000M2–M4, 1000S1–S3 and RNASim-1K (1000 sequences each);
  - 8 BAliBASE protein sets and 16S.M;
  - stretch goals: 16S.3, 16S.T, RNASim-10K.
- *Settings:* MAGUS(Fast) with the paper's flags (25 subsets, 10 backbones × 200 sequences,
  MCL inflation 4, minclusters; 100 subsets for 16S.3/16S.T, as in the published runs). The
  new method reuses **the same subsets and backbones** per replicate, so only the merge
  changes (paired design).
- *Criteria:* SPFN, SPFP and their average via FastSP [8], as in the paper; also total-column
  score, alignment length relative to the reference (MAGUS is reported to over-inflate
  columns), and running time.
- *Statistics:* paired Wilcoxon signed-rank tests across replicates for each model condition.
  The *m* rule is tuned on replicates R0–R1 and tested on R2–R4.
- *Baselines:* the authors' published MAGUS(Fast), MAGUS(Slow) and PASTA alignments
  (rescored); our MAGUS rerun; GCM(fm+opt) [6].

**Validation already done.** Rescoring the authors' published alignments reproduces the
paper's figures to within 0.1 percentage points on every dataset checked (e.g. 1000L3: MAGUS
11.1% vs 11.0% in the paper, PASTA 17.5% vs 17.5%; 16S.T: 9.9% vs 10.0%). VALIDATION_RERUN_SENTENCE

PILOT_TABLE

**Four-week plan.**

| Week | Work |
|---|---|
| 1 | Finish baseline MAGUS reruns (5 replicates × 10 conditions + BAliBASE + 16S.M; scripted, about 40 CPU-hours); oracle error decomposition for every condition. |
| 2 | Sweep *m*, grouping and evidence weighting on R0–R1; fix one rule. |
| 3 | Held-out evaluation on R2–R4, BAliBASE and 16S.M; statistical tests; 16S.T/16S.3 if time allows. |
| 4 | Analysis (which conditions benefit: rate of evolution, gappiness, column inflation) and write-up. |

*Scope fallback:* if splitting does not help, the deliverable is a per-condition breakdown of
MAGUS's error (subset alignment vs backbone vs merge), which is still new and informative.

**Related work.**

- [1] Smirnov & Warnow, MAGUS, *Bioinformatics* 37:1666 (2021).
- [2] Mirarab et al., PASTA, *J Comput Biol* 22:377 (2015); Liu et al., SATé, *Science* 324:1561 (2009).
- [3] Park et al., UPP2, *Bioinformatics* 39:btad007 (2023).
- [4] Tseng, Walia & Turakhia, TWILIGHT, *Bioinformatics* 41:i332 (2025).
- [5] Warnow, "Advances in large-scale MSA" (talk, 2025).
- [6] Zaharias, Smirnov & Warnow, MWT-AM, *IEEE/ACM TCBB* 20:1700 (2023).
- [7] Liu & Warnow, WITCH-NG, *Bioinform Adv* 3:vbad024 (2023); Shen, Park & Warnow, WITCH, *J Comput Biol* 29:782 (2022).
- [8] Mirarab & Warnow, FastSP, *Bioinformatics* 27:3250 (2011).
- [9] Smirnov, Recursive MAGUS, *PLOS Comput Biol* 17:e1008950 (2021).
- [10] Kececioglu, maximum weight trace, CPM (1993).
- [11] Wheeler & Kececioglu, Opal, *Bioinformatics* 23:i559 (2007).
- [12] Edgar, Muscle5, *Nat Commun* 13:6968 (2022).
- [13] Becker & Stanke, learnMSA2, *Bioinformatics* 40:ii79 (2024).
- [14] Gudyś et al., FAMSA2, *Nat Biotechnol* (2026).
