---
title: "CS581 Project Proposal: Support-Aware Merging in MAGUS — Trusting Only the Evidence the Backbones Agree On"
author: "Evan Lin · [NetID] · October 2026 (draft; pilot numbers provisional)"
---

**Problem.** MAGUS (Smirnov & Warnow, *Bioinformatics* 2021) aligns large datasets by divide-and-conquer:
1. Split the sequences into subsets and align each with MAFFT L-INS-i.
2. Merge the subset alignments with the Graph Clustering Merger (GCM).

GCM learns how the subsets line up from ten *backbone* alignments of 200 sequences each. Every pair of residues that a backbone aligns becomes an edge between subset-alignment columns, weighted by how many backbones contain it. MCL clusters this graph, and a trace turns the clusters into the final alignment.

GCM treats every piece of backbone evidence the same, whether all ten backbones agree on it or only one does. The MAGUS authors list "modifying how GCM defines the weights on the pairs of columns" as future work (Zaharias, Smirnov & Warnow, TCBB 2022). The course project list asks to "improve MAGUS if possible", including by modifying GCM. We ask whether GCM is hurt by low-support evidence, and whether a simple, data-type-independent rule that keeps only evidence the backbones agree on gives better alignments and better trees.

**Pilot findings (done).**

(1) *Mechanism.* We scored GCM's cross-subset edges against the reference alignment. Edges contributed by only 1–3 of the 10 backbones are almost all wrong: on BAliBASE, 2–10% of their residue pairs are correct, vs 76–88% for edges with support ≥ 4. MAGUS still uses all of them.

(2) *Recipe.* The recipe was chosen on 10 training replicates, with the train/held-out split fixed before any results. It has two parts:
- delete cross-subset edges supported by fewer than 4 of the 10 backbones;
- down-weight (×0.03) backbone residue pairs that a fast second MAFFT run (FFT-NS-2) of the same backbone does not confirm. This costs about 1% of the backbone time.

Paired against MAGUS on identical subsets and backbones, over 24 datasets (one MAGUS run each):
- proteins −2.76 SP-error points (9/1/0, p = 0.002): AliSim simulations −5.3, BAliBASE −1.1;
- DNA/RNA (ROSE, RNASim, 16S) +0.01 (6/4/4, worst +0.47);
- held-out datasets alone: proteins −4.17 (4/0/0), DNA/RNA +0.06.

For scale, MAGUS's own gain over PASTA on the paper's 1,000-sequence data is about 2.7 points.

(3) *Replication.* An independent, pre-registered study on 28 held-out protein sets (simulated, HomFam, 10AA, fresh BAliBASE draws) found the support threshold alone (≥ 5 of 10, no second aligner) gives −1.90 (p = 0.004) and −0.03 on nucleotides.

(4) *Limits.*
- The gain is mostly on simulated proteins and through recall (fewer missed homologies); HomFam is flat (−0.1).
- *Trees do not improve.* A pre-registered test on 12 simulated protein sets (FastTree, nRF vs the true tree) found:
  - recipe minus MAGUS −0.15 RF points (6/2/4, p = 0.64);
  - on the hard sets −0.40, although the recipe lowers SP error by about 5 points there;
  - support threshold alone −0.03.
- *Where the tree headroom is.* On full-length DNA, even the *true* alignment improves ML trees by only about 1 FN point. On moderate proteins, MAGUS trees already match true-alignment trees. Only on hard proteins is there room: MAGUS trees have 9.3–12.1% RF error vs 5.2–9.7% with the true alignment.
- *Prior course project.* A Spring 2025 CS581 project (P. Srinivasan, "Evaluating Tree Estimation Error of MAGUS Alignments") measured FastTree error of MAGUS vs PASTA alignments on DNA only (ROSE 1000M1/M4, RNASim 10K). It found that a larger MAGUS subset size slightly improved trees while lowering alignment accuracy, so alignment and tree accuracy can disagree. We extend that question to proteins and to a change in GCM itself.

(5) *Other MAGUS components (the instructor's suggestions), same paired protocol.*
- Replacing MCL with other clustering (Leiden, Louvain, label propagation, connected components, constrained agglomeration; 22 datasets) gives nothing on MAGUS's raw graph. Once the support filter is applied, the clustering method matters ≤ 0.5 points. So the damage is in the graph, not in MCL.
- The filter also makes GCM's trace 2–5× faster.
- Swapping the subset aligner (MUSCLE5, ProbCons, FAMSA2, Clustal Omega) helps BAliBASE slightly (−0.3 to −0.5) but hurts simulated proteins and HomFam (+1 to +9.5). A Spring 2025 CS581 project (P. Srinivasan, "Evaluating Tree Estimation Error of MAGUS Alignments") measured FastTree error of MAGUS vs PASTA alignments on DNA only (ROSE 1000M1/M4, RNASim 10K). It found that a larger MAGUS subset size slightly improved trees while lowering alignment accuracy, so alignment and tree accuracy can disagree. We extend that question to proteins and to a change in GCM itself.

**Research questions.**
1. Does support-aware merging improve MAGUS alignments across the MAGUS paper's benchmarks and HomFam, with several MAGUS runs per dataset, measured end to end (accuracy and runtime)?
2. Does it improve maximum-likelihood trees (FastTree, IQ-TREE) where there is headroom (hard proteins), and if not, why? Measure which errors it fixes (SPFN vs SPFP), and whether those errors fall in phylogenetically informative columns.
3. Why does low-support evidence hurt, for example by fragmenting MCL clusters? Can a principled weight (each edge weighted by a function of its support) replace the hard threshold and the second aligner?
4. Does the rule carry over to larger inputs (Recursive-MAGUS / 10,000+ sequences) and to more backbones?

**Approach and evaluation.** The implementation is a thin wrapper around MAGUS's own GCM code, already built and validated: our MAGUS runs reproduce the paper's published error rates (median difference 0.02 points). Data:
- the MAGUS paper's ROSE 1000-sequence conditions, RNASim and 16S;
- BAliBASE RV100, HomFam (Homstrad-seed scoring) and 10AA;
- AliSim protein simulations with true trees.

Metrics: SP error, SPFN, SPFP and TC (FastSP); wall-clock and CPU time; tree error (normalized RF) of FastTree and IQ-TREE trees vs the true tree on simulated data. All comparisons are paired on identical subsets and backbones, with several MAGUS runs per dataset (MAGUS's run-to-run noise is about 0.5 points), W/T/L and Wilcoxon tests. Baselines: MAGUS, MAGUS with 20 backbones, PASTA, and MAFFT L-INS-i where feasible. For reproducibility, every dataset will be cited by its public source (Illinois Data Bank DOIs), and every table
regenerated by one script from raw outputs in the project repository.

**Timeline.**
- Week 1: replicate the recipe with 3 MAGUS runs per dataset on all benchmarks; end-to-end timing.
- Week 2: tree accuracy (FastTree + IQ-TREE) on simulated data, and the error analysis (RQ2).
- Week 3: mechanism and principled support weighting (RQ3); larger datasets (RQ4).
- Week 4: write-up.

**Risks.**
- (a) Trees did not improve in the pilot (−0.15 RF, n.s.). The contribution is therefore better alignments plus an explanation of why a ~5-point alignment gain does not reach the trees, which is itself an open question in the MSA lectures ("Which alignment criteria are predictive of tree accuracy?"). It is not a claim of better trees.
- (b) The gain is concentrated on simulated proteins (BAliBASE −1.1, HomFam flat); we will report per-benchmark results, not only averages.
- (c) Novelty is partial. Filtering by agreement is an old idea in consistency-based alignment (T-Coffee/M-Coffee, TCS). What is new is applying it to GCM's alignment graph with one data-type-independent rule, and the measured mechanism.

**Course connection.** MAGUS, GCM and MCL, and SP-error evaluation are taught in the MSA lectures. The instructor's project list has a "Projects related to MAGUS" section ("The goal here is to improve MAGUS if possible"), whose suggestions include the boldfaced "Modify the MAGUS design by using a different clustering method instead of the Markov Clustering algorithm within GCM". We modify GCM's input graph instead of its clustering method. The method keeps MAFFT as the only aligner.
