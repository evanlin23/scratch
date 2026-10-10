---
title: "CS581 Project Proposal: Consensus Evidence for MAGUS — Cleaner Alignment Graphs for Protein Alignment"
author: "Evan Lin · [NetID] · October 2026 (draft, pilot numbers provisional)"
---

**Problem.** MAGUS (Smirnov & Warnow, *Bioinformatics* 2021) aligns large datasets by divide-and-conquer:
it splits the sequences into subsets, aligns each with MAFFT L-INS-i, and merges the subset alignments with
the Graph Clustering Merger (GCM). GCM's only information about how the subsets relate comes from ten
*backbone* alignments of 200 sequences (eight per subset), also aligned with L-INS-i: every pair of residues
that a backbone aligns becomes a weighted edge in an alignment graph, which MCL clusters and a trace turns
into the final alignment. GCM trusts every backbone pair equally, and the MAGUS authors list "modifying how GCM
defines the weights on the pairs of columns" as future work (Zaharias, Smirnov & Warnow, TCBB 2022). We ask
whether GCM is limited by the *quality of its evidence* and whether cheap, reference-free consensus between
alignments can clean it, keeping MAFFT as the only aligner.

**Pilot findings (done).** (1) *Mechanism.* On the MAGUS paper's BAliBASE protein sets, most of MAGUS's error
is false positives, and across 194 merge variants the change in evidence precision on cross-subset residue
pairs predicts the change in false-positive rate (Spearman ρ = −0.66). Residue pairs made by only one aligner
are mostly wrong on 7 of 8 sets (4-29% correct). (2) *Method.* Keeping only the backbone pairs that both
L-INS-i and MAFFT FFT-NS-2 (`--op 3`, ~1.5 CPU-s per backbone) align, and merging as usual, lowers MAGUS's
error by 1.83 points on BAliBASE (7/1/0 sets, one MAGUS run each, paired, Wilcoxon p = 0.016); masking columns
the ten L-INS-i backbones disagree on (no new alignment at all) gives −1.73 (7/1/0, p = 0.008). For scale,
MAGUS's own gain over PASTA on the paper's 1,000-sequence data is 2.7 points. (3) *Held-out test (pre-registered,
interim).* On simulated proteins with true alignments (AliSim LG+G4 with indels, 1,000 sequences, two divergence
levels, 8 replicates) the MAFFT-only intersection lowers error on 8/8 replicates (mean −2.27 points, range −0.89
to −4.17; p = 0.008); the masking variant is weaker (−0.84, 5/8). On HomFam (10 families, Homstrad-seed scoring)
it is flat (−0.15, 6/1/3, n.s.; −1.60 to +2.09), where MAGUS is limited by recall rather than precision; over all 20
held-out sets −0.97 (15/1/4, p = 0.007). (4) *Limits.* Hard filtering hurts on nucleotide data (ROSE, RNASim: +7 to
+31), where pairs found by only one alignment are mostly correct, and alignment gains do not yet show up in FastTree
trees (8 simulated sets, mean RF −0.1, within noise). Adding Clustal Omega backbones hurts simulated proteins,
and a pre-registered reference-free gate did not transfer. (5) *A recipe for all data types (chosen on 10
training sets, no data-type switch).* Cross-subset GCM edges supported by fewer than 4 of the 10 backbones are
almost all wrong on every data type (2-10% correct vs 76-88% at support ≥ 4). Down-weighting unconfirmed pairs to
0.03 and deleting edges with support < 4 (recipe fixed on 10 training replicates) gives, over all 24 replicates,
proteins −2.76 (9/1/0, p = 0.002; BAliBASE −1.08, AliSim −5.28) and DNA/RNA +0.01 (6/4/4, worst +0.47); on the
held-out replicates alone, proteins −4.17 (4/0/0) and DNA/RNA +0.06 (4/3/3). It costs one FFT-NS-2 run per
backbone (~1% of L-INS-i). The support threshold alone, with no second aligner, gives proteins −2.03 (6/3/1)
and DNA/RNA −0.03 (worst +0.05), but ≈ 0 on BAliBASE. On DNA/RNA the recipe is safe but does not help, because
MAGUS's DNA error is recall-limited.

**Research questions.** (1) Does consensus evidence improve MAGUS across protein benchmarks (BAliBASE with
fresh draws, HomFam, 10AA, simulated proteins), and on which data does it hurt? (2) What explains when it helps:
measured evidence precision vs recall per dataset, and can a reference-free statistic (backbone agreement,
TCS-style consistency, MUMSA overlap) predict it, including for nucleotides? (3) Which consensus is best per unit
of runtime: intersection with one cheap MAFFT mode, edge-support thresholds (≥ k of 10 backbones), consistency
masking, or more backbones? (4) Do the better alignments give better ML trees on simulated data?

**Approach and evaluation.** Implementation is a thin wrapper around MAGUS's own code (already built: backbone
aligner swap, merge-only reruns on identical subsets and backbones, consistency masking). Data: the MAGUS paper's
BAliBASE RV100 sets, HomFam (Homstrad-seed scoring), 10AA, AliSim protein simulations with true trees, plus
ROSE/RNASim/16S nucleotide controls. Metrics: SP error (FastSP), SPFN and SPFP separately, TC, wall-clock and
CPU end to end, and FastTree/IQ-TREE tree error on simulations. Every comparison is paired (same subsets and
backbones, several MAGUS draws per dataset; MAGUS's draw-to-draw noise is ~0.5 points), with Wilcoxon tests.
Baselines a reviewer will ask for: MAGUS default and MAGUS(Slow), more or larger backbones, a global edge-support
threshold, TCS-filtered backbones, M-Coffee and MAFFT L-INS-i on the full data where feasible, PASTA.

**Timeline.** Week 1: finish the held-out test (HomFam, fresh BAliBASE draws, trees); end-to-end runtime.
Week 2: the when-does-it-help analysis and a reference-free switch (protein vs nucleotide, or agreement-based).
Week 3: alternatives and baselines (edge-support threshold, TCS, more backbones, M-Coffee). Week 4: write-up.

**Risks.** (a) The effect is data-dependent: large on simulated and BAliBASE proteins, flat on HomFam, harmful on
DNA when filtering is hard. Unless a weighting or switch works across all of them, the claim is "precision-limited
protein data", and the when-does-it-help analysis (RQ2) becomes the core of the project. (b) Novelty is partial: consistency scores (GUIDANCE, trimAl, TCS) and
multi-aligner consensus (M-Coffee) are known; what is new is using them to clean GCM's alignment-graph evidence,
and the precision mechanism, which we will state only for the data where we measure it. (c) MCL already exploits
transitive consistency, so gains could shrink with tuned MCL; we test that directly.

**Course connection.** MAGUS, GCM and divide-and-conquer alignment are core course material, the instructor's
project list suggests "study MAGUS with different base methods", and consistency-based alignment (T-Coffee) and
SP-score evaluation are covered in the MSA lectures. The method keeps MAFFT as the only aligner.
