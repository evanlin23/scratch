# Forest+DTM pilot: do DMR forest components + GTM beat NJ / FastME at short sequence lengths?

CS581 (Fall 2026) project pilot, candidate idea B from `cs581/literature/wide_scan.md` §2.1.
Everything here is reproducible from `cs581/forest/code/` (see "Reproducing" at the end).

## 1. The question and where it comes from

* Divide-and-conquer deck (`581-Divide-and-Conquer-trees-2023`), summary slide: *"The Guide Tree
  Merger (GTM) is the current leading DTM technique, based on empirical performance … However, GTM
  does NOT allow blending … Open problem: Develop a better DTM approach."*
* Lectures 1–2 (statistical consistency, sequence-length requirements, absolute fast converging
  methods): NJ can need exponentially long sequences; AFC methods (short quartets, DCM-NJ) only
  trust short distances.
* Kim, Lokhov, Vuffray, Romero-Severson & Goldberg, *BMC Bioinformatics* 27:186 (2026),
  doi:10.1186/s12859-026-06488-y, implemented the Daskalakis–Mossel–Roch forest algorithm
  (*SIAM J. Discrete Math.* 2011, doi:10.1137/09075576X) in `lanl/distphylo`. Forest returns only the
  reliably reconstructable part of the tree as leaf-disjoint components; it beat NJ only in limited
  cases. Their Discussion proposes, as future work, feeding forest components to a DTM.

**Pilot question.** If we take the forest components as constraint trees and merge them with GTM
under an NJ or FastME guide tree (a modern DCM-NJ), do we get a *full* tree that is more accurate
than NJ / FastME at short sequence lengths, in particular in hard (long-branch, deep) regimes?

## 2. Prior art (≈20 min search; details and DOIs in `lit/prior_art.md`)

* **Not done.** Forest+DTM appears only as a suggestion in Kim et al. 2026 (0 citing works found on
  Semantic Scholar / OpenAlex). None of ~50 works citing DMR 2011 merges forest components with a
  DTM or supertree method.
* Closest ancestors: **DCM-NJ / DCM1** (Huson, Nettles & Warnow, *JCB* 1999,
  doi:10.1089/106652799318337): threshold graph → overlapping small-diameter subsets → NJ per
  subset → strict-consensus merge; **DCM2** (Huson, Vawter & Warnow, ISMB 1999; no DOI found);
  **DCM-NJ+MP** (Nakhleh et al., *Bioinformatics* 2001, doi:10.1093/bioinformatics/17.suppl_1.S190);
  **Short Quartet Method** (Erdős, Steel, Székely & Warnow 1999,
  doi:10.1002/(SICI)1098-2418(199903)14:2<153::AID-RSA3>3.0.CO;2-R and 10.1016/S0304-3975(99)00028-6).
  These use *overlapping* subsets and consensus/dyadic-closure merges, not a disjoint forest + DTM.
* Other "reliable partial tree" methods, none of which merges components: Daskalakis et al. RECOMB
  2006 (doi:10.1007/11732990_24), Mossel TCBB 2007 (doi:10.1109/TCBB.2007.1010), Gronau, Moran &
  Snir RSA 2012 (doi:10.1002/rsa.20372), Mihaescu, Hill & Rao Algorithmica 2013
  (doi:10.1007/s00453-012-9644-4), Brown & Truszkowski AMB 2012 (doi:10.1186/1748-7188-7-32).
* DTMs: NJMerge (doi:10.1186/s13015-019-0151-x), TreeMerge (doi:10.1093/bioinformatics/btz344),
  Constrained-INC (doi:10.1186/s13015-019-0136-9), GTM (Smirnov & Warnow, *BMC Genomics* 2020,
  doi:10.1186/s12864-020-6605-1), Park, Zaharias & Warnow, *Algorithms* 2021 (doi:10.3390/a14050148).

So novelty is not the risk; usefulness is.

## 3. Reproduction / validation of the Forest implementation

`lanl/distphylo` (commit cloned 2026-10-09) is a research script: it hard-codes cluster paths,
shells out to R for scoring, and loops in pure Python with networkx graph copies (8–15 s per
(m, M, τ) point at n = 128). I re-implemented the algorithm (`code/forest_fast.py`) with the same
steps (clustering graph d < m → MiniContractor per leaf pair with 2τ phi-gaps → Extender →
tree popping), using integer bitmasks, one sort per pair shared by all τ, and an exact pairwise
compatibility check.

1. **Exact equivalence on small data** (`code/validate_small.py`, `results/validation_small.tsv`):
   distphylo's own `mini_contractor`/`extender`/`get_unique_trees` (imported unmodified) vs mine on 4
   simulated datasets (n = 24–30, short / long / deep / rate-heterogeneous), 102 (m, M, τ) points,
   64 of them with the ball smaller than the component (Extender active) and 66 with several
   components: **split sets identical on 102/102**, including the validity flag.
2. **distphylo's shipped output** (n = 128, k = 500 alignment `aln_128_500_1.fa`,
   `grid_summary_ntips128_1_k500_sorted.tsv`; `results/validation_distphylo_*.tsv`):

   | (m, M, τ) | published (#splits / false / #comp) | distphylo code rerun | ours |
   |---|---|---|---|
   | (0.8, 2.3, 0.10) | 22 / 0 / 1 | 22 / 0 / 1 | 22 / 0 / 1 |
   | (0.4, 1.0, 0.045) | 49 / 0 / 15 | 49 / 0 / 15 | 49 / 0 / 15 |
   | (0.4, 1.0, 0.03) | 77 / 0 / 15 | 77 / 0 / 15 | 77 / 0 / 15 |
   | (0.5, 1.2, 0.035) (their best) | 99 / 0 / 5 | — | 99 / 0 / 5 |
   | (0.6, 2.0, 0.05), (0.8, 2.0, 0.08) | absent (no forest) | invalid | invalid |
   | (0.6, 1.4, 0.05), (0.45, 1.0, 0.03) | 92 / 1 / 2, 92 / 0 / 9 | — | invalid (93 splits, incompatible) |

   The last row is the only divergence: distphylo's randomized tree popping (10 shuffles of
   `Tree.from_split_bitmasks`) silently drops one conflicting split, whereas I reject incompatible
   split sets. Speed: 0.02–0.17 s vs 0.4–15 s per grid point.
3. **Paper claims reproduced** (n = 128, k = 500 alignment): NJ, BIONJ, FastME and FastTree all
   recover all 125 splits (paper: "NJ … yields all 125 correct splits"); the Forest is
   less resolved — their dense, true-tree-informed grid gives at most 99/125 correct splits;
   my true-tree-free grid search (below) gives 95/125, all correct, in 4 s.
4. **Kim et al.'s simulation result**, short-branch regime (birth–death topology, branch lengths
   U[0.005, 0.05], JC, n = 50): see the "Kim et al. style" table in §5.1 — Forest returns no
   incorrect split more often than NJ at short k, but needs much longer sequences to be fully
   resolved, which is the qualitative finding of their Figs. 3–5.

## 4. Method

**Simulation.** Birth–death topologies (λ = 1, μ = 0.5; dendropy), unrooted; sequences with AliSim
(IQ-TREE 3.1.4). Regimes (n = 100 unless stated, 20 replicates per cell, a new tree per replicate):

| regime | branch lengths | model | why |
|---|---|---|---|
| short (n = 50) | U[0.005, 0.05] | JC | Kim et al. "short" |
| long | U[0.05, 0.1] | JC | Kim et al. "long" |
| deep | U[0.1, 0.4] | JC | many saturated pairs (16–37 % at k ≤ 3000) |
| het | ultrametric BD tree, height 2, × lognormal(0, 1) rate per branch | K2P (κ = 4), K2P distances | long branches next to short ones |

k ∈ {100, 300, 1000, 3000, 10 000, 100 000} (short: 100–5000). Large trees: n = 500, deep and het,
k ∈ {300, 3000}, 5 replicates, reduced grid.

**Distances.** JC (or K2P) from the alignment. Undefined (saturated) distances are replaced by
2 × the largest finite distance for the distance methods (cap2); the Forest uses the raw matrix
(saturated = ∞, never below m or M). A follow-up (§5.3) checks this choice.

**Methods.** NJ (FastME `-m N`), BIONJ (`-m I`), FastME (`-m B -n B -s`, balanced ME + NNI + SPR),
FastTree 2.1.11 (`-nt`, `-gtr` for K2P; skipped at k = 100 000 because it takes ~4 min per run),
and:

* **Forest** with hyperparameters chosen without the true tree, following Kim et al.'s Discussion:
  τ = ½ × {2, 5, 10, 15, 20, 30, 40, 50, 60, 70}th percentiles of the NJ tree's internal branch
  lengths; m = b × {1.0001, 0.9, 0.8, 0.7, 0.6, 0.45, 0.3}, where b is the bottleneck edge of the
  minimum spanning tree (b·1.0001 = smallest m giving one component); M ∈ {2m + 3τ_max + 10⁻³,
  1.5 × that} (Theorem 1 needs M > 2m + 3τ, m > 3τ). Among valid forests (all components
  compatible), keep the one with the most splits. Scored as a partial tree: FN = 1 − (correct
  splits)/(n − 3), FP = false/(forest splits).
* **Forest+GTM(guide)**, guide ∈ {NJ, FastME}: each component with ≥ 4 leaves keeps all its forest
  splits; its polytomies are refined with the guide method re-run on the component's sub-matrix
  ("sub", DCM-style; topped up with induced guide splits), or with the full guide tree's induced
  splits ("ind"). (GTM's loader resolves polytomies in constraint trees arbitrarily, so they must be
  binary first.) Singletons, pairs and triples are passed as they are. GTM (`vlasmirnov/GTM`,
  default `convex` mode) then merges them under the guide tree.
* **Controls** that separate the contribution of the Forest:
  *Forest comps only + GTM* = the same components, forest splits dropped, refined by subset NJ;
  *Centroid dec. + GTM* = the standard DTM pipeline with no Forest: centroid-edge decomposition of
  the guide tree into subsets of ≤ 25 leaves, subset trees by the same method, GTM.

**Statistics.** Per cell, paired by replicate: mean difference in FN rate, W/T/L with tie band
|Δ| < 0.5/(n − 3) (less than one split), two-sided Wilcoxon signed-rank p (zero differences
dropped; p = 1 when all differences are zero). No multiple-testing correction; with 20 replicates
p < 0.01 is the bar I treat as a real difference.

RESULTS_PLACEHOLDER
