# GTM-Blend: validation of the published GTM numbers and a pilot of a blending disjoint tree merger

CS581 (Warnow, Fall 2026) project de-risking. All code is in `cs581/gtm/code/` and all result tables are in `cs581/gtm/results/`.
Large data and intermediate trees were kept outside the repository (`/opt/gtmdata`). `code/fetch_data.sh` re-downloads them.

## 0. Summary

| Question | Answer |
|---|---|
| Can we reproduce the published GTM-pipeline numbers? | **Yes.** Rescoring the published trees reproduces 29/30 table entries to ≤0.05 points (FN %); the 30th is a FastTree run (51.4 vs 50.9). Rerunning GTM ourselves on the published guide and subset trees gives the published GTM tree exactly (RF = 0) in 39/40 replicate×guide cases. |
| Is there room for a better DTM on the published conditions? | **Very little.** On RNASim1000 and Cox1-HET, GTM is within 0.1–0.5 points of the (optimistic) lower bound for *any* merger of those subset trees, blended or not. On 1000M1-HF, oracle-guided blending heuristics that know the true tree improve GTM by only 1.4 points (FastTree guide) / 0.6 points (IQ-TREE guide). |
| Does blending help at all? | **Yes, when the decomposition subsets are not clades of the true tree.** This happens when the guide tree is poor. In simulation (200 taxa, short internal branches, exact subset trees), the best *unblended* merge has 21–25% error, while blended merges with 0% error exist. |
| Does our blending merger (GTM-Blend-ML) beat GTM? | **In simulation, yes, in all 3 conditions, 60/60 replicates, every p < 1e-4.** Yule with exact subset trees: 28.0% → 13.7% FN. Caterpillar with exact subset trees: 26.6% → 16.5%. Yule with IQ-TREE subset trees: 50.2% → 44.1%. **On published 1000M1-HF:** FastTree guide 42.45% → 41.42% (4 better / 0 worse / 1 tied), IQ-TREE guide unchanged; pooled 5/0/5, p = 0.0625 (the minimum attainable p with 5 non-zero pairs). It never made a tree worse. |
| GTM vs TreeMerge (instructor's question)? | Published data: GTM is never worse. Simulation: on **caterpillar** trees TreeMerge trends better than GTM (−1.6 points, 11/5/4, p = 0.083, *not significant*); on Yule trees they tie. GTM-Blend-ML beats TreeMerge 20/20 in every simulated condition. |
| Does parsimony work as the score? | **No.** Pooled over the 40 published cases, constrained parsimony search lowers the parsimony score but makes the tree significantly *worse* than GTM: +0.21 points FN (31 of 40 worse), p = 2e-4. Parsimony insertion is worse still: +0.54 points, p = 1e-5. Maximum likelihood is required. |
| Caveat | At 200 taxa with 50-taxon subsets, full-data IQ-TREE (33.1%) beats every DTM pipeline that uses *estimated* subset trees (GTM 50.2%, GTM-Blend-ML 44.1%). The simulations isolate *merger* error; they are not a case for DTMs at this scale. |
| Viable 4-week project? | **Yes, if it is framed as "blending fixes GTM when the guide or decomposition is poor"**, with simulations as the main evidence and the published data as a "does no harm / bounded headroom" check. It is *not* viable as "beat GTM on the published RNASim1000/Cox1-HET conditions": the ceiling analysis shows that cannot happen. See §5. |

## 1. Which paper and data

The Data Bank deposit doi:10.13012/B2IDB-7008049_V1 belongs to **Park, Zaharias & Warnow, "Disjoint Tree Mergers for Large-Scale Maximum Likelihood Tree Estimation", *Algorithms* 14(5):148, 2021, doi:10.3390/a14050148**.
That paper evaluates GTM, TreeMerge and Constrained-INC inside ML pipelines.
It is *not* the original GTM paper (Smirnov & Warnow, *BMC Genomics* 21(S2):235, 2020, doi:10.1186/s12864-020-6605-1), which is about species trees.
The scan in `cs581/literature/wide_scan.md` conflates the two. All validation below is against Park et al. 2021, Tables 3–5 and A2.

Conditions used (the three small ones):

| Condition | Taxa | Sites | Reps | True trees | Notes |
|---|---|---|---|---|---|
| RNASim1000 | 1000 | 21,946 | 5 | in `RNASim1000.tar.gz` (`true.tt`) | 3 subsets of ≤500 |
| Cox1-HET | 2341 | 658 | 10 | in `Cox1-HET.zip` (`true-tree.tre`) | heterotachy; 7–8 subsets |
| 1000M1-HF | 1000 | 3,880 | 5 | **not in the deposit**: taken from the SATé data doi:10.13012/B2IDB-5139418_V1 (`sate_journal/1000M1.tar.bz2`, replicates R0–R4, `rose.tt`) | fragmentary sequences; 3 subsets |

`rose.tt` (the binary true tree) reproduces the published numbers, `rose.mt` is about 0.15 points off. So the paper evidently scored against `rose.tt`.

The pipeline in every case (paper §3.1):
- the starting tree (FastTree 2 or IQ-TREE 2) is the guide tree;
- the guide tree is split by centroid decomposition into subsets of at most 500 taxa;
- IQ-TREE 2 estimates the subset trees on the true alignment;
- a DTM merges them (TreeMerge/Constrained-INC use topological distances from the guide).

The error criterion is the FN (missing-branch) rate against the true tree.

## 2. Validation

### 2.1 Rescoring the published trees

`code/validate.py` uses our own bipartition code (`code/phylo.py`, cross-checked against DendroPy unrooted FP/FN on a sample). Output: `results/rescore_published.tsv` and `results/rescore_summary.txt`.

| Condition | Method | Published FN % | Our rescoring FN % |
|---|---|---|---|
| RNASim1000 | FastTree | 14.9 | 14.94 |
| | IQ-TREE | 15.1 | 15.09 |
| | RAxML-NG | 15.1 | 15.15 |
| | Constrained-INC / FT | 15.1 | 15.09 |
| | **GTM / FT** | 14.7 | **14.70** |
| | TreeMerge / FT | 14.7 | 14.70 |
| | Constrained-INC / IQ | 14.5 | 14.46 |
| | **GTM / IQ** | 14.4 | **14.36** |
| | TreeMerge / IQ | 14.4 | 14.36 |
| | TreeMerge-PAUP* (Table A2) | 14.7 | 14.70 |
| Cox1-HET | FastTree | 23.9 | 23.95 |
| | IQ-TREE | 19.6 | 19.58 |
| | RAxML-NG | 18.2 | 18.23 |
| | Constrained-INC / FT | 18.9 | 18.87 |
| | **GTM / FT** | 18.9 | **18.85** |
| | TreeMerge / FT | 18.9 | 18.89 |
| | Constrained-INC / IQ | 19.7 | 19.67 |
| | **GTM / IQ** | 18.7 | **18.71** |
| | TreeMerge / IQ | 18.7 | 18.72 |
| | TreeMerge-PAUP* (Table A2) | 18.9 | 18.92 |
| 1000M1-HF | FastTree | 50.9 | 51.43 ⚠ |
| | IQ-TREE | 30.2 | 30.18 |
| | RAxML-NG | 24.9 | 24.90 |
| | Constrained-INC / FT | 42.4 | 42.39 |
| | **GTM / FT** | 42.4 | **42.45** |
| | TreeMerge / FT | 42.5 | 42.55 |
| | Constrained-INC / IQ | 28.6 | 28.61 |
| | **GTM / IQ** | 28.4 | **28.37** |
| | TreeMerge / IQ | 28.5 | 28.51 |
| | TreeMerge-PAUP* (Table A2) | 42.5 | 42.55 |

⚠ The deposit contains four different full-data FastTree trees for 1000M1-HF, one in each `CreateConstraintTrees/{FastTree,IQTree2}/{120,500}` directory. Their errors are 51.43, 50.58, 51.05 and 50.68%. The published 50.9% lies inside that range; we could not tell which run the table used.
All 29 other entries agree to within the table's rounding (≤0.05 points).

### 2.2 Rerunning GTM ourselves

`code/rerun_gtm.py` runs github.com/vlasmirnov/GTM (HEAD 18e3bc9, default `convex` mode) on the published guide trees and the published IQ-TREE subset trees. Output: `results/gtm_rerun.tsv`.

- **39/40** replicate×guide cases give *exactly* the published GTM tree (RF = 0 to the published tree).
- The exception is 1000M1-HF R0 with the FastTree guide: RF 45 to the published tree, FN 43.07% (ours) vs 42.47% (published). The likely cause is that this run used a different FastTree guide tree, the same ambiguity as the ⚠ row above.
- GTM's `old` mode reproduces only 19/40, so the paper used the current (`convex`) mode.
- GTM takes 0.3–0.6 s per merge.

So the baseline is fully reproducible: same inputs, same code, same trees, same numbers.

### 2.3 Statistics on the published trees (paired across replicates)

Source: `results/published_paired.txt`. Wilcoxon signed-rank, two-sided.

- **No published condition has GTM worse than TreeMerge or Constrained-INC.**
- Pooled over all 40 replicate×guide pairs, TreeMerge is worse than GTM by +0.04 points (4 better / 12 worse / 24 tied, p = 0.046), and Constrained-INC by +0.33 points (p = 0.031).
- The "condition where GTM is less accurate than TreeMerge" is therefore not present in the published data.

> **Statistical caveat for the project:** with 5 replicates (RNASim1000, 1000M1-HF), the smallest attainable two-sided exact Wilcoxon p-value is 2/2⁵ = 0.0625. **No 5-replicate condition can ever reach p < 0.05 on its own.** Use ≥10 replicates per condition (Cox1-HET has 10; SATé 1000M1 has 20 replicates of which only 5 have the HF analysis), or pool conditions in a pre-registered way.

## 3. Pilot: a blending DTM

### 3.1 How much can *any* merger improve on GTM? (ceiling analysis)

**Recoverable true branches.** A true bipartition can appear in a merged tree T* with T*|S_i = T_i only if its restriction to every subset is trivial or a bipartition of T_i (`code/ceiling.py`). This gives a lower bound ("floor") on the FN rate of *any* disjoint tree merger of the given subset trees, blended or not. It is optimistic: the bound ignores that the recoverable branches must also fit together.

**Best unblended merge.** `code/oracle_gtm.py` runs GTM with the *true* tree as the guide. GTM returns the unblended merger closest to its guide, so this is the best unblended merge.

**Achievable blended merge.** `code/oracle_pub.py` runs two oracle-scored constrained searches: constrained insertion and constrained SPR, both scored by RF distance to the true tree. These are *heuristics* given the true tree, **not upper bounds**: GTM-Blend-ML with the ML score did better than them on one replicate (1000M1-HF R2, FastTree guide: 40.6% vs 41.2%).

| Condition / guide | GTM | best unblended | oracle-guided blended heuristic | floor (any DTM, optimistic) | recoverable missing branches per rep |
|---|---|---|---|---|---|
| RNASim1000 / FT | 14.70 | 14.56 | – | 14.50 | 2.0 |
| RNASim1000 / IQ | 14.36 | 14.28 | – | 14.28 | 0.8 |
| Cox1-HET / FT | 18.85 | 18.65 | – | 18.36 | 11.5 |
| Cox1-HET / IQ | 18.71 | 18.60 | – | 18.34 | 8.6 |
| 1000M1-HF / FT | 42.45 | 42.27 | **41.06** | 30.28 | 120.8 |
| 1000M1-HF / IQ | 28.37 | 28.29 | **27.74** | 26.57 | 17.8 |

Interpretation:
- **On RNASim1000 and Cox1-HET, the subset trees determine the error.** No merger can improve on GTM by more than 0.1–0.5 points.
- On 1000M1-HF with the FastTree guide, the per-branch floor (30.3%) looks like 12 points of headroom. But oracle-guided blended searches reach only 41.1% (−1.4 points). Most of that floor is probably not jointly realizable with these subset trees; we have no exact optimum, since the blended problem is NP-hard.
- **A better DTM, blended or not, cannot show a meaningful accuracy gain on the published small conditions.** That is a strong negative result in itself, and it is exactly the kind of finding the course values.

### 3.2 The method: GTM-Blend-ML (constrained SPR from GTM, ML-scored)

Code: `code/blend.py` (search and feasibility), `code/mlspr.py` (likelihood).

**Feasibility.** An SPR move prunes subtree X at node p and regrafts it on a target edge. Let e₁…e_m be the path edges from p, with far sides F₁ ⊇ … ⊇ F_m.
- The move removes the bipartitions F₁…F_{m−1} and adds F₂∪X … F_m∪X.
- Because each T_i is binary, the new tree still has T'|S_i = T_i iff every restriction (F_j∪X)|S_i, for j = 2..m, is trivial or already a split of T_i.
- This is an O(m·k) bitset test per move.
- We verified it against explicit restriction on random moves (0/60 mismatches). Every final tree is also re-checked: all trees reported here keep every subset tree induced.

**Score.**
- Every feasible move within radius 4 is scored by a lazy-SPR log-likelihood.
- The likelihood is GTR+Γ4, using Felsenstein pruning with directional partial likelihoods. Only the pendant branch and the split of the target branch are re-fitted, on a 3×3 grid; model parameters come from RAxML-NG.
- The top 8 distinct topologies are re-scored by RAxML-NG `--evaluate` with full branch-length optimization.
- The best move is accepted if it raises logL; the search repeats until no move improves.
- The lazy score tracks the full score closely: it is typically within 1–10 logL units of RAxML-NG and is a lower bound.

**Why this works.** Feasible moves are rare: about 400–1,400 per round, versus about 50,000 unconstrained moves at the same radius. So exhaustive ML scoring of the feasible neighbourhood is affordable.

**Reachability.** With an oracle score (RF to the true tree), constrained SPR from GTM reaches the **true tree (0% FN) in 1–13 moves** on every simulated case with exact subset trees that we tested. The search space is therefore adequate; the question is whether ML ranks the moves correctly.

### 3.3 Negative results (kept on purpose)

1. **Parsimony is the wrong criterion.** Constrained parsimony SPR from GTM (`run_blend.py`) and parsimony constrained insertion + SPR (`run_insert.py`) lower the parsimony score in every case, but they *increase* FN. Source: `results/pilot_published.tsv`.

   | Method | Condition | n | GTM FN % | method FN % | mean diff (95% CI) | better/worse/tie | Wilcoxon p |
   |---|---|---|---|---|---|---|---|
   | Parsimony constrained SPR (r = 6) | Cox1-HET (FT + IQ) | 20 | 18.78 | 19.06 | +0.29 | 1/18/1 | <0.004 each guide |
   | | RNASim1000 (FT + IQ) | 10 | 14.53 | 14.94 | +0.40 | 0/10/0 | 0.0625 each guide (n = 5 floor) |
   | | 1000M1-HF FT | 5 | 42.45 | 42.14 | −0.30 (−0.75, +0.16) | 3/1/1 | 0.375 |
   | | 1000M1-HF IQ | 5 | 28.37 | 28.41 | +0.04 | 1/2/2 | 0.5 |
   | | **pooled** | 40 | 21.88 | 22.08 | **+0.21 (+0.10, +0.30)** | 5/31/4 | **1.7e-4** |
   | Parsimony insertion + SPR | Cox1-HET | 20 | 18.78 | 19.61 | +0.83 | 0/20/0 | 0.002 each guide |
   | | RNASim1000 | 10 | 14.53 | 15.06 | +0.52 | 0/10/0 | 0.0625 each guide |
   | | 1000M1-HF FT | 5 | 42.45 | 42.00 | −0.44 (−0.72, −0.18) | 4/0/1 | 0.125 |
   | | 1000M1-HF IQ | 5 | 28.37 | 28.77 | +0.40 | 1/4/0 | 0.125 |
   | | **pooled** | 40 | 21.88 | 22.41 | **+0.54 (+0.37, +0.70)** | 5/34/1 | **1.0e-5** |

   The one exception is 1000M1-HF with the poor FastTree guide, where parsimony gives a small, non-significant gain; that is the only published condition where the merge, rather than the subset trees, limits accuracy (§3.1).

   In simulation the true tree had a *worse* parsimony score than GTM's tree (17,170 vs 17,167 on a 200-taxon caterpillar), and a parsimony move increased FN from 1.5% to 5.1%.

2. **Greedy insertion-based blending** (Constrained-INC style, start from the largest subset tree and insert taxa one by one on feasible edges) is unreliable:
   - with parsimony placement it is ≈ GTM on 1000M1-HF;
   - with ML placement (EPA-ng, `code/mlinsert.py`) it was much better or much worse than GTM depending on the replicate. On a caterpillar case it gave 54.8% vs GTM's 12.7%.

   Starting from GTM's tree and searching locally is much more robust than building a blended tree from scratch.

3. **Unconstrained-looking gains on the published data are small.** On 1000M1-HF R1 (FastTree guide), GTM-Blend-ML raised logL by 406 units in 12 moves, yet FN did not change. Better likelihood under the constraints does not imply better topology when the subset trees themselves are 30% wrong.

### 3.4 Simulation: a condition where GTM must lose

**Design** (`code/sim.py`, `code/sim_ml.py`, `code/sim_batch.sh`).
- 200-taxon model trees: Yule or caterpillar, with internal branch lengths ~Exp(mean 0.005 or 0.002) and pendant ~Exp(mean 0.1).
- AliSim GTR+Γ4 alignments, 1000 sites.
- A FastTree guide tree.
- Centroid decomposition of the guide into subsets of ≤50 taxa (k = 4–7).
- Subset trees are either *exact* (true tree restricted, which isolates merger error) or estimated by IQ-TREE 2 on each subset (realistic).
- Mergers compared: GTM (guide = FastTree tree), GTM with the true tree as guide (best unblended), oracle insertion (shows that blended trees with 0% merger error exist), and GTM-Blend-ML.

All paired across the same replicates; results in `results/sim_results.tsv` and `results/sim_summary.txt`.

Extra arms:
- **TreeMerge**: the original github.com/ekmolloy/treemerge with PAUP* 4a169 for pairwise branch lengths. Distances are topological (edge-count) distances on the FastTree guide, as in Park et al. Runner: `code/sim_treemerge.py`; our only change to TreeMerge is `list(graph.neighbors(...))` for networkx ≥ 2.
- **IQ-TREE 3 on the full alignment** (GTR+G, Yule replicates only, for lack of time on the caterpillar ones).

Mean FN % over 20 replicates (each condition: 200 taxa, 1000 sites, max subset 50, k = 4–7 subsets):

| Condition | guide (FastTree) | subset trees | GTM | best unblended (GTM, true guide) | TreeMerge | **GTM-Blend-ML** | IQ-TREE full data |
|---|---|---|---|---|---|---|---|
| A. Yule, internal mean 0.005, exact subset trees | 57.6 | 0.0 | 28.02 | 24.8 | 28.05 | **13.68** | 33.05 |
| B. Yule, internal mean 0.005, IQ-TREE subset trees | 57.6 | 32.7 | 50.20 | 48.5 | 50.58 | **44.11** | 33.05 |
| C. Caterpillar, internal mean 0.002, exact subset trees | 76.3 | 0.0 | 26.57 | 21.1 | 24.95 | **16.52** | (1 rep only) |

Paired tests (two-sided Wilcoxon signed-rank; difference = B − A, in FN points):

| Condition | Comparison | mean diff (95% bootstrap CI) | better/worse/tie | Wilcoxon p | sign-test p |
|---|---|---|---|---|---|
| A | GTM-Blend-ML vs GTM | **−14.34 (−17.41, −11.57)** | 20/0/0 | **8.8e-5** | 1.9e-6 |
| A | TreeMerge vs GTM | +0.03 (−0.61, +0.74) | 7/6/7 | 0.92 | 1 |
| A | GTM-Blend-ML vs TreeMerge | −14.37 | 20/0/0 | 8.7e-5 | 1.9e-6 |
| B | GTM-Blend-ML vs GTM | **−6.09 (−7.64, −4.72)** | 20/0/0 | **8.8e-5** | 1.9e-6 |
| B | TreeMerge vs GTM | +0.38 (−0.20, +1.07) | 8/7/5 | 0.41 | 1 |
| B | GTM-Blend-ML vs TreeMerge | −6.47 | 20/0/0 | 8.7e-5 | 1.9e-6 |
| B | GTM-Blend-ML vs IQ-TREE (full data) | **+11.07 (+9.09, +13.05)** (worse) | 0/20/0 | 8.8e-5 | 1.9e-6 |
| C | GTM-Blend-ML vs GTM | **−10.05 (−15.10, −5.89)** | 20/0/0 | **8.8e-5** | 1.9e-6 |
| C | **TreeMerge vs GTM** | **−1.62 (−3.71, +0.33)** | 11/5/4 | **0.083** | 0.21 |
| C | GTM-Blend-ML vs TreeMerge | −8.43 | 20/0/0 | 8.8e-5 | 1.9e-6 |

What this shows:
- **When the subsets are not clades of the true tree, unblended merging is the bottleneck,** even with perfect subset trees: GTM 27–28% error, best unblended 21–25%, best blended 0%.
- GTM-Blend-ML removes half or more of that merger error in every replicate.
- In all 60 replicates every subset tree is still exactly induced, and logL rises (mean +150 to +280 units).
- The ML search does not reach the true tree. The true tree's logL is still 47–507 units higher, so a larger radius or better moves could improve further.
- **Caterpillar trees are where TreeMerge comes closest to beating GTM** (p = 0.083 with 20 replicates). This is the instructor's suggested condition. A 4-week project should run more replicates and caterpillar variants there.
- **Honest limitation:** at this scale, full-data ML (IQ-TREE, 33.1%) is much more accurate than any DTM pipeline with *estimated* 50-taxon subset trees (GTM 50.2%, GTM-Blend-ML 44.1%). Small subsets on a hard, short-branch tree give bad subset trees. The simulations demonstrate the merger effect, not a practical pipeline win. A project should test larger n and larger subsets, where full ML becomes expensive and the DTM regime is realistic.
- Condition A's IQ-TREE comparison (13.7% vs 33.1%) is not a fair comparison, because GTM-Blend-ML received exact subset trees there.

### 3.5 GTM-Blend-ML on the published 1000M1-HF data

This uses the published guide trees, IQ-TREE subset trees and alignments. The search starts from the published GTM tree (radius 4, top-8 RAxML-NG confirmation, at most 30 rounds). R2 with the FastTree guide was started before the cap was added and was stopped at round 65. Single thread. Source: `results/pilot_published.tsv`.

| Guide | Replicate | GTM FN % | GTM-Blend-ML FN % | moves | ΔlogL | minutes |
|---|---|---|---|---|---|---|
| FastTree | R0 | 42.47 | **40.86** | 26 | +686 | 87 |
| | R1 | 38.55 | 38.55 | 12 | +406 | 27 |
| | R2 | 43.00 | **40.58** | 63 | +1306 | ~250 (stopped) |
| | R3 | 38.97 | 38.87 | 17 | +459 | 55 |
| | R4 | 49.24 | **48.23** | 30 (cap) | +725 | 119 |
| | **mean** | **42.45** | **41.42** | | | |
| IQ-TREE | R0 | 26.20 | 26.10 | 1 | +115 | 9 |
| | R1–R4 | 28.91 (mean of 25.33, 28.30, 28.20, 33.80) | 28.91 (unchanged) | 0–4 | +0–51 | 3–13 |
| | **mean** | **28.37** | **28.35** | | | |

Paired:
- FastTree guide: −1.03 points (95% CI −1.81, −0.26), 4/0/1, Wilcoxon p = 0.125.
- IQ-TREE guide: −0.02, 1/0/4, p = 1.
- Pooled 10: −0.52 (−1.08, −0.05), 5/0/5, p = 0.0625. That is the smallest p attainable with 5 non-zero pairs, so this is **not significant at α = 0.05 by construction**; more replicates are required.

The gain matches the bounded headroom of §3.1: about 1 point where the guide is poor, nothing where it is good. Every output keeps all subset trees induced and has higher likelihood than GTM's tree. We did not run Cox1-HET or RNASim1000: §3.1 shows ≤0.5 points of headroom there, and RNASim1000's 21,946-site alignment makes each RAxML-NG evaluation slow.

## 4. Answers to the instructor's prompts

- *"GTM does not allow blending… develop a better DTM that allows blending."*
  - GTM-Blend-ML is a blending DTM: it outputs a tree that induces every subset tree exactly, with subsets interleaved where the data support it.
  - Its search space is provably the full set of constraint-preserving SPR neighbours.
  - It is a post-processor on top of GTM, so it never needs a worse starting point than GTM.
- *"Find a condition where GTM is less accurate than TreeMerge."* and *"Test GTM and TreeMerge under balanced vs. caterpillar model trees."*
  - **Published data:** none. GTM is never worse than TreeMerge; pooled, TreeMerge is +0.04 points worse (p = 0.046) (§2.3).
  - **Simulation (§3.4):**
    - On Yule (balanced-ish) trees GTM and TreeMerge tie (p = 0.92 and 0.41).
    - On **caterpillar** trees with a poor guide (FastTree 76% error), TreeMerge is better than GTM by 1.6 points on average (11 better / 5 worse / 4 tied, Wilcoxon p = 0.083, n = 20).
    - That is a trend in the predicted direction, not yet a significant finding. A project should replicate it with more replicates and caterpillar variants.
  - The mechanism is explicit: decomposition subsets that are not clades of the true tree force any unblended merger to err. TreeMerge can blend, but its blending is driven by guide-tree distances, so it recovers little. An ML-scored blending search recovers much more (GTM-Blend-ML beats TreeMerge 20/20 in every condition).
- *"Test GTM in a pipeline with methods other than maximum likelihood (e.g., NJ or FastME)."*
  - Covered by the sibling pilot on branch `claude/cs581-forest` (`cs581/forest/REPORT.md`). We did not redo it.
  - That pilot ran GTM with NJ/FastME guide and subset trees (centroid decomposition, and Forest components). It found that centroid decomposition + GTM(NJ) beats NJ by 8–12 FN points on deep trees, and that GTM(FastME) ties FastME.
  - Our parsimony results (§3.3) add a related point. A non-ML criterion used *inside the merger* (to decide blending) is harmful here: parsimony picks trees with lower parsimony but higher error.

## 5. Is this a viable 4-week project?

**Yes, with the right framing.** The hypothesis worth testing: *a likelihood-scored, constraint-preserving SPR search starting from GTM corrects GTM's unblended merging errors when the decomposition subsets are not clades of the true tree, and costs nothing when they are.*

### Plan

| Week | Work | Deliverable |
|---|---|---|
| 1 | Reproduce (done here: validation table, GTM rerun; TreeMerge+PAUP* already runs on new inputs). Install Constrained-INC (github.com/steven-le-thien/INC). Speed up GTM-Blend-ML (accept several non-overlapping improving moves per round). | Table 2.1, all baselines runnable |
| 2 | Simulation grid: guide-tree quality (internal branch mean 0.002/0.005/0.01), tree shape (Yule/caterpillar), subset trees exact vs IQ-TREE, n = 200 and 1000, max subset 50/100/500. 20 replicates each. GTM vs TreeMerge vs CINC vs GTM-Blend-ML. | Main figure with Wilcoxon tests |
| 3 | Published data: GTM-Blend-ML on 1000M1-HF (10 cases), Cox1-HET (20 cases), RNASim1000 (10 cases), reporting no-harm and bounded-headroom results against the ceiling analysis. Optionally RNASim10k (10 reps) with the FastTree guide, where GTM's own paper used it at scale. | Table vs published numbers |
| 4 | Runtime engineering (accept several non-overlapping moves per round; cap rounds), ablations (radius, top-K, parsimony vs ML), write-up. | Report |

### Data (all public)
- doi:10.13012/B2IDB-7008049_V1: Park et al. 2021 inputs/outputs (RNASim1000, Cox1-HET, 1000M1-HF analyses, RNASim10k/50k).
- doi:10.13012/B2IDB-5139418_V1: SATé 1000M1 true trees (needed for 1000M1-HF), 20 replicates.
- Simulated data: generated by `code/sim.py` (AliSim from IQ-TREE 3; deterministic seeds).

### Success criteria (pre-register them)
- Primary: FN rate, paired against GTM on the same replicates, two-sided Wilcoxon signed-rank, α = 0.05, ≥10 replicates per condition (20 preferred).
- Report the mean difference with a bootstrap 95% CI, and better/worse/tie counts.
- Secondary: ML score (RAxML-NG, fixed model) and running time.
- Report the ceiling (§3.1) next to every published-data result so that "no improvement" is interpretable.

### Baselines
GTM (convex mode), best unblended merge (GTM with the true guide; oracle), TreeMerge, Constrained-INC, and the full-data ML methods (IQ-TREE 2/3, RAxML-NG, FastTree).

### Risks
1. **Headroom on published data is tiny (shown in §3.1).** Mitigation: make simulation the primary evidence and use published data for "does no harm".
2. **Runtime.** About 4–10 min per 200-taxon replicate; 3–120 min per 1000-taxon case (single thread, 30-round cap; one uncapped case ran about 4 h). Mitigation: batch non-conflicting moves, smaller radius, multi-threaded RAxML-NG.
3. **Estimated subset trees with high error limit gains,** because the constraints themselves are wrong. This is real and should be reported, not hidden.
4. **TreeMerge setup.** It runs here (Python 2.7 env + dendropy 4.3 + networkx, PAUP* 4a169 test build with libgfortran 4, one networkx-2 patch); see `code/sim_treemerge.py`. The PAUP* test build expires, so pin a copy.
5. **Novelty check before writing.**
   - Blending DTMs exist: NJMerge, TreeMerge and Constrained-INC all blend.
   - ML local search under *multiple partial* constraint trees is, as far as we found, not in RAxML-NG or IQ-TREE, which each take a single constraint tree.
   - Re-check the 2025–2026 literature (e.g. "constrained SPR" with "multiple constraint trees").
6. **Full-data ML wins at small n (§3.4).** Any project claim must be about the DTM regime (n ≥ 1000, where full ML is slow or fails), or explicitly about merger error.

### Recommendation
Do it, framed as above, and lead with the simulations. It is low-risk because the hard parts already work:
- validation, a feasibility test with proof sketch, a working ML search, and a simulation pipeline;
- a strongly significant effect in the condition where theory says GTM must fail.

Being honest about the bounded headroom on the published benchmarks is itself a finding worth reporting.

## 6. Reproducing

```bash
bash cs581/code/setup.sh                       # FastTree, IQ-TREE 3, MAFFT, Python deps
MAMBA_ROOT_PREFIX=/opt/mm/root /opt/mm/micromamba install -y -n bio -c conda-forge -c bioconda epa-ng raxml-ng
bash cs581/gtm/code/fetch_data.sh               # Data Bank + SATé 1000M1 + GTM source -> /opt/gtmdata
cd cs581/gtm/code
python3 validate.py /opt/gtmdata ../results/rescore_published.tsv
python3 rerun_gtm.py /opt/gtmdata/GTM_src /opt/gtmdata/rerun ../results/gtm_rerun.tsv
python3 ceiling.py ../results/ceiling.tsv; python3 oracle_gtm.py /opt/gtmdata/GTM_src /opt/gtmdata/rerun ../results/oracle_unblended.tsv
python3 oracle_pub.py 1000M1-HF ../results/oracle_blend_1000M1.tsv
python3 published_paired.py ../results/rescore_published.tsv
./pilot_all.sh /opt/gtmdata/pilot                         # parsimony pilot (negative)
python3 run_mlspr_pub.py 1000M1-HF R0 FT /opt/gtmdata/mlpub 4   # GTM-Blend-ML on published data
./sim_batch.sh yule 0.005 exact 1 20 3                    # simulation (also: yule 0.005 iqtree, cat 0.002 exact)
# TreeMerge arm: Python 2.7 env + PAUP* (test build) + original TreeMerge (one networkx-2 patch, see §3.4)
MAMBA_ROOT_PREFIX=/opt/mm/root /opt/mm/micromamba create -y -n tm27 -c conda-forge python=2.7 numpy networkx "libgfortran-ng=7"
/opt/mm/root/envs/tm27/bin/pip install dendropy==4.3.0
#   PAUP*: https://phylosolutions.com/paup-test/paup4a169_ubuntu64.gz -> /opt/tools/paup/paup
#   TreeMerge: git clone https://github.com/ekmolloy/treemerge; copy python/ to /opt/tools/treemerge/python;
#   sed -i -E 's/([A-Za-z_]+)\.neighbors\(([^)]*)\)/list(\1.neighbors(\2))/g' treemerge.py njmerge2.py
./tm_batch.sh
for r in $(seq 1 20); do python3 sim_fullml.py /opt/gtmdata/simq/yule_n200_m50_i0.005/$r; done   # full-data IQ-TREE baseline
python3 sim_summary.py /opt/gtmdata/simq ../results/sim_results.tsv
python3 pilot_summary.py /opt/gtmdata/pilot /opt/gtmdata/mlpub ../results/pilot_published.tsv
```
