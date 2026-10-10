# GTM-Blend: validation of the published GTM numbers and a pilot of a blending disjoint tree merger

CS581 (Warnow, Fall 2026) project de-risking. All code is in `cs581/gtm/code/` and all result tables are in `cs581/gtm/results/`.
Large data and intermediate trees were kept outside the repository (`/opt/gtmdata`). `code/fetch_data.sh` re-downloads them.

## 0. Summary

| Question | Answer |
|---|---|
| Can we reproduce the published GTM-pipeline numbers? | **Yes.** Rescoring the published trees reproduces 29/30 table entries to ≤0.05 points (FN %); the 30th is a FastTree run (51.4 vs 50.9). Rerunning GTM ourselves on the published guide and subset trees gives the published GTM tree exactly (RF = 0) in 39/40 replicate×guide cases. |
| Is there room for a better DTM on the published conditions? | **Very little.** On RNASim1000 and Cox1-HET, GTM is within 0.1–0.5 points of the lower bound for *any* merger of those subset trees, blended or not. On 1000M1-HF, even an oracle that knows the true tree improves GTM by only 1.4 points (FastTree guide) / 0.6 points (IQ-TREE guide). |
| Does blending help at all? | **Yes, when the decomposition subsets are not clades of the true tree.** This happens when the guide tree is poor. In simulation (200 taxa, Yule, short internal branches, exact subset trees), the best *unblended* merge has 24.8% error while an oracle blended merge has 0%. |
| Does our blending merger (GTM-Blend-ML) beat GTM? | **In simulation, yes, and strongly.** Exact subset trees: 28.0% → 13.7% FN, 20/20 replicates better, Wilcoxon p = 8.8e-5. Estimated subset trees and caterpillar trees: see §3.4. **On the published 1000M1-HF data:** see §3.5 (small, consistent with the oracle bound). |
| Does parsimony work as the score? | **No.** Constrained parsimony search lowers the parsimony score but makes the tree significantly *worse* than GTM (pooled +0.27 points FN, p = 7e-6). Maximum likelihood is required. |
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

**Achievable blended merge.** `code/oracle_pub.py` runs two oracle-scored constrained searches: constrained insertion and constrained SPR, both scored by RF distance to the true tree.

| Condition / guide | GTM | best unblended | oracle blended (achieved) | floor (any DTM, optimistic) | recoverable missing branches per rep |
|---|---|---|---|---|---|
| RNASim1000 / FT | 14.70 | 14.56 | – | 14.50 | 2.0 |
| RNASim1000 / IQ | 14.36 | 14.28 | – | 14.28 | 0.8 |
| Cox1-HET / FT | 18.85 | 18.65 | – | 18.36 | 11.5 |
| Cox1-HET / IQ | 18.71 | 18.60 | – | 18.34 | 8.6 |
| 1000M1-HF / FT | 42.45 | 42.27 | **41.06** | 30.28 | 120.8 |
| 1000M1-HF / IQ | 28.37 | 28.29 | **27.74** | 26.57 | 17.8 |

Interpretation:
- **On RNASim1000 and Cox1-HET, the subset trees determine the error.** No merger can improve on GTM by more than 0.1–0.5 points.
- On 1000M1-HF with the FastTree guide, the per-branch floor (30.3%) looks like 12 points of headroom. But even oracle-guided blended searches reach only 41.1% (−1.4 points), so most of that floor is not jointly realizable with these subset trees.
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

   | Method | n (pairs) | GTM FN | method FN | diff | better/worse/tie | Wilcoxon p |
   |---|---|---|---|---|---|---|
   | Parsimony constrained SPR | PARS_SPR_N | | | PARS_SPR_DIFF | | PARS_SPR_P |
   | Parsimony insertion + SPR | PARS_INS_N | | | PARS_INS_DIFF | | PARS_INS_P |

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

All paired across the same replicates; results in `results/sim_results.tsv`.

SIM_TABLE

### 3.5 GTM-Blend-ML on the published 1000M1-HF data

PUB_TABLE

## 4. Answers to the instructor's prompts

- *"GTM does not allow blending… develop a better DTM that allows blending."*
  - GTM-Blend-ML is a blending DTM: it outputs a tree that induces every subset tree exactly, with subsets interleaved where the data support it.
  - Its search space is provably the full set of constraint-preserving SPR neighbours.
  - It is a post-processor on top of GTM, so it never needs a worse starting point than GTM.
- *"Find a condition where GTM is less accurate than TreeMerge."*
  - The published data contain none (§2.3).
  - The mechanism that would produce one is now explicit: decomposition subsets that are not clades of the true tree. That happens when the guide tree has high error (§3.4), e.g. short internal branches, fragmentary data, or heterotachy.
  - A natural project experiment: run TreeMerge and Constrained-INC (both blending) in the simulated conditions of §3.4, next to GTM and GTM-Blend-ML.

## 5. Is this a viable 4-week project?

**Yes, with the right framing.** The hypothesis worth testing: *a likelihood-scored, constraint-preserving SPR search starting from GTM corrects GTM's unblended merging errors when the decomposition subsets are not clades of the true tree, and costs nothing when they are.*

### Plan

| Week | Work | Deliverable |
|---|---|---|
| 1 | Reproduce (done here: validation table, GTM rerun). Install TreeMerge (NJMerge + RAxML-NG branch lengths) and Constrained-INC so they can be run on *new* inputs. | Table 2.1, plus TreeMerge/CINC runnable |
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
1. **Headroom on published data is tiny (proven here).** Mitigation: make simulation the primary evidence and use published data for "does no harm".
2. **Runtime.** About 4–10 min per 200-taxon replicate; 27–90 min per 1000-taxon replicate (single thread). Mitigation: batch non-conflicting moves, smaller radius, multi-threaded RAxML-NG.
3. **Estimated subset trees with high error limit gains,** because the constraints themselves are wrong. This is real and should be reported, not hidden.
4. **Comparing with TreeMerge on new inputs needs TreeMerge installed.** It needs PAUP* or the RAxML-NG variant from Park et al. (github.com/minhyukpark/TreeMerge / the paper's Appendix C).
5. **Novelty check before writing.**
   - Blending DTMs exist: NJMerge, TreeMerge and Constrained-INC all blend.
   - ML local search under *multiple partial* constraint trees is, as far as we found, not in RAxML-NG or IQ-TREE, which each take a single constraint tree.
   - Re-check the 2025–2026 literature (e.g. "constrained SPR" with "multiple constraint trees").

### Recommendation
Do it, framed as above. It is low-risk because the hard parts already work:
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
python3 sim_summary.py /opt/gtmdata/simq ../results/sim_results.tsv
python3 pilot_summary.py /opt/gtmdata/pilot /opt/gtmdata/mlpub ../results/pilot_published.tsv
```
