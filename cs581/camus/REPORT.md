# CAMUS pilot: better base tree, adaptive quartet filter, scalability

*Pilot for one CS581 (Fall 2026) project idea. About 4.5 h of wall-clock on 4 cores and 15 GB RAM. Branch `claude/cs581-camus`.*
*All numbers come from `results/SUMMARY.md`, `results/per_rep.csv` and `results/SUMMARY_gtrue.md`, which `code/aggregate.py` regenerates.*

## 1. Question

The networks lecture slide "CAMUS: Summary and Future Work" (`notes/lectures/CS581-intro-to-phylogenetic-networks.txt`) lists:
- "Improve scalability to more than 200 species";
- "Develop other techniques to estimate the underlying tree T";
- "Extend technique to estimate more complex networks";
- and notes "No guarantees of statistical consistency, even if the true network is level-1".

This item is open problem 11 in `literature/slide_open_problems.md` and entry K in `literature/wide_scan.md` §2.7.

CAMUS pipeline:
1. ASTRAL species tree T.
2. Quartets from the gene trees.
3. Quartet filter at threshold t.
4. Exact DP for the level-1 network containing T that maximizes the number of satisfied quartets.

I piloted three ideas:
- **(a)** A better base tree T: other estimators, a "swap" to the other displayed tree, choosing T by the CAMUS objective, and a true-tree oracle that measures the headroom.
- **(b)** A sample-size-aware (adaptive) quartet filter in place of the fixed ratio threshold t.
- **(c)** Profiling where CAMUS's time and memory go.

## 2. Prior art (details and DOIs in `prior_art.md`)

- **CAMUS paper.** Willson & Warnow, *Bioinformatics* 42(S1) btag245 (ISMB 2026), doi:10.1093/bioinformatics/btag245. It already lists "CAMUS inside a heuristic search that explores different constraint trees" as future work, and t = 0.5 was tuned on one 26-taxon condition. So idea (a) is *announced* by the authors' group, not open ground.
- **NetCS.** Dai & Molloy, WABI 2026, doi:10.4230/LIPIcs.WABI.2026.1. Level-1 reconstruction from quartets via majority vote. Near-perfect given the true tree of blobs, and it runs 200 taxa × 1000 genes in minutes. **This is a direct, fast competitor**, and its error also comes from the tree step.
- **Per-quartet hypothesis tests.** NANUQ (doi:10.1186/s13015-019-0159-2) and TINNiK/NANUQ+ (doi:10.1186/s13015-024-00266-2) already classify 4-taxon sets with per-quartet tests at levels α and β. Idea (b) is therefore new *for CAMUS* but conceptually anticipated.
- **Other methods.** SNaQ / SNaQ.jl (doi:10.1371/journal.pgen.1005896; bioRxiv doi:10.1101/2025.11.17.688917), PhyloNet-MPL (doi:10.1186/1471-2164-16-S10-S10, unverified), Squirrel (doi:10.1093/molbev/msaf067), PhyNEST (doi:10.1093/sysbio/syae054), InPhyNet (bioRxiv doi:10.1101/2025.05.05.652278). I found no CAMUS follow-ups.

## 3. Data, metric and baseline reproduction

**Data.** Illinois Data Bank doi:10.13012/B2IDB-6892704.
- **V2** (IDB-0062483, 2026-08-23) adds `inferred-networks.tar.xz` (217 KB). It holds the **published CAMUS, ASTRAL, SNaQ and PhyloNet networks** for 20 replicates per condition, plus a t sweep (f0, f02, f05, f08) on n25 FastTree. V1 has no CAMUS outputs.
- I streamed the 4.3 GB archive and extracted only the trees and networks (0.9 GB, no sequences), in `/opt/data/camus`.

**Conditions used** (n = ingroup taxa, plus OUT; 1000 gene trees each):

| condition | gene trees | replicates | role |
|---|---|---|---|
| n25 FastTree | g_500.nwk | 20–39 | **training** (the condition the paper tuned t on) |
| n15 IQ-TREE | iqtree_500.nwk | 00–19 | held out |
| n25 IQ-TREE | iqtree_500.nwk | 00–19 | held out |
| n50 FastTree | g_500.nwk | 00–19 (all variants on 00–11; pre-registered subset on 12–19) | held out |
| n25 true gene trees | g_true.nwk | 20–39 | control |

**Metric.**
- `code/netmetric.py` computes softwired clusters as the union of the clusters of all 2^r displayed trees, excluding trivial clusters.
- FN = |C(true) \ C(est)| / |C(true)| and FP = |C(est) \ C(true)| / |C(est)|.
- It matches **PhyloNet `CmpNets -m cluster` to all digits on 5/5 test pairs**. That is the metric the paper used; PhyloNet took more than 5 h per network at 51 taxa, while ours takes milliseconds.
- Error = (FN+FP)/2 of the **1-reticulation network**, which is the paper's Exp. 3 protocol ("exactly one reticulation"). I also report the network with the true number of reticulations; conclusions are the same.

**Reproduction.**
- CAMUS was built from github.com/jsdoublel/camus @ fb6a287 with defaults: `-t 0.5 -q 2`, base tree = the published ASTRAL tree rooted at OUT.

| condition | k=1 network identical to published | published FN / FP | ours FN / FP |
|---|---|---|---|
| n15 IQ-TREE | 20/20 | 0.159 / 0.112 | 0.159 / 0.112 |
| n25 FastTree | 20/20 | 0.110 / 0.082 | 0.110 / 0.082 |
| n50 FastTree | NN50/NN50 | P50 | O50 |
| n25 IQ-TREE | 10/20 | 0.131 / 0.116 | 0.137 / 0.110 |

- The n25 IQ-TREE mismatch is in the input, not CAMUS. The published `camus_f05_astral_iqtree` runs used as T the ASTRAL tree built from the **FastTree** gene trees (2/2 replicates checked), not the shipped `astral-iqtree.nwk`. I kept the consistent pipeline (ASTRAL on the same gene trees). CAMUS releases v1.0.0–v1.0.2 give identical results.
- The paper reports no numeric FN/FP values, only figures, so "reproduce a published number" means exact network identity plus the same mean FN/FP from the published networks.

## 4. Methods tested

Code is in `code/`. `run_rep.py` runs every variant on one replicate and is restartable; `aggregate.py` does the stats.

- **(b) Adaptive filter.** `camus-ztest.patch` adds environment-variable switches to CAMUS's `Threshold.Keep`.
  - Counts are sorted c1 ≥ c2 ≥ c3. Default CAMUS keeps the second topology iff `floor(t(c2+c3)) < c2 − c3`. This is a ratio rule that ignores how many gene trees support the minor topologies, so with c2 = 3, c3 = 0 a noise pattern passes.
  - Variants: `zXonly` keeps it iff (c2−c3−1)/√(c2+c3) > X, a continuity-corrected sign test of H0: c2 = c3 (equal minor topologies under ILS alone). `zXand` and `zXor` combine the test with the t = 0.5 rule.
  - Also run: a t sweep (0, 0.2, 0.3, 0.5, 0.7, 0.8) and a per-replicate oracle t.
- **(a) Base tree.**
  - Estimators: ASTRAL (default), TREE-QMC (`tqmc`) and wASTRAL (`wastral`), all from ASTER/bioconda. ASTRID did not install from bioconda and was skipped.
  - `swap`: the other displayed tree of default's 1-reticulation network.
  - Score selection: pick, among {default, tqmc, wastral, swap}, the 1-reticulation network with the largest CAMUS objective (filtered quartets displayed). This is my own numpy reimplementation, `code/quartets.py`.
  - `true_major` (oracle): the true network's major displayed tree.
- Variants for the held-out conditions were fixed after the training condition, in `results/PREREGISTRATION.md`: `z3only` (primary), `z5only`, `tqmc`, and score selection.

## 5. Results

Paired by replicate. diff = variant − default (negative = better). W/T/L uses a tie band of |diff| < 1e-9; errors are discrete fractions, so ties are exact. p = two-sided Wilcoxon signed-rank with zeros dropped. Full tables are in `results/SUMMARY.md`.

### 5.1 Base tree (idea a)

RESULTS_A

### 5.2 Adaptive filter (idea b)

RESULTS_B

### 5.3 True gene trees (control, n25, reps 20–39)

RESULTS_C

### 5.4 Runtime and memory (idea c)

RESULTS_D

## 6. Verdict

VERDICT
