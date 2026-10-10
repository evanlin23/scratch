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
| n50 FastTree | 20/20 | 0.109 / 0.081 | 0.109 / 0.081 |
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

Error of the 1-reticulation network (mean over replicates):

| variant | n15 IQ (20) | n25 FT, train (20) | n25 IQ (20) | n50 FT (20) | held-out pooled diff (60) | W/T/L | p |
|---|---|---|---|---|---|---|---|
| default (ASTRAL) | 0.135 | 0.096 | 0.124 | 0.095 | — | — | — |
| tqmc (TREE-QMC) | 0.133 | 0.095 | 0.111 | 0.094 | −0.0056 | 16/33/11 | 0.27 |
| wastral | 0.189 | 0.110 | 0.154 | 0.105 | **+0.0314** | 15/18/27 | 0.004 |
| swap | 0.143 | 0.098 | 0.119 | 0.090 | −0.0008 | 14/35/11 | 0.36 |
| **true_major (oracle)** | **0.035** | **0.055** | **0.054** | **0.052** | **−0.0710** | **49/6/5** | **1e-9** |

Score-based selection among {default, tqmc, wastral, swap}, using no truth:

| condition | default | selected | diff | W/T/L | p | best of the 4 (oracle) |
|---|---|---|---|---|---|---|
| n15 IQ | 0.135 | 0.134 | −0.0010 | 4/13/3 | 0.94 | 0.112 |
| n25 FT (train) | 0.096 | 0.098 | +0.0013 | 5/11/4 | 0.82 | 0.079 |
| n25 IQ | 0.124 | 0.115 | −0.0090 | 5/13/2 | 0.22 | 0.098 |
| n50 FT | 0.095 | 0.085 | −0.0105 | 11/4/5 | 0.036 | 0.072 |

Findings:
- **Headroom is large.** With the true major tree as T, error falls by about 60% (0.118 → 0.047 pooled), in every condition.
- **The ASTRAL tree is almost never a displayed tree of the true network.** That holds in 2/20 (n15: 4/20), 1/20, 2/20 and 0/20 replicates. The mean distance to the nearest displayed tree is 1.5–2.3 bipartitions. Even with 1000 **true** gene trees it holds in only 2/20. Because CAMUS can only add edges to T, every wrong bipartition of T is permanent.
- **Off-the-shelf estimators don't close the gap.** TREE-QMC is a wash (−0.006, p = 0.27) and wASTRAL is significantly worse. Swapping to the other displayed tree does nothing.
- **The CAMUS objective is not aligned with accuracy.** On n25 FastTree, the default network scores *higher* than the true network on the objective in 12/20 replicates, and in 12/20 even with true gene trees (`results/objective_*.tsv`).
  - So a search over T that maximizes the CAMUS score has no reason to find the true tree.
  - Selection among 4 candidates helps a little at n50 (−0.0105, p = 0.036) but not elsewhere, and recovers only a small part of the oracle-best-of-4 gap.

### 5.2 Adaptive filter (idea b)

| variant | n15 IQ | n25 FT (train) | n25 IQ | n50 FT | held-out pooled diff (n=60) | W/T/L | p |
|---|---|---|---|---|---|---|---|
| default t=0.5 | 0.135 | 0.096 | 0.124 | 0.095 | — | — | — |
| **z3only** (pre-registered primary) | 0.130 | 0.087 | 0.097 | 0.091 | −0.0119 | 26/21/13 | 0.10 |
| **z5only** (pre-registered secondary) | 0.125 | 0.088 | 0.095 | 0.093 | **−0.0136** | 21/33/6 | **0.013** |
| z3and (test AND t=0.5) | | | | | −0.0090 (n=52) | 6/46/0 | 0.031 |
| best fixed t (sweep 0–0.8) | 0.5 | 0.2/0.5 | 0.3 | 0.2–0.5 | t=0.3: −0.0009 | | |
| per-replicate oracle t | 0.112 | 0.081 | 0.094 | 0.089 | (upper bound of choosing t per replicate) | | |

The z-test variants are on all held-out replicates; the t sweep and `z3and` are on n15, n25 IQ and n50 reps 00–11 only.

Findings:
- **The adaptive filter gives a small, fairly consistent gain.** It is about 1.2–1.4 points of cluster error (≈10% relative) and mostly comes from n25 IQ-TREE (z5only −0.028, p = 0.007). At n50 the effect is ≈0 (−0.002, p = 0.44).
- **With the true number of reticulations** instead of k = 1, z5only gives −0.0163 (p = 0.003) and z3only −0.0149 (p = 0.025).
- **The z3 vs z5 ordering flipped** between training and held-out, so "pick z by cross-validation" is not yet shown to beat a fixed z.
- **Per-replicate t gains are small.** No fixed t beats 0.5 in pooled held-out data. Choosing t per replicate by oracle gains only 0.01–0.03, roughly what the z-test already captures.
- **With true gene trees** (control), z3only gives −0.020 (p = 0.014) and z5only −0.016 (p = 0.016). The gain is not an artefact of gene-tree error.

### 5.3 True gene trees (control, n25, reps 20–39)

| variant | error (k=1) | diff vs default | W/T/L | p |
|---|---|---|---|---|
| default | 0.103 | — | — | — |
| z3only | 0.083 | −0.020 | 9/10/1 | 0.014 |
| z5only | 0.087 | −0.016 | 7/13/0 | 0.016 |
| tqmc | 0.105 | +0.002 | 4/12/4 | 0.84 |
| true_major (oracle) | 0.056 | −0.047 | 15/3/2 | 0.0004 |

Default CAMUS is **no more accurate with true gene trees** (0.103) than with FastTree gene trees (0.096) on the same replicates. Gene-tree estimation error is not the bottleneck; the fixed base tree and the filter are.

### 5.4 Runtime and memory (idea c)

- **n15–n50 timing.** Single-threaded, one CAMUS run takes 1.2 s at n15, 7.8 s at n25 and **188 s at n50**. At n50, **184 of 185 s is quartet extraction** from the 1000 gene trees: CAMUS enumerates every quartet of every gene tree into hash maps. The DP takes less than 0.1 s, because only a few hundred non-tree quartets survive the t = 0.5 filter.
- **n100 profile** (1 thread, replicate 00, first k gene trees, from `getrusage` and CAMUS log timestamps):

| gene trees | quartet extraction | DP etc. | total | peak RSS | non-tree quartets kept |
|---|---|---|---|---|---|
| 100 | 608 s | 480 s | 1088 s | 0.98 GB | 492,221 |
| 200 | 1225 s | 262 s | 1487 s | 1.22 GB | 281,442 |

- **Extrapolation.** Extraction is linear in gene trees, about 6 s per tree at n101, which gives roughly 100 min of single-thread extraction for 1000 trees. That fits the paper's ~20 min on 32 cores. With few genes, many second topologies pass the filter, so the DP grows (8 min at 100 genes); the paper's 156 GB at n201 must come from the O(n⁴) quartet maps.
- **Dense counting is cheaper.** My numpy dense counting (`code/quartets.py`, one 4-set per row, four-point condition) takes 0.034 s per tree at n51 (≈5× faster than CAMUS) and 0.78 s per tree at n101 (≈8× faster), single-threaded and under load.
  - A dense uint16 [C(n,4) × 3] table needs 25 MB at n101, 0.4 GB at n201 and 6 GB at n401.
  - So a dense or streamed quartet counter (or counting only the 4-sets that can pass the filter) would very likely let CAMUS run past 200 taxa on a laptop.
  - That is engineering rather than new science, and NetCS already does 200 taxa in minutes.
- I did not run CAMUS at 150 or 200 taxa: at 1000 genes the expected time is several hours per run on 4 cores, and the paper reports 156 GB.

## 6. Verdict

**Overall: unclear, leaning promising for a narrower project. Not promising for the original "use a different tree estimator" idea.**

Evidence summary:
- Baseline reproduction is exact: 70/80 published networks are identical, and the 10 mismatches trace to a different input tree in the published runs. The metric equals PhyloNet's.
- **(a) Base tree.**
  - The oracle shows the base tree is *the* dominant error source: −0.071 pooled, 49/6/5, p = 1e-9, and ASTRAL's T is a displayed tree in under 10% of replicates.
  - Every practical fix I tried fails or is negligible: TREE-QMC, wASTRAL, swap, and selection by CAMUS score.
  - The CAMUS objective prefers the CAMUS network over the true network in about 60% of replicates, so score-driven tree search, the paper's own future-work idea, is unlikely to work as is.
  - This is an informative negative result rather than a method.
- **(b) Adaptive filter.** A one-line, sample-size-aware sign test gives a real but small gain: −0.014 pooled held-out, p = 0.013, W/T/L 21/33/6, about 10% relative, holding with true gene trees. It is cheap, novel for CAMUS, and anticipated by NANUQ/TINNiK.
- **(c) Scalability.** The bottleneck is quartet extraction, not the DP. A dense counter is a clear 5–8× speed and >100× memory win, but it is engineering, and NetCS competes there.

**A 4-week project that would be worth doing** (my recommendation):
1. **Week 1.** Make the filter a proper test: sign test, multiple-testing control (Bonferroni or BH over 4-sets), and z chosen by leave-one-condition-out validation. Run all six conditions. Also add the n100 condition using our fast counter.
2. **Week 2.** "Network-aware base tree." Diagnose *why* ASTRAL's T is not a displayed tree: which bipartitions are wrong, and whether they sit next to the reticulation. Then test a cheap repair: contract T's low-support or near-reticulation edges and let CAMUS (or a small enumeration) resolve them using a criterion that is aligned with accuracy, e.g. a likelihood under the NMSC or a held-out quartet score rather than the raw CAMUS count. The oracle shows ~0.07 of headroom.
3. **Week 3.** Dense quartet counter in Go (a patch to CAMUS) and profiling at n100–n200. Compare with NetCS (code in TREE-QMC) on the same data.
4. **Week 4.** Write-up.

**Risks:**
- The (b) effect is small and condition-dependent (≈0 at n50), and the z3/z5 ranking flipped between training and held-out.
- (a) may not yield a practical method: the objective misalignment is a structural obstacle, and the authors' group is already working on tree search.
- NetCS is a fast, accurate competitor that may make CAMUS-specific improvements less interesting.
- Only 20 replicates per condition, and an ILS/reticulation mix of a single simulation design.

**Limitations of this pilot:**
- n50 reps 12–19 ran only the pre-registered variants (CPU budget).
- ASTRID was not tested.
- 3–7 n15 replicates had been aggregated before the pre-registration was written.
- n100+ accuracy was not tested.
