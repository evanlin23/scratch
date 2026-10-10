# Pilot: better quartet amalgamation under ILS + HGT (CS581, Fall 2026)

Status: pilot of about 4 hours, run 2026-10-09. Nothing is decided. Every number below
comes from `results/tables.md`, which `code/analyze.py` regenerates from the per-replicate
`result.json` files.

## 1. Questions and where they come from

| Question | Source |
|---|---|
| "Can we design better quartet amalgamation methods?" | Phylogenomics part 2, slide 35, "(Some) Open Questions" (`notes/lectures/CS581-phylogenomics-2026-part2.txt`) |
| "Address gene tree heterogeneity due to multiple causes (including horizontal gene transfer)" | First-day slide "Open problems (and possible course projects)" (`notes/lectures/CS581-Fall2026-firstday.txt`) |
| Background: the dominant quartet is the species-tree quartet under ILS, GDL and bounded random HGT (Allman et al.; Legried et al.; Roch & Snir; Daskalakis & Roch), so "Finding MQSST (within a constrained search space) [is] a powerful approach that maintains consistency" and is "fairly robust" | Phylogenomics part 2, slide 34 |

The pilot tests three candidate ideas:
- **(a)** SPR hill-climbing on the explicit quartet score, starting from each heuristic's tree.
- **(b)** Harvesting clusters from several heuristics and solving ASTRAL's DP over the enlarged cluster set.
- **(c)** Downweighting quartets or genes that are likely affected by HGT.

## 2. Prior art

The full notes, with DOIs and verified/memory flags, are in `results/prior_art.md`.

- **TREE-QMC**: Han & Molloy, *Genome Res* 2023, doi:10.1101/gr.277629.122. **Weighted TREE-QMC**: *Syst Biol* 2025, doi:10.1093/sysbio/syaf009. No post-hoc local search, no HGT weighting.
- **wQFM-TREE**: Rafi et al., *Bioinf Adv* 2025, doi:10.1093/bioadv/vbaf053. FM-style moves happen only inside each bipartition step.
- **wASTRAL**: Zhang & Mirarab, *MBE* 2022, doi:10.1093/molbev/msac215. Weights model gene-tree estimation error (support, branch length), not HGT.
- **ASTER / ASTRAL-IV**: Zhang, Nielsen & Mirarab, *MBE* 2025, doi:10.1093/molbev/msaf172. CASTER: *Science* 2025.
  - Each round does random-order placement, then **NNI refinement**, then a **DP over tripartitions harvested from all rounds**.
  - It runs 4 rounds plus subsampling rounds until there is no improvement. `-g` injects user guide trees.
  - So ASTER already contains a version of (a) and a version of (b).
- **Q-SPR**: Arasti & Mirarab, *Algorithms Mol Biol* 2024, doi:10.1186/s13015-024-00257-3.
  - Optimal-SPR hill-climbing on the quartet score, from ASTRAL-III and ASTER starts.
  - The score improved in a minority of cases, by less than 0.5%. RF did not improve consistently (45 better, 84 worse).
  - **(a) is essentially published.**
- **FASTRAL**: Dibaeinia, Tabe-Bordbar & Warnow, *Bioinformatics* 2021, doi:10.1093/bioinformatics/btab093. Its X comes from ASTRID runs on gene subsets, a single-method harvest that *shrinks* X.
- **ASTRAL-III** `-e`/`-f` and **ASTRAL-MP** (randomized X expansion) already provide the machinery for (b).
- **HGT**: Davidson et al. 2015, doi:10.1186/1471-2164-16-S10-S1; Roch & Snir 2013; Daskalakis & Roch 2016.
  - The MSC-deviation statistic exists in MSCquartets and is used in NANUQ (doi:10.1186/s13015-019-0159-2) and TINNiK (doi:10.1186/s13015-024-00266-2), but only for network/blob inference.
  - **No paper found uses an HGT or MSC-deviation signal to downweight the input to a quartet amalgamation search.** (c) looks open.

The search took 30 minutes, so 2025–26 preprints may be missed.

## 3. Data, tools, and baseline reproduction

### Data

**HGT+ILS data (Davidson et al. 2015), doi:10.13012/B2IDB-6670066_V1.**
- 51 taxa: 50 ingroup taxa plus an outgroup. 6 HGT rates: 0, 0.08, 0.2, 0.8, 8 and 20 expected HGT events per gene.
- 50 replicates × 1000 genes; we used replicates 01–20.
- Gene trees are FastTree-2 estimates ("est", which carry SH-like support for wASTRAL) or the true gene trees.
- *k*-gene conditions use the first *k* genes. We ran 50 and 200 genes for all rates. We ran 1000 genes for the two highest rates: replicates 01–10 at 8/gene and 01–20 at 20/gene.
- Error is computed over all 51 leaves, unrooted.

**ILS-only data with published species trees (Nute et al. 2018), doi:10.13012/B2IDB-7735354_V1.**
- "full" (no missing data); 26 taxa (25 ingroup taxa plus an outgroup); RAxML gene trees; 20 replicates.
- Used for (i) validation, and (ii) a second, harder testbed: the three highest-ILS conditions (500K-1E-6, 500K-1E-7, 2M-1E-7) with 50 and 200 genes.

### Tools

| tool | version |
|---|---|
| ASTRAL-IV (`astral4`) and wASTRAL (`wastral`, hybrid weighting) | ASTER v1.25.4.8 and v1.25.3.8 (bioconda) |
| ASTRAL-III | 5.7.8 (bioconda `astral-tree`) |
| TREE-QMC | 3.0.4 (bioconda) |
| wQFM-TREE | GitHub `abdur-rafi/wQFM-TREE` @7d7c24b, using its own `run.sh` (PAUP* greedy consensus + jar) |
| ASTRID | Our own re-implementation: mean topological internode distance + FastME 2.1.6.3 (balanced-ME taxon addition + SPR/NNI). ASTRID-2 is not on bioconda. |

### Our code
- `code/qtool.cpp`: builds the explicit quartet table Q[q][3] from gene trees, then does scoring, SPR hill-climbing and per-gene support. It uses all C(n,4) quartets: 249,900 for n = 51.
- `code/run_rep.py`: runs every method on one replicate. It is restartable, writing `result.json` per method.
- `code/analyze.py`: produces the tables.
- `code/validate_published.py`: the validation below.

### Validation
- **`qtool`'s score matches ASTRAL's printed score exactly.** We rescored all 360 published ASTRAL trees (ASTRAL 4.10.5; 18 conditions × 20 replicates). Our normalized quartet score equals the "Normalized quartet score" in each published log up to the log's printed precision: max |diff| = 5.0e-9 over 360 trees. This also confirms that the published *k*-gene runs used the first *k* RAxML gene trees.
- **Published accuracy, as mean FN rate over 20 replicates** (`results/validation_published.tsv`):

| condition | genes | ASTRAL (RAxML gene trees) | ASTRID | MP-EST | ASTRAL (true gene trees) | SVDquartets |
|---|---|---|---|---|---|---|
| 500K-1E-6 | 1000 | 0.161 | 0.150 | 0.176 | 0.061 | 0.211 |
| 500K-1E-6 | 50 | 0.354 | 0.374 | 0.457 | 0.252 | 0.513 |
| 2M-1E-6 | 1000 | 0.028 | 0.022 | 0.033 | 0.013 | 0.072 |
| 10M-1E-6 | 1000 | 0.011 | 0.011 | 0.015 | 0.000 | 0.039 |

  The ordering matches the paper's qualitative findings: SVDquartets is the worst summary method, ASTRAL and ASTRID are close, and MP-EST is behind. All 18 conditions are in the TSV.
- **Our ASTRAL-IV rerun vs the published ASTRAL 4.10.5 trees** (same RAxML gene trees; 120 replicate×condition pairs in the ILS testbed):
  - ASTRAL-IV's normalized quartet score is ≥ the published one in 120/120 pairs, and strictly higher in 25.
  - Yet its mean FN rate is **+0.5 pp higher**.
  - So the reproduction is consistent, and it is a first sign that a higher quartet score does not mean a more accurate tree.
- **Davidson et al. 2015**: we did not try to match their figures numerically. They report ASTRAL-2 and wQMC; we run ASTRAL-IV. Our ASTRAL-IV error rises from about 2% (no HGT, 200 genes) to about 5–12% at 20 HGT/gene. That is the same qualitative trend as their figure ("ASTRAL robust until the highest HGT rate").

## 4. Methods tested

All candidates share the MQSST objective: maximize Σ_q Q[q][T|q] over the gene-tree quartet table.

**(a) `ls[X]`: SPR hill-climbing from X ∈ {ASTRAL-IV, TREE-QMC, wQFM-TREE, ASTRID}.**
- First-improvement, random order; each candidate is rescored exactly against the quartet table.
- Regraft radius is ≤ 4 edges, a superset of NNI. Full-radius SPR found the same optima in all 5 spot checks and is about 4× slower.
- The search stops when a full neighborhood has no strict improvement.

**(b) Cluster harvesting:**
- `harvest`: ASTRAL-IV with `-g` guide trees from TREE-QMC, wQFM-TREE, ASTRID and wASTRAL.
- `astral3+e`: ASTRAL-III's exact DP over its default X plus all bipartitions of those trees (`-e`; polytomies resolved arbitrarily). It is compared against plain `astral3`.

**(c) HGT-aware objectives** (decided before seeing results; no tuning):
- `capminor`: per quartet, cap the larger minority count at the smaller one. Under the MSC the two minority topologies have equal expected frequency, so any directional excess is attributed to HGT or error and removed.
- `vote`: each quartet votes only for its dominant topology, weight 1. This is a pure "dominant quartet" consensus, the quantity that the bounded-HGT theorems are about.
  - Both `capminor` and `vote` optimize the modified objective by SPR local search from the ASTRAL-IV tree.
- `reweight`: per-gene weight = (fraction of the gene's quartets agreeing with the ASTRAL-IV tree / max)^4. It is implemented as 0–4 copies of each gene tree, and ASTRAL-IV is rerun on the weighted multiset. This is a robust-regression style downweighting of outlier (HGT-like) genes.

**Protocol.**
- Replicates 01–05 are development replicates; we inspected them while debugging.
- Replicates 06–20 (06–10 at 1000 genes for 8/gene) are held out, and every reported test uses only them.
- Paired by replicate. Error = normalized RF.
- Two-sided Wilcoxon signed-rank on per-replicate error differences, with zero differences dropped.
- W/T/L uses a tie band of exactly equal RF; errors are discrete, so any nonzero difference is ≥ 1 edge.
- Runtime is single-threaded wall-clock on a 4-core VM with 4 jobs in parallel.

## 5. Results

All tables use the held-out replicates only. Full per-condition tables (mean error, mean score, per-condition tests, dev-set numbers) are in `results/tables.md`. Δerr is the candidate's normalized RF minus the baseline's, in percentage points; negative means the candidate is better. "est" = FastTree gene trees, "true" = true gene trees, "raxml" = the ILS testbed.

### 5.1 Headroom: ASTRAL-IV already finds the best quartet score anyone finds

"Best" means the highest normalized quartet score among all 15 trees per replicate: 5 baselines, ASTRAL-III, ASTRAL-III+e, 4 local searches, harvest and the 3 HGT variants.

| data | genes | n | ASTRAL-IV is best | max gap (×100) | err ASTRAL-IV % | err of per-replicate best method (oracle) % |
|---|---|---|---|---|---|---|
| HGT | 50 est | 90 | 89/90 | 0.002 | 7.38 | 5.49 |
| HGT | 200 est | 90 | 90/90 | 0 | 4.00 | 2.55 |
| HGT | 1000 est | 20 | 20/20 | 0 | 1.46 | 1.04 |
| HGT | 50 true | 90 | 90/90 | 0 | 4.35 | 3.36 |
| ILS | 50 raxml | 45 | 43/45 | 0.037 | 28.21 | 22.42 |
| ILS | 200 raxml | 45 | 45/45 | 0 | 15.27 | 11.50 |

**ASTRAL-IV attains the best score found in 377 of 380 held-out replicates.** The only gaps are on the 26-taxon, 50-gene, very high ILS data, and they are tiny.

There is accuracy headroom: the oracle column picks the best method per replicate after the fact. But no method reaches it by scoring higher.

### 5.2 Is the MQSST objective aligned with accuracy? (`results/alignment_summary.tsv`)

We scored the **true species tree** and compared it with the ASTRAL-IV tree:

| data | genes | n | true tree scores lower | equal | higher |
|---|---|---|---|---|---|
| HGT | 50 est | 90 | 85 | 5 | 0 |
| HGT | 200 est | 90 | 82 | 8 | 0 |
| HGT | 1000 est | 20 | 12 | 8 | 0 |
| HGT | 50 true | 90 | 73 | 17 | 0 |
| ILS | 50 raxml | 45 | 45 | 0 | 0 |
| ILS | 200 raxml | 45 | 42 | 3 | 0 |

**The true tree never scores higher than ASTRAL-IV's tree, under any condition.** It does not even with true gene trees. Whenever ASTRAL-IV is wrong, the error comes from the objective at finite sample size, not from the search.

The same holds under the two HGT-aware objectives. Under `capminor` the true tree is higher in only 10 of 380 held-out replicates, and under `vote` in only 9 of 380 (per-condition counts in the alignment file).

### 5.3 Idea (a): SPR local search on the quartet score

| data | genes | candidate vs baseline | n | Δerr (pp) | W/T/L | p | score ↑/↓ |
|---|---|---|---|---|---|---|---|
| HGT | 50 est | ls[ASTRAL-IV] vs ASTRAL-IV | 90 | 0.00 | 0/90/0 | n/a | 0/0 |
| HGT | 200 est | ls[ASTRAL-IV] vs ASTRAL-IV | 90 | 0.00 | 0/90/0 | n/a | 0/0 |
| ILS | 50 raxml | ls[ASTRAL-IV] vs ASTRAL-IV | 45 | 0.00 | 0/45/0 | n/a | 1/0 |
| HGT | 200 est | ls[TREE-QMC] vs TREE-QMC | 90 | +0.32 | 11/58/21 | 0.12 | 40/0 |
| HGT | 50 est | ls[TREE-QMC] vs TREE-QMC | 90 | −0.09 | 21/52/17 | 0.42 | 49/0 |
| ILS | 50 raxml | ls[ASTRID] vs ASTRID | 45 | +1.74 | 12/12/21 | 0.19 | 44/0 |
| ILS | 200 raxml | ls[wQFM-TREE] vs wQFM-TREE | 45 | −0.77 | 17/15/13 | 0.15 | 38/0 |

- From an ASTRAL-IV start, SPR moves the tree in 1 of 380 held-out replicates.
- From TREE-QMC, wQFM-TREE or ASTRID starts, SPR always raises the score (up to the ASTRAL-IV level). Accuracy changes in both directions and is never significantly better.
- This matches Q-SPR (Arasti & Mirarab 2024).

### 5.4 Idea (b): harvesting clusters from several heuristics

| data | genes | candidate vs baseline | n | Δerr (pp) | W/T/L | p | score ↑/↓ |
|---|---|---|---|---|---|---|---|
| HGT | 50 / 200 / 1000 est | harvest (ASTRAL-IV `-g`) vs ASTRAL-IV | 90 / 90 / 20 | +0.02 / 0 / 0 | 0/89/1, 0/90/0, 0/20/0 | n/a | 0/0 in all three |
| HGT | 50 est | ASTRAL-III+e vs ASTRAL-III | 90 | 0.00 | 3/85/2 | 1.0 | 7/0 |
| ILS | 50 raxml | ASTRAL-III+e vs ASTRAL-III | 45 | −0.10 | 4/38/3 | 1.0 | 11/0 |
| ILS | 50 raxml | ASTRAL-III+e vs ASTRAL-IV | 45 | −1.06 | 6/37/2 | 0.14 | 0/18 |

- Enlarging X helps ASTRAL-III's score a little: higher in 11 of 45 ILS 50-gene replicates. It never beats ASTRAL-IV's score. ASTER's own randomized rounds already cover the clusters the other heuristics contribute.
- Side observation: on the hardest ILS data, ASTRAL-III(+e) has a *lower* score than ASTRAL-IV but slightly *lower* error (−1.0 pp, p = 0.14). That is consistent with 5.2.

### 5.5 Idea (c): HGT-aware objectives

Pooled over the 6 HGT rates:

| data | genes | candidate vs ASTRAL-IV | n | Δerr (pp) | W/T/L | p |
|---|---|---|---|---|---|---|
| HGT | 50 est | capminor | 90 | +0.02 | 8/71/11 | 0.81 |
| HGT | 200 est | capminor | 90 | −0.21 | 11/75/4 | 0.19 |
| HGT | 1000 est | capminor | 20 | 0.00 | 1/18/1 | 1.0 |
| HGT | 50 true | capminor | 90 | 0.00 | 12/66/12 | 0.46 |
| HGT | 50 est | vote | 90 | +0.28 | 4/75/11 | 0.12 |
| HGT | 200 est | vote | 90 | −0.09 | 7/78/5 | 0.58 |
| HGT | 50 est | reweight | 90 | +0.46 | 19/44/27 | 0.15 |
| HGT | 200 est | reweight | 90 | +0.16 | 13/58/19 | 0.097 |
| ILS | 50 raxml | reweight | 45 | **+2.80** | 9/17/19 | **0.007** |
| ILS | 200 raxml | reweight | 45 | **+1.64** | 6/22/17 | **0.009** |

Highest HGT rate (20 events/gene) only:

| genes | capminor Δerr, W/T/L, p | vote | reweight |
|---|---|---|---|
| 50 est | −0.28, 3/9/3, 0.88 | +0.56, 2/8/5, 0.39 | +1.11, 4/5/6, 0.50 |
| 200 est | **−0.83, 4/11/0, 0.12** | −0.42, 2/13/0, 0.50 | +0.14, 2/8/5, 0.81 |
| 1000 est | +0.14, 0/14/1, 1.0 | +0.14, 0/14/1, 1.0 | +0.56, 1/10/4, 0.31 |
| 50 true | +0.14, 4/6/5, 0.50 | +0.28, 2/10/3, 0.44 | +0.69, 3/7/5, 0.41 |

- No HGT-aware variant is significantly better anywhere.
- The one hint, `capminor` at 20 HGT/gene with 200 genes (4 wins, 0 losses), does not replicate with 50 or 1000 genes, or with true gene trees. capminor alone has 20 per-condition tests, so one 4/0 split is expected by chance.
- Gene reweighting is significantly *worse* on ILS-only data, as expected: deep coalescence looks like an outlier gene.

**Context: the HGT in this simulation barely hurts ASTRAL.** At 1000 genes ASTRAL-IV's error is 1.4% even at 20 HGT/gene; at 200 genes it is 2.2% (no HGT) vs 5.6% (20/gene). So the room for an HGT fix on these data is a few percentage points at most.

Other baselines, for reference:
- wASTRAL vs ASTRAL-IV on 50 FastTree genes: −0.63 pp, 26/50/14, p = 0.051.
- At 20 HGT/gene and 200 genes, TREE-QMC has 4.4% error vs ASTRAL-IV's 5.6%, at a lower quartet score.

## 6. Runtime

Mean seconds per replicate, single-threaded, with 4 jobs sharing 4 cores:

| data | genes | ASTRAL-IV | ASTRAL-III | wASTRAL | TREE-QMC | wQFM-TREE | ASTRID | ls[·] (SPR) | harvest | capminor / vote | reweight |
|---|---|---|---|---|---|---|---|---|---|---|---|
| HGT, 51 taxa | 50 | 1.2 | 2.0 | 1.6 | 1.7 | 1.8 | 0.2 | 2.9–3.7 | 1.3 | 3.1 | 2.7 |
| HGT, 51 taxa | 200 | 4.7 | 4.6 | 6.6 | 2.8 | 3.8 | 0.5 | 3.0–3.6 | 5.0 | 3.3 | 13.3 |
| HGT, 51 taxa | 1000 | 51.2 | 91.3 | 91.0 | 9.4 | 20.4 | 2.7 | 4.3–4.7 | 52.2 | 4.5 | 58.7 |
| ILS, 26 taxa | 50 | 0.4 | 3.6 | – | 0.6 | 1.3 | 0.1 | 0.1–0.2 | 0.4 | 0.2 | 0.4 |

- `qtool` keeps the full C(n,4) table: 249,900 quartets for n = 51. A table build takes about 1.3 s at 1000 genes, and each SPR neighbor costs about 1.3 ms.
- This does not scale beyond about 100 taxa without incremental scoring, as in Q-SPR or ASTRAL's tree-based weights.
- Each local-search time above includes building the table.

## 7. Verdict

**(a) SPR local search: not promising.**
- It is already published (Q-SPR, 2024), and ASTER already does NNI rounds.
- In 380 held-out replicates, SPR changed the ASTRAL-IV tree once.
- From weaker starts it raises the score but not the accuracy.

**(b) Multi-method cluster harvesting: not promising.**
- It is mechanically available (`-e`, `-g`).
- It never improved ASTRAL-IV's score or accuracy (0/200 score improvements on HGT data).
- It helps ASTRAL-III's score slightly, which is moot because ASTRAL-IV already does better.

**(c) HGT-aware down-weighting: unclear, leaning not promising on these data.**
- It is the only one of the three that is novel.
- Three simple a-priori variants give no significant gain. Gene reweighting hurts under ILS.
- More fundamentally, ASTRAL's error at 20 HGT/gene and 1000 genes is only 1.4%.
- The true tree is out-scored by the MQSST tree under all three objectives, so reweighting quartet counts, in the ways tried, does not fix the finite-sample problem.

**What the pilot does show.** For "Can we design better quartet amalgamation methods?" the bottleneck at these scales (26–51 taxa) is **not the optimizer**. ASTRAL-IV hits the best known score in 377/380 cases, and the true tree never scores higher. Better methods would have to change the *estimator*, for example:
- shrinkage or regularization of the quartet table;
- using quartet-level uncertainty (wASTRAL-like weights were the only baseline near significance, p = 0.051);
- a different statistic.

That is a statistics question, not the "better search" question the candidate ideas assumed.

### If this were the 4-week project anyway

The most defensible version is a reframed (c): **"When does HGT break quartet amalgamation, and can MSC-deviation weighting fix it?"**

- **Week 1:** simulate harder HGT with SimPhy: highway or recurrent transfers between fixed lineages, and non-uniform receptor choice. Bounded random HGT (as in Davidson) is provably fine for quartet methods, so it is the wrong regime. Reuse `qtool`, `run_rep.py`, `analyze.py` and the validated baselines.
- **Week 2:** implement a principled per-quartet weight in place of the ad-hoc cap. For example, use the MSCquartets T3 or star-test p-value or a likelihood ratio of "MSC vs MSC+directional excess", and optimize it with ASTRAL's DP rather than SPR. That needs either ASTER source changes or an explicit-quartet DP for 50 or fewer taxa.
- **Week 3:** held-out comparison against ASTRAL-IV, wASTRAL, TREE-QMC and wQFM-TREE on both highway-HGT and plain-ILS data. A no-harm check under ILS is mandatory, because `reweight` failed it.
- **Week 4:** write-up. The negative results in 5.1 and 5.2 are already a solid section.

### Risks
1. Under bounded random HGT there is little to gain. ASTRAL error ≤ 1.5% at 1000 genes, and the theory says quartet majority is consistent.
2. Any reweighting that helps under HGT can hurt under ILS. Gene reweighting cost +1.6 to +2.8 pp here.
3. A highway-HGT simulator has to be built and justified, which takes time from the method.
4. An explicit quartet table is O(n⁴) and limits experiments to about 100 taxa unless the weights are folded into ASTER's tripartition scoring.
5. A good result may still look incremental next to wASTRAL and weighted TREE-QMC.

### Limitations of this pilot
- One HGT simulator and 26–51 taxa.
- `capminor` and `vote` were optimized by local search from the ASTRAL-IV tree rather than by DP.
- ASTRID is our own re-implementation.
- Replicates 01–05 were used for debugging; all claims use 06–20.
- 1000-gene runs cover only the 8/gene and 20/gene rates.

## 8. Reproduce

```bash
# tools: ASTER, TREE-QMC, ASTRAL-III, FastME via micromamba/bioconda (see top of code/run_rep.py for paths);
# wQFM-TREE from GitHub (abdur-rafi/wQFM-TREE)
g++ -O3 -march=native -o /opt/runs/qtool cs581/quartets/code/qtool.cpp
python3 cs581/quartets/code/validate_published.py           # needs Nute et al. trees under /opt/data/miss
python3 cs581/quartets/code/run_rep.py <model> <rep> <ngenes> est|true /opt/runs/hgt   # one replicate
python3 cs581/quartets/code/run_rep.py --files <genes> <true_tree> <ngenes> <outdir>   # any dataset
python3 cs581/quartets/code/analyze.py && python3 cs581/quartets/code/alignment.py
```
