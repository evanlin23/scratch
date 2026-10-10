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
