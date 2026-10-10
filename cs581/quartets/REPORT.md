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
