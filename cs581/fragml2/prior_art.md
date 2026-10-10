# Prior art: placement + graft + ML polish on fragmentary alignments

Compiled 2026-10-10 for the CS581 pilot. Citations marked [v] were checked online during this search. Citations marked [m] come from memory, so check the DOI before citing them.

## 1. Smirnov & Warnow (the direct baseline)

Smirnov V, Warnow T. *Phylogeny estimation given sequence length heterogeneity.* Syst Biol 70(2):268–282 (2021; online 2020). doi:10.1093/sysbio/syaa058 [v]. PMC7875441.

- **Data.** ROSE 1000M1–M4 (from SATe; 1000 taxa, about 1000 nt, 20 replicates), RNASim and RNASim2 (1000 taxa), and the biological sets 16S.M (901 taxa) and 23S.M (278 taxa).
  - **LF:** 25% of sequences become fragments with mean length 50% of the median length.
  - **HF:** 50% of sequences become fragments with mean length 25% of the median length.
  - Fragment lengths are normal with SD 60 (as in Nguyen et al. 2015, UPP). Each fragment is a random substring.
- **Protocol 1: align then ML.** Alignment with PASTA, SEPP or UPP (backbone tree from FastTree or RAxML), then a tree from FastTree or RAxML-NG under GTRGAMMA. RAxML-NG used 1 start on the simulated data and 5 on the biological data.
- **Protocol 2: placement.** PASTA backbone alignment and tree on the full-length sequences. Fragments are added to the alignment by SEPP or UPP, then placed with pplacer, pplacer(c) (restricted to the SEPP-chosen subtree) or APPLES. Fragments that land on the same edge form a polytomy. **There is no ML polishing after placement.** The authors only say that placement trees "would need additional techniques to refine". EPA-ng is mentioned but not evaluated.
- **HF results (missing-branch rate, FN):**

| | 1000M1 | M2 | M3 | M4 | RNASim | 16S.M |
|---|---|---|---|---|---|---|
| UPP(R)-RAxML | **0.370** | 0.304 | 0.237 | 0.167 | 0.377 | 0.340 |
| UPP(R)-pplacer | 0.488 | 0.437 | 0.380 | 0.320 | 0.507 | 0.496 |
| PASTA-RAxML | 0.765 | 0.616 | 0.355 | 0.164 | 0.436 | 0.409 |

  LF on 1000M1: UPP-RAxML 0.157, UPP-pplacer 0.215, PASTA-RAxML 0.246.
- **Other findings:** pplacer beats APPLES, and unconstrained pplacer beats pplacer(c). FastTree is poor whenever fragments are present. Whether the backbone tree came from FastTree or RAxML barely matters.
- **Time (1000M2-HF):** placement pipelines take about 30–35 min, mostly alignment. RAxML-NG alone takes 40–60 min. UPP-RAxML takes 70–80 min in total.
- **Implication for the pilot:** on 1000M1-HF, placement alone leaves about 12 FN points on the table compared with UPP+RAxML. The open question is whether a short ML polish starting from the grafted tree recovers them.

**Data.** The paper cites Dryad `10.5061/dryad.8pk0p2nj8`. **That DOI returns 404 at doi.org and in the Dryad API.** The real deposit is **doi:10.5061/dryad.95x69p8h8** [v]. Its Dryad metadata links it to syaa058 and to a Zenodo supplement, 10.5281/zenodo.7688198. The licence is CC0.

| file | size | versioned download (current v6; v3 identical, same sha256) |
|---|---|---|
| README.txt | 1,596 B | https://datadryad.org/downloads/file_stream/2134136 (v3: .../399220) |
| UnalignFragTree.zip | 1,774,981,278 B | https://datadryad.org/downloads/file_stream/2134135 (v3: .../399219) |

- The API endpoints are `https://datadryad.org/api/v2/files/2134135/download` and `.../2134136/download`. The whole dataset is at `https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.95x69p8h8/download`.
- The zip's sha256 is `381bff069562aa5fcc88f92879dd12fbf9725d0ac3b40d736f349bca2e991042`.
- From this sandbox, `file_stream` returns 403 (the AWS ELB blocks it) and `api/v2/.../download` returns 401. Download it from a browser or a campus machine.
- The supplement PDF is downloadable: https://zenodo.org/api/records/7688198/files/Supplement.pdf/content.
- Another source for 1000M1(HF) analyses is the Park/Warnow DTM dataset on the Illinois Data Bank, IDB-7008049 (https://databank.illinois.edu/datasets/IDB-7008049). It is described as containing "analyses of RNASim1000, Cox1-Het, and 1000M1(HF)". It also returned 403 from here.

## 2. uDance

Balaban M, Jiang Y, Zhu Q, McDonald D, Knight R, Mirarab S. *Generation of accurate, expandable phylogenomic trees with uDance.* **Nat Biotechnol** 42(5):768–777 (2024). doi:10.1038/s41587-023-01868-8 [v]. It is not in Genome Research.

- **Pipeline:**
  1. Place the queries on the backbone with APPLES-2.
  2. Partition the placement tree by size and diversity, adding outgroups.
  3. Run an **unconstrained** ML gene tree per partition (RAxML-NG by default, or IQ-TREE 2 or RAxML-8).
  4. Build an ASTRAL-constrained species tree per partition. It must match the backbone on the outgroups; elsewhere it is free by default.
  5. Stitch the partitions together.

  This is "place then local re-estimation", applied to multi-gene species trees.
- **Fragment handling:**
  - Before inference it drops a gene from a partition if fewer than 100 sites remain.
  - It drops sequences shorter than 75 bp per partition per gene. The justification cited is Sayyari 2017.
  - It runs TreeShrink on FastTree gene trees.
  - Backbone QC uses leave-one-out APPLES-2 reinsertion.

  So fragments are **filtered out** rather than placed and kept. There is no fragmentary-data benchmark of the kind used by the Warnow lab.

## 3. Graft/placement-then-search practice

- **EPA.** Berger SA, Krompass D, Stamatakis A. Syst Biol 60(3):291–302 (2011). doi:10.1093/sysbio/syr010 [m]. **EPA-ng:** Barbera P et al. Syst Biol 68(2):365–369 (2019). doi:10.1093/sysbio/syy054 [m]. EPA-ng places onto a fixed tree and writes jplace for gappa [v].
- **pplacer.** Matsen FA, Kodner RB, Armbrust EV. BMC Bioinf 11:538 (2010). doi:10.1186/1471-2105-11-538 [m]. Includes `guppy tog`, which grafts placements onto the tree.
- **gappa.** Czech L, Barbera P, Stamatakis A. Bioinformatics 36(10):3263 (2020). doi:10.1093/bioinformatics/btaa070 [m]. `gappa examine graft` makes the grafted tree.
- **PhyloMagnet.** Schön ME et al. Bioinformatics 36(6):1718 (2020) [m]. Uses graft for taxonomic screening, not tree inference.
- **SCAMPP.** Wedell E, Cai Y, Warnow T. IEEE/ACM TCBB (2023) [v]. **BSCAMPP:** Wedell, Shen, Warnow, WABI 2023 (doi:10.4230/LIPIcs.WABI.2023.3) and TCBB 2025 [v]. These scale pplacer/EPA-ng to trees of 200k leaves, using subtrees of 2000 by default for BSCAMPP+EPA-ng. The largest gains are on fragmentary queries. These papers measure placement error, not grafted-tree-then-ML.
- **SCAMPP+FastTree.** Chu & Warnow, Bioinf Adv 3:vbad008 (2023) [v].
- **APPLES.** Balaban M et al. Syst Biol 69(3):566 (2020) [m].
- **INC / Constrained-INC.** Zhang Q, Rao S, Warnow T. *Constrained incremental tree building…* Algorithms Mol Biol 14:2 (2019). doi:10.1186/s13015-019-0136-9 [v].
  - Taxa are inserted by quartet voting, guided by constraint trees.
  - Le T et al. (AlCoB 2019; TCBB 2021) found raw INC inaccurate. INC-ML, which uses ML subset trees, did better.
  - These insert full-length taxa, not fragments.
- **GTM / DTM.** Park M, Zaharias P, Warnow T. *Disjoint tree mergers for large-scale ML tree estimation.* Algorithms 14(5):148 (2021) [v].
  - Includes **1000M1-HF**.
  - The pipeline is a FastTree starting tree, then a decomposition, then RAxML-NG/IQ-TREE on the subsets, then GTM.
  - FastTree 2 is less accurate than RAxML on alignments with many fragments, and GTM can match or beat FastTree and IQ-TREE.
- **Review.** Zaharias P, Warnow T. *Recent progress on methods for estimating and updating large phylogenies* (2021/22 preprint) [v]. Says placement is "also useful when the input exhibits sequence length heterogeneity". No graft+polish experiment.
- **PUmPER.** Izquierdo-Carrasco F, Cazes J, Smith SA, Stamatakis A. Bioinformatics 30(10):1476 (2014) [v]. Extends old trees with new taxa (RAxML-Light), then continues the ML search. This is the closest "extend then ML polish" tool in the RAxML lineage. Not evaluated on fragments.
- **UShER + matOptimize.** UShER: Turakhia Y et al. Nat Genet 53:809 (2021) [m]. matOptimize: Ye C et al. Bioinformatics 38(15):3734 (2022) [m]. Greedy MP placement followed by SPR optimisation, explicitly because "stepwise addition… sometimes places new sequences suboptimally" [v].
  - **Thornlow B et al.**, *Online phylogenetics with matOptimize produces equivalent trees and is dramatically more efficient… than de novo and ML implementations* (bioRxiv 10.1101/2021.12.02.471004; PMID 35611334) [v]. This is the strongest quantitative "place then optimise ≈ full search, much cheaper" result. It is on dense SARS-CoV-2 data with near-complete genomes, not fragments.
- **TIPars.** Ye et al. (2022) [m]. Parsimony-based insertion into a reference tree.
- **Constraints.**
  - RAxML-NG `--tree-constraint` accepts multifurcating or **incomplete** constraints. Taxa not in the constraint are "free to move" [v; Krumlov 2024 slides; RAxML Google group]. A resolved backbone used as a constraint therefore fixes the induced backbone topology while the fragments are placed by SPR. That is the variant in the pilot.
  - IQ-TREE `-g` behaves similarly [m].
  - RAxML-NG `--tree <file>` / IQ-TREE `-t` accept the grafted tree as a starting tree.

## 4. Evidence that constrained or placement-initialised ML helps on fragments

**I found none that targets fragments.** The nearest evidence:

- **uDance:** divide-and-conquer RAxML-NG beats FastTree-2 on the full data.
- **Thornlow et al.:** online placement plus SPR equals de novo ML at a fraction of the cost, on non-fragmentary viral data.
- **Park et al. 2021:** a FastTree-initialised D&C ML pipeline on 1000M1-HF.
- **Starting-tree literature:** RAxML-NG speedups depend on the start type, being higher with parsimony starts (e.g. *The Free Lunch is not over yet*, PMC10518076 [v]). A good starting tree mainly saves time.

## 5. Known fragment issues in ML (brief)

- **Sayyari E, Whitfield JB, Mirarab S.** MBE 34(12):3279–3291 (2017). doi:10.1093/molbev/msx261 [v/m]. Fragmentary sequences reduce gene-tree and species-tree accuracy; filtering them helps.
- **Wiens JJ.** Syst Biol 52(4):528 (2003) [m]. Very incomplete taxa can be placed accurately if enough characters exist overall. **Wiens & Morrill**, Syst Biol 60(5):719 (2011) [m].
- **Sanderson MJ, McMahon MM, Steel M.** BMC Evol Biol 10:155 (2010) [m]. Terraces and incomplete coverage.
- **Simmons MP.** Cladistics 28:208 (2012) [m]. Missing data can mislead likelihood analyses through branch-length effects.
- **TreeShrink.** Mai U, Mirarab S. BMC Genomics 19:272 (2018). doi:10.1186/s12864-018-4620-2 [m]. Removes long-branch outliers, which are often fragments.
- **Rogue taxa.** RogueNaRok, Aberer, Krompass & Stamatakis, Syst Biol 62:162 (2013) [m]. Fragments often behave as rogues.

## What is new and what is not

**Known practice:**

- **Backbone tree on full-length sequences plus placement of fragments.** Smirnov & Warnow 2021, SEPP/UPP, SCAMPP/BSCAMPP. On 1000M1-HF it has a measured FN of about 0.49, against about 0.37 for UPP+RAxML-NG.
- **Grafting placements into a tree.** `guppy tog`, `gappa graft`.
- **Placement or insertion followed by tree rearrangement.** Done in UShER+matOptimize (parsimony, SARS-CoV-2), PUmPER (RAxML-Light), uDance (local unconstrained RAxML-NG plus ASTRAL) and INC-ML/GTM (D&C ML).
- **Incomplete constraint trees in RAxML-NG and IQ-TREE.** Standard features.
- **Fragment-induced ML error.** Sayyari 2017; FastTree in Smirnov 2021 and Park 2021.

**Apparently untested (no paper found):**

1. A head-to-head comparison on the Warnow-lab fragmentary benchmarks (1000M1–M4 HF/LF, RNASim, 16S.M) of four approaches, measuring FN/FP, ML score and wall/CPU time:
   - (a) placement-graft alone;
   - (b) graft followed by a *short unconstrained* RAxML-NG/IQ-TREE polish from the grafted tree;
   - (c) RAxML-NG with the backbone as a non-comprehensive `--tree-constraint`;
   - (d) a full RAxML-NG search on the same UPP alignment.

   Smirnov & Warnow stop at (a) versus (d), and none of the D&C papers use placement to seed the search.
2. Whether polishing closes the roughly 12-point FN gap between (a) and (d) on 1000M1-HF, and how many SPR rounds that takes.
3. Whether the backbone constraint in (c) hurts when the backbone itself is wrong. Smirnov found the backbone method barely matters for placement, but a constraint freezes backbone errors.
4. The **EPA-ng premasking/rate-scaler bug** found by the previous pilot, which affects fragment placement on trees with more than 2000 tips. No literature mentions it. Its downstream effect on a graft+polish pipeline is unstudied: does polishing repair the misplacements, or do they cost ML search time? The default BSCAMPP subtree size (2000) may hide it. This is plausibly the most novel angle.
5. Cost accounting that charges the alignment step (UPP) to every arm. In Smirnov & Warnow, alignment dominated the time of the placement pipelines.
