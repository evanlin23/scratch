# CS581 project scan: "easy new method" opportunities outside MSA merging

*Prepared 2026-10-09 for the CS581 (Fall 2026, Warnow) course project. MAGUS-style alignment merging is deliberately excluded (see `sota_review.md`).*

**How this was checked.**
- Every journal DOI below was resolved against Crossref on 2026-10-09.
- Every Illinois Data Bank DOI, file list and size was read from the Data Bank's own JSON (`https://databank.illinois.edu/datasets/IDB-xxxxxxx.json`). This needs a browser User-Agent; plain clients get HTTP 403, which is the same workaround `code/setup.sh` uses.
- README files were downloaded and read for the TreeMerge, DTM-ML, DISCO, GDL-comparison, FastMulRFS, missing-data, DISCO+QR, CAMUS and TIPP3 datasets.
- Lecture content comes from the text extracted into `notes/lectures/*.txt`.
- Spring 2025 student project slides were downloaded from tandy.cs.illinois.edu and read, to see what has already been tried.
- **(unverified)** marks anything I could not confirm.

---

## TL;DR

| Rank | Opportunity (new method) | Swappable component | Key dataset (with published baseline outputs?) |
|---|---|---|---|
| **1** | **GTM-Blend:** a disjoint tree merger that allows blending (constrained local search started from GTM) | The merge step of DTM pipelines (GTM) | doi:10.13012/B2IDB-7008049_V1: constraint trees, guide trees, and GTM/TreeMerge/Constrained-INC/IQ-TREE/RAxML-ng outputs (**yes**) |
| 2 | **ASTRID-Pro:** an internode distance corrected for gene duplication and loss (GDL), using orthologous pairs and speciation nodes only | The distance matrix in ASTRID / NJst | doi:10.13012/B2IDB-5721322_V1 (FastMulRFS data with ASTRAL/STAG/DupTree/MulRF/FastMulRFS species trees; **yes**) |
| 3 | **DISCO-R:** root and tag gene-family trees by reconciliation with a first-pass species tree, then decompose | DISCO's root-and-tag step | doi:10.13012/B2IDB-4050038_V1 (DISCO data; baselines must be rerun, which takes minutes) |
| 4 | **Forest+DTM:** merge the "reliable forest" from the Daskalakis–Mossel–Roch fast-converging method with GTM | The final step of a distance pipeline | Simulator and code at github.com/lanl/distphylo (Kim et al. 2026); NJ/FastME baselines rerun in seconds |
| 5 | **CAMUS base tree / quartet filter:** a better underlying tree or quartet filter for level-1 network estimation | Step 1 or step 3 of the CAMUS pipeline | doi:10.13012/B2IDB-6892704 (V1 and V2; V2 not inspected): true networks plus 1000 gene trees per replicate; no CAMUS outputs |

**Recommendation: #1, GTM-Blend.** It is the tree-estimation analogue of what you already did with MAGUS's merge step:
- The instructor lists it verbatim as an open problem. The divide-and-conquer lecture says: *"Open problem: Develop a better DTM approach that allows blending."* The projects page says: *"Develop a new Disjoint Tree Merger … compare to GTM"*.
- The published dataset ships the exact inputs to the merge step (subset trees, guide trees, distance matrices) together with the published outputs of GTM and its competitors.
- One GTM merge takes under a second, so the whole study fits on a laptop.
- #2 is the best backup. It is also instructor-listed, it is even cheaper to run, and published baseline species trees exist.

---

## 1. Course pages and what the course covers

**URLs (fetched 2026-10-09):**
- Main page: https://tandy.cs.illinois.edu/CS581-Fall2026.html
- Lectures (the page the user gave): https://tandy.cs.illinois.edu/CS581-Fall2026-lectures.html
- Homework and readings: https://tandy.cs.illinois.edu/CS581-Fall2026-hw.html
- Project suggestions: https://tandy.cs.illinois.edu/CS581-Sp2025-projects.html. The page itself is headed "CS 581 (Fall 2026) Final Project Suggestions" and is linked from the Fall 2026 page.
- Syllabus: https://tandy.cs.illinois.edu/581-syllabus.html
- Previous offering, including student project slides: https://tandy.cs.illinois.edu/CS581-Sp2025-lectures.html. The full archive of offerings is at https://tandy.cs.illinois.edu/warnow-teaching.html.

**Project rules that matter now:**
- HW #7 is a 2-page proposal *"describing a project feasible in about 4 weeks"*, due **Mon Oct 12, 2026**.
- An in-person meeting with the instructor is required.
- The paper must be at least 3,000 words, in Bioinformatics style but single-column. It is due **Dec 8**, with presentations Dec 1–3.
- Grading weights writing at 30% and content at 70%, and reproducibility is required.
- The suggestions page notes that Bayesian projects are probably too computationally heavy.

**Fall 2026 lecture topics (in order):**
1. Course overview.
2. Newick strings, additive matrices, the Naive Quartet Method and distance methods.
3. Rooted trees.
4. Trees from subtrees (supertrees and DTMs) and maximum parsimony.
5. Maximum likelihood and Bayesian estimation.
6. MSA: Needleman-Wunsch, profile HMMs, and MSA methods in practice.
7. Species trees under the MSC (ASTRAL, ASTRID, MP-EST, NJst, SVDquartets).
8. Phylogenetic networks (PhyloNet, SNaQ, CAMUS).
9. Midterm (Oct 8).
10. Deep learning in phylogenetics (guest lecture, Oct 15). Readings: Sapoval et al. 2022; Braichenko et al. 2025; Zhang et al., MBE 2025; Tang et al. 2024.
11. Language phylogenies.
12. Student paper presentations (Oct 22–Nov 19).

"Time permitting" topics are divide-and-conquer tree estimation (DACTAL, NJMerge/TreeMerge/GTM) and phylogenomics part 2 (GDL: ASTRAL-Pro, DISCO, HGT).

The textbook is *Computational Phylogenetics*, Chapters 1–10. Earlier homework readings included Liu et al., Science 2009 (SATé), the MAGUS paper, Mirarab et al. 2014 (statistical binning), Maddison 1997, and DEPP.

**The instructor's own list of open problems and possible projects** (first-day slides, `CS581-Fall2026-firstday.txt`):
- **MSA:** merging two alignments; consensus alignments.
- **Supertrees:** scale to 100,000+ species with high accuracy; approximation algorithms.
- **Species trees and networks:**
  - handle gene tree heterogeneity from multiple causes, including HGT;
  - theoretical guarantees under stochastic models;
  - scale to large numbers of species with whole genomes.
- **Statistical models of evolution:**
  - new models of tree shape and sequence evolution;
  - new theory on identifiability and consistency;
  - evaluating methods under these new models.
- **Deep learning** to improve large-scale tree and network estimation.
- **Applications:** linguistics, microbiome, protein structure and function, tumor evolution.

**Further lecture-specific open problems:**
- *Phylogenomics part 2, slide 35:*
  - Which other species tree methods are consistent under GDL or DLCOAL?
  - Is ASTRAL-Pro consistent under a random model of rooting and tagging error?
  - What is the sample complexity of ASTRID and NJst?
  - *Can we find a distance correction for GDL so that ASTRID/NJst are consistent?*
  - *Can we improve DISCO?*
  - *Better quartet amalgamation methods?*
  - How do we scale concatenation?
  - What about networks?
- *Divide-and-conquer lecture:* *"GTM does NOT allow blending, and so should be able to be improved. Open problem: Develop a better DTM approach that allows blending."*
- *Networks lecture, CAMUS future work:*
  - scale beyond 200 species;
  - *"Develop other techniques to estimate the underlying tree T"*;
  - more complex networks.

**Spring 2025 projects in the same course (useful to avoid repeating):**

| Student | Project | Outcome |
|---|---|---|
| Utkarsh Sharma | GTM with NJ-LogDet guide and/or subset trees on 1000M1/1000M4 | An evaluation, not a new merger. Subset-tree quality mattered about 9× more than guide-tree quality. |
| Xinyu Gu | Cross-validation-weighted HMM ensembles inside WITCH | **Did not improve** over WITCH. |
| The-Anh Vu | DISCO+QR-STAR (rooting) | — |
| Boyang Sun | learnMSA2 as the MAGUS base method | Much worse. |
| Ian Chen | MAGUS iteration | — |
| Prathik Srinivasan | Tree error of MAGUS alignments | — |
| Sean Liu | CAMUS on linguistic data | — |
| Cindy Zeng | learnMSA2 pretraining | — |

Prior course projects that became papers include FASTRAL, Quintet Rooting, pplacerDC, pplacer-XR (which became SCAMPP), the GDL method comparison, and BAli-Phy alignment.

---

## 2. Opportunities by area

Each entry gives: the problem; the current best methods; the data; the metric; the idea; feasibility; and an honest chance of success with a prior-art check.

### 2.1 Divide-and-conquer and large-scale ML tree estimation

**Background on the methods:**
- **NJMerge:** Molloy & Warnow, *AMB* 14:14 (2019), doi:10.1186/s13015-019-0151-x.
- **TreeMerge:** Molloy & Warnow, *Bioinformatics* 35(14):i417 (2019), doi:10.1093/bioinformatics/btz344.
- **Constrained-INC / INC-ML:** Le, Sy, Molloy, Zhang, Rao & Warnow, *IEEE/ACM TCBB* (2021), doi:10.1109/tcbb.2020.2990867.
- **GTM:** Smirnov & Warnow, *BMC Genomics* 21(S2):235 (2020), doi:10.1186/s12864-020-6605-1.
- **DTM pipelines for ML trees:** Park, Zaharias & Warnow, *Algorithms* 14(5):148 (2021), doi:10.3390/a14050148.
- **ML codes:**
  - FastTree 2: doi:10.1371/journal.pone.0009490.
  - RAxML-NG: doi:10.1093/bioinformatics/btz305.
  - IQ-TREE 2: doi:10.1093/molbev/msaa015.
  - VeryFastTree: doi:10.1093/bioinformatics/btaa582.

What these show:
- GTM solves the *unblended* "DTM with guide tree" problem exactly in O(N²) time. The blended version is NP-hard.
- In the GTM paper, the GTM step itself took about 0.4–0.55 s per merge. GTM matched or beat TreeMerge and NJMerge.
- Per the lecture slides, GTM had the most accurate topology on RNASim10k. On RNASim50k, IQ-TREE failed and RAxML had nearly 100% error.

#### A. GTM-Blend: a blending disjoint tree merger (instructor-listed)

**Problem.** The input is disjoint subset trees T₁…T_k plus auxiliary information (a guide tree and/or the alignment). The output is a tree T* with T*|S_i = T_i for every i. "Blending" means the subsets need not form separate regions of T*.

The GTM paper gives an 8-leaf caterpillar for which *no* unblended merge is correct. Its future-work list includes "allow GTM to blend".

**Prior art (checked):**
- NJMerge, TreeMerge and Constrained-INC do blend, but they were less accurate than GTM.
- SDSR (Reshef et al., arXiv:2603.10215, March 2026) is a spectral divide-and-conquer method that merges without blending (its outgroup join).
- No 2023–2026 blending merger that beats GTM was found.
- The Spring 2025 GTM project (Sharma) only varied which methods built the guide and subset trees.

**Idea.** Start from the published GTM tree. Run a hill-climbing search with SPR/NNI moves.
- *Feasibility check.* Keep a move only if every T_i is still induced. Store the bipartitions of each T_i as bitsets and check T*|S_i after each move. The check is O(n·k) per move, which is trivial at n = 1000.
- *Candidate moves.* Restrict to "boundary" moves: subtrees adjacent to the edges that GTM added, regrafted within radius r. These are exactly the moves that interleave neighbouring subsets.
- *Score.* Agreement with the guide tree is not enough. Subsets are carved out of the guide tree itself, so the guide rarely asks for blending. Use a data-driven score instead:
  - (a) Fitch parsimony on the alignment, as a cheap filter;
  - (b) ML log-likelihood, evaluated with RAxML-NG or IQ-TREE on the top candidates only;
  - (c) balanced minimum evolution on LogDet or JC distances, as a variant.
- *Stretch goal (unverified feasibility).* Teach IQ-TREE's constraint check to accept several *partial* constraint trees. Today RAxML-NG and IQ-TREE each take one constraint tree, and the T_i cannot be encoded as one tree without forcing them to be clades.

**Data** (doi:10.13012/B2IDB-7008049_V1, about 34 GB in total; the small conditions are what you need). Download pattern: `https://databank.illinois.edu/datafiles/<web_id>/download`.

| File | Size | Web id |
|---|---|---|
| `RNASim1000.tar.gz` | 9.5 MB | eu0zs |
| `RNASim1000_Analysis.tar.gz` | 338 MB | vyuj4 |
| `1000M1_HF_Analysis.tar.gz` | 203 MB | p4e01 |
| `Cox1-HET.zip` | 3.4 MB | 6446x |
| `Cox1-Het_Analysis.tar.gz` | 267 MB | 1k83z |
| `RNASim10k_R0x_Analysis.tar.gz` | about 0.63 GB per replicate | — |

What the analysis tarballs contain, per the README:
- the centroid decompositions (subset sizes 120 and 500);
- IQ-TREE 2 subset trees;
- FastTree and IQ-TREE starting (guide) trees;
- node and branch-length distance matrices;
- outputs of GTM, Constrained-INC, TreeMerge (PAUP* and RAxML-NG variants), IQ-TREE, RAxML-NG and FastTree.

Two gaps:
- **(unverified):** the 1000M1 true tree may need to come from the SATé data, doi:10.13012/B2IDB-5139418_V1.
- **(unverified):** the replicate counts for RNASim1000 and Cox1-HET.

For a species-tree version, the TreeMerge dataset (doi:10.13012/B2IDB-9570561_V1, 740 MB) has ASTRAL-III subset trees, NJ starting trees, and published NJMerge, NJMerge-2 and TreeMerge species trees. It covers 1000 taxa, 20 replicates, and 10M/500K tree heights.

**Metric.** FN/RF rate against the true tree, plus ML score. Use a paired comparison against the published GTM tree for the same replicate.

**Feasibility.**
- Week 1: rescore the published trees, i.e. reproduce the paper's figures.
- Weeks 2–3: implement the search in Python (dendropy or treeswift plus bitsets) and run it.
- Runtime: minutes per replicate with a parsimony filter and ML scoring of the top candidates. GTM itself takes under a second.
- Add a small caterpillar-tree simulation, as the projects page suggests ("look for a condition where GTM is less accurate than TreeMerge").

**Chance.**
- About 35% that the topology improves clearly (paired) over GTM in some condition.
- About 60% that you can show blending fixes GTM on caterpillar or unbalanced trees.
- Under 15% that it beats IQ-TREE or RAxML-NG at 1,000 taxa.
- The ML score will almost surely improve, since the search starts from GTM and only accepts uphill moves.
- **Not already published (as far as I found).**

#### B. Forest+DTM: an absolute-fast-converging-inspired distance pipeline (course topics: statistical consistency, AFC methods, DTMs)

**Background.** Kim, Lokhov, Vuffray, Romero-Severson & Goldberg, *BMC Bioinformatics* 27:186 (2026), doi:10.1186/s12859-026-06488-y, implemented the "forest" algorithm of Daskalakis, Mossel & Roch (*SIAM J. Discrete Math.* 2011, doi:10.1137/09075576x).
- The forest algorithm returns only the reliably reconstructable components.
- Forest beat NJ only in limited cases, and it returns a forest rather than a tree.
- Their Discussion explicitly suggests using DTMs to merge the forest components. These components are leaf-disjoint, which is exactly what a DTM takes as input.

**Idea.** Take the forest components as constraint trees and merge them with GTM (or GTM-Blend), using an NJ or FastME guide tree. This is a modern DCM-NJ in spirit (Huson, Nettles & Warnow, *JCB* 1999, doi:10.1089/106652799318337).

**Data and code.** github.com/lanl/distphylo contains both the simulator and the Forest code. The data are JC sequences on trees of 10–100 leaves, with 100–100,000 sites.

**Metric.** FN and FP rates against the true tree, as a function of sequence length.

**Feasibility.** Very small. All of it is Python and runs in seconds.

**Chance.** About 35% to beat NJ or FastME at short sequence lengths. **Not done yet**; it is listed only as future work in that paper. The baselines must be rerun, but they are trivial.

### 2.2 Species trees under GDL (phylogenomics part 2)

**Current methods:**
- ASTRAL-Pro: Zhang, Scornavacca, Molloy & Mirarab, *MBE* 37:3292 (2020), doi:10.1093/molbev/msaa139.
- ASTRAL-Pro 2: *Bioinformatics* (2022), doi:10.1093/bioinformatics/btac620. ASTRAL-Pro3 is in the ASTER repository; I found no paper for it **(unverified)**.
- DISCO: Willson, Roddur, Liu, Zaharias & Warnow, *Syst Biol* 71:610 (2022), doi:10.1093/sysbio/syab070.
- FastMulRFS: Molloy & Warnow, *Bioinformatics* 36(S1) (2020), doi:10.1093/bioinformatics/btaa444.
- GDL method comparison (ASTRID-multi vs ASTRAL-Pro, etc.): Willson, Roddur & Warnow, AlCoB 2021, doi:10.1007/978-3-030-74432-8_8.
- wQFM-DISCO: Hakim, Ratul & Bayzid, *Bioinf Adv* (2024), doi:10.1093/bioadv/vbae189.
- A 2025 quartet method for orthologs plus paralogs: Rafi et al., bioRxiv doi:10.1101/2025.04.04.647228.

ASTRID-DISCO is about as accurate as ASTRAL-Pro and sometimes better.

#### C. ASTRID-Pro: a GDL-corrected internode distance (instructor-listed)

**Problem.** As far as I could find, ASTRID-multi (I understand it averages internode distances over all copy pairs; not checked in the code) has no consistency proof under GDL. The lecture asks for "a distance correction for GDL" that would make ASTRID/NJst consistent.

**Idea.**
1. Root each gene-family tree and tag its nodes. A node is a duplication if its two children's species sets overlap; this is ASTRAL-Pro/DISCO style.
2. For each species pair (A, B), take only orthologous copy pairs, i.e. pairs whose LCA is a speciation node.
3. Count only the **speciation** nodes on the path.
4. Average over genes and run FastME.

Variants:
- the closest-copy distance (STAG-like: Emms & Kelly, bioRxiv doi:10.1101/267914, which uses closest copies per gene and then a consensus);
- support weighting, as in weighted ASTRID (Liu & Warnow, *AMB* 18:6, 2023, doi:10.1186/s13015-023-00230-6).

Theory angle: under GDL alone (no ILS), the speciation-only path length between orthologs equals the species-tree path minus loss-pruned nodes. Whether the average is additive is a clean proof question.

**Data** (all with true species trees):
- **FastMulRFS data**, doi:10.13012/B2IDB-5721322_V1, about 2.9 GB.
  - 100 taxa; DL rates 0 to 5×10⁻¹⁰; population sizes 10M and 50M; 10 replicates per condition.
  - RAxML gene trees from 25–250 bp alignments; 25–500 genes.
  - **Published species trees** from ASTRAL-multi, STAG, DupTree, MulRF and FastMulRFS. Its ASTRID-multi trees are flagged as buggy; ignore them.
- **DISCO data**, doi:10.13012/B2IDB-4050038_V1. `trees.tar.gz` is 5.7 GB.
  - 100 species × 1000 genes as the default; also varying GDL, ILS, number of genes, missing data, and a 1000-species condition.
  - Estimated gene trees from 50/100/500 bp. Mean gene tree estimation error is 0.19–0.56 (from `stats.txt`).
- **GDL comparison data**, doi:10.13012/B2IDB-2418574_V1: exp2 is 18 MB and exp3 is 700 MB. Includes ASTRID-multi mapping files.

**Baselines.** ASTRID-multi, ASTRID-DISCO, ASTRAL-Pro 2 and FastMulRFS. Each runs in seconds to minutes per replicate.

**Metric.** Species-tree RF error.

**Feasibility.** Very easy. It is about 300 lines of Python.

**Chance.**
- About 50% to beat ASTRID-multi.
- About 25% to beat ASTRAL-Pro 2 or ASTRID-DISCO.
- A theorem or counterexample is a credible deliverable on its own.
- **Not found published.**

#### D. DISCO-R: species-tree-guided root-and-tag for DISCO (instructor-listed: "Can we improve DISCO?" and "rooting gene trees relative to species trees … NOTUNG")

**Problem.** DISCO roots and tags each tree in isolation, by minimising duplications plus losses with ASTRAL-Pro's heuristic, and then repeatedly cuts off the smallest subtree below each duplication. Errors in rooting and tagging propagate.

**Idea.** Use two passes:
1. Run ASTRAL-Pro or ASTRID-DISCO to get a species tree S₀.
2. Re-root and re-tag each gene-family tree by DL-parsimony reconciliation against S₀. Use LCA mapping over all rootings, which is the classic NOTUNG-style approach. For a non-binary S₀, use only the highly supported branches.
3. Run DISCO's decomposition.
4. Re-estimate the species tree.

You can measure tagging accuracy directly against SimPhy's true gene and locus trees.

**Data.** The DISCO data above, and the DISCO+QR data (doi:10.13012/B2IDB-5748609_V1: 21, 51 and 101 species; six duplication rates; low and high ILS; 2.3 GB).

**Baselines.** DISCO v1.4.1 (Python; `--no-decomp` writes the rooted trees), ASTRAL-Pro 2 and ASTRAL-Pro3.

**Chance.** About 35% for a gain in high-GDL or high-ILS conditions.
- **Prior art:** DISCO+QR and DISCO+QR-STAR (The-Anh, Spring 2025) root the *species* tree, not the gene trees.
- I found no two-pass DISCO.
- ASTRAL-Pro3's internals are **(unverified)**.

### 2.3 Species trees under ILS, supertrees, and quartet amalgamation

**Current methods:**
- ASTRAL-III: doi:10.1186/s12859-018-2129-y.
- wASTRAL: Zhang & Mirarab, *MBE* 2022, doi:10.1093/molbev/msac215.
- TREE-QMC: Han & Molloy, *Genome Res* 2023, doi:10.1101/gr.277629.122.
- wQFM-TREE: *Bioinf Adv* 2024, doi:10.1093/bioadv/vbaf053.
- ASTRID: doi:10.1186/1471-2164-16-S10-S3.
- Weighted ASTRID: see entry C.
- FASTRAL: Dibaeinia, Tabe-Bordbar & Warnow, *Bioinformatics* 2021, doi:10.1093/bioinformatics/btab093.
- FastRFS: doi:10.1093/bioinformatics/btw600.
- Exact-RFS-2 / GreedyRFS: Yu et al., *AMB* 2021, doi:10.1186/s13015-021-00189-2.
- QMC: doi:10.1109/tcbb.2008.133.
- wQMC: doi:10.1093/sysbio/syu087.

#### E. Weighted quartet amalgamation by constrained local search (instructor-listed)

**Idea.** Start from the TREE-QMC or wQMC tree and hill-climb with SPR on the weighted quartet score of an explicit quartet set (from SVDquartets, gene trees, or a deep-learning quartet classifier). You could also add ASTRAL-style dynamic programming over clusters harvested from several heuristic runs.

**Data:**
- HGT+ILS data (Davidson et al. 2015), doi:10.13012/B2IDB-6670066_V1. 2.4 GB; 50 taxa; 6 HGT rates × 50 replicates; true and FastTree gene trees.
- ILS data with missing taxa, doi:10.13012/B2IDB-7735354_V1. 25 taxa; includes **published ASTRAL/ASTRID/MP-EST/SVDquartets species trees**.

**Chance.** About 20%. The field is crowded (TREE-QMC, wQFM-TREE, wASTRAL).

#### F. ASTRID with a correction for missing data

**Background.** Rhodes, Nute & Warnow (arXiv:2001.07844; no venue found) show that ASTRID and NJst are *inconsistent* under i.i.d. taxon deletion. They give no corrected method.

**Idea.** Under i.i.d. deletion with rate p, a node on the path whose off-path subtree has k leaves survives with probability 1−p^k. Re-weight each observed node on the path by the inverse of this survival probability. Alternatively, complete each gene tree with OCTAL (doi:10.1186/s13015-018-0124-5) against a first-pass ASTRID tree.

**Data.** Nute et al. 2018 data, doi:10.13012/B2IDB-7735354_V1 (2 GB; published ASTRID trees).

**Chance.** About 15% to make an empirical difference, because 25-taxon random deletion barely hurts ASTRID. The value is mainly theoretical.

#### G. A supertree-based divide-and-conquer pipeline versus GTM (instructor-listed)

**Idea.** Decompose into overlapping subsets, estimate subset trees, and merge them with GreedyRFS or FastRFS. Compare against GTM on the TreeMerge species-tree data.

**Chance.** About 20%. It is also more engineering work than entry A.

### 2.4 Phylogenetic placement

**Current methods:**
- pplacer: doi:10.1186/1471-2105-11-538.
- EPA-ng: doi:10.1093/sysbio/syy054.
- APPLES: doi:10.1093/sysbio/syz063.
- APPLES-2: *Mol Ecol Resour* 2021, doi:10.1111/1755-0998.13527.
- SCAMPP (pplacer-SCAMPP and EPA-ng-SCAMPP): Wedell, Cai & Warnow, *TCBB* 20(2):1417 (2023), doi:10.1109/TCBB.2022.3170386.
- SCAMPP+FastTree: doi:10.1093/bioadv/vbad008.
- BSCAMPP: Wedell, Shen & Warnow, *IEEE TCBB* 22(4):1593 (2025), doi:10.1109/TCBBIO.2025.3562281.
- DEPP: doi:10.1093/sysbio/syac031.
- uDance: doi:10.1038/s41587-023-01868-8.

**How BSCAMPP works.** Each query votes for its v = 5 closest backbone leaves by Hamming distance. The most-voted leaves seed subtrees of up to B = 2000 leaves. EPA-ng then places the queries on each subtree.

**Headline numbers** (50K-leaf backbone, 10K full-length queries):

| Method | Delta error | Runtime |
|---|---|---|
| SCAMPP(p) | 0.39 | over 27 h |
| SCAMPP(e) | 0.43 | 1,421 min |
| BSCAMPP(e) | 0.45 | about 7 min |

#### H. Better subtree selection in BSCAMPP

**Idea.** Change how the placement subtree is chosen, for example with a union of seeds, k-mer or HMM similarity, or a size that adapts to fragment length.

**Data:**
- SCAMPP data, doi:10.13012/B2IDB-9257957_V1. 70 MB; contains 16S.B.ALL, LTP_s128_SSU, green85 and nt78 with backbones and queries.
- RNASim fragments, doi:10.13012/B2IDB-8788479_V1.
- APPLES data, Dryad doi:10.5061/dryad.78nf7dq.

**Metric.** Delta error.

**Chance.** About 20%. Errors are already below 0.5 edges for full-length queries. Gains are more plausible for fragments. No published placement outputs exist, so the baseline must be rerun with `pip install bscampp`.

### 2.5 Metagenomic profiling (TIPP family)

**Current methods:**
- TIPP: doi:10.1093/bioinformatics/btu721.
- TIPP2: doi:10.1093/bioinformatics/btab023.
- TIPP3 and TIPP3-fast: Shen, Wedell, Pop & Warnow, *PLOS Comput Biol* 21(4):e1012593 (2025), doi:10.1371/journal.pcbi.1012593.
- TIPP-SD (species detection): ACM-BCB 2025, doi:10.1145/3765612.3767205; *PLOS Comput Biol* 2026, doi:10.1371/journal.pcbi.1014347.

**How the TIPP3 variants differ:**

| | Binning | Alignment | Placement | Support threshold | Runtime (16 cores) |
|---|---|---|---|---|---|
| TIPP3 | BLAST | WITCH | pplacer-taxtastic | 90% | 101–312 h |
| TIPP3-fast | BLAST | BLAST | BSCAMPP(e) | 95% | 0.5–5.2 h |

The reference package has 38 marker genes with about 55k sequences each. The metric is normalized Hellinger distance per taxonomic rank. The paper says WITCH, not the placer, is the biggest accuracy factor.

#### I. TIPP3-fast with BSCAMPP(p) and a support-aware label rule (instructor-listed)

**Idea.** Swap pplacer into BSCAMPP (`--placement-method pplacer`). Add a new rule for inheriting labels from the placement, for example likelihood-weighted LCA over the placement mass instead of a fixed threshold. This is the instructor's "taxon identification from a tree in which all but one leaf is labeled" idea.

**Data:**
- Reference package, doi:10.13012/B2IDB-4931852_V2 (2.2 GB; this size is from V1, and I did not check V2's size).
- Benchmark, doi:10.13012/B2IDB-5467027_V1. Simulated reads are 23 GB; genomes and taxonomy are 308 MB.

**Chance.** About 20%. pplacer is only marginally better than EPA-ng, and the compute and data are heavy for a laptop. Treat this as a fallback, not a pick.

### 2.6 Adding sequences to alignments (UPP, UPP2, WITCH, WITCH-NG, EMMA, HMMerge)

**Current methods:**
- UPP: doi:10.1186/s13059-015-0688-z.
- UPP2: doi:10.1093/bioinformatics/btad007.
- WITCH: doi:10.1089/cmb.2021.0585.
- WITCH-NG: doi:10.1093/bioadv/vbad024.
- EMMA: doi:10.1186/s13015-023-00247-x.
- HMMerge: doi:10.1093/bioadv/vbad052.

**Data:**
- ROSE-HF and ROSE-LF: doi:10.13012/B2IDB-6128941_V1 (388 MB).
- 16S.B.ALL-100-HF: doi:10.13012/B2IDB-6604429_V1 (34 MB).
- 5000-het: doi:10.13012/B2IDB-3974819_V1 (1.3 GB).
- EMMA data: doi:10.13012/B2IDB-2567453_V1 (600 MB).
- INDELible length-heterogeneity data: doi:10.13012/B2IDB-0900513_V1 (38 MB).
- UPP data: doi:10.13012/B2IDB-3174395_V1 (3.8 GB).

#### J. Query-specific HMMs for adding sequences

**Idea.** For each query, build one HMM from its k nearest backbone sequences (by k-mer similarity), instead of WITCH's weighted ensemble. Use SPFN and SPFP via your existing FastSP harness.

**Chance.** About 15–20%. The Warnow lab has mined this area heavily, and the Spring 2025 cross-validation-ensemble project failed. The main upside is that your harness already exists.

### 2.7 Phylogenetic networks

#### K. CAMUS: better underlying tree or quartet filter

**Background.** CAMUS: Willson & Warnow, *Bioinformatics* 42(S1):btag245 (ISMB 2026), doi:10.1093/bioinformatics/btag245. Code: Go module github.com/jsdoublel/camus.

The pipeline is:
1. Build an ASTRAL species tree T.
2. Extract the quartets.
3. Filter them with an ILS test at threshold t.
4. Run the dynamic program to get the best level-1 network on T.

CAMUS is close to SNaQ and PhyloNet-MPL at 16–26 species and better than PhyloNet-MPL(FT) at 51.

**Idea.** Either:
- choose the threshold t adaptively per four-taxon set; or
- re-estimate T from quartets after removing those that a first-pass network explains via reticulation, then iterate.

**Data.** doi:10.13012/B2IDB-6892704_V1 (4.3 GB `.tar.xz`). There are 6 conditions × 50 replicates, with the true network, 1000 true gene trees, and FastTree or IQ-TREE gene trees. **No CAMUS outputs** are included. A V2 (2026-08-31) also exists; its API returned an error, so its contents are **(unverified)**.

**Metric.** Cluster FN and FP rates. You must implement network cluster comparison yourself.

**Chance.** About 25%. There is a moving competitor (NetCS, WABI 2026).

### 2.8 Deep learning (the second half of the course)

Data from Zaharias et al., *JCB* 2022 (doi:10.1089/cmb.2021.0383) is at doi:10.13012/B2IDB-8921156_V1 (3.6 GB). A 4-week, laptop-scale *new method* that beats ML or quartet baselines is unlikely (under 10%). **Not recommended**, except as the quartet source for entry E.

---

## 3. Ranking

Scores run from 1 to 5, and higher is better. "Instr." means the topic is on the instructor's list.

| # | Opportunity | Ease | Chance of clear gain | Novelty | Instr. | Published baseline outputs | Total |
|---|---|---|---|---|---|---|---|
| A | GTM-Blend | 5 | 3 | 5 | 5 | Yes (GTM, TreeMerge, C-INC, IQ-TREE, RAxML-NG) | **18** |
| C | ASTRID-Pro (GDL distance) | 5 | 3 | 4 | 5 | Yes (ASTRAL-multi, FastMulRFS, STAG, DupTree, MulRF) | **17** |
| D | DISCO-R (guided root and tag) | 4 | 3 | 4 | 5 | Inputs only; rerun takes minutes | **16** |
| B | Forest+DTM | 4 | 3 | 4 | 3 | Simulator only | **14** |
| K | CAMUS base tree / filter | 3 | 2 | 4 | 4 | Inputs only | **13** |
| G | Supertree divide-and-conquer vs GTM | 3 | 2 | 3 | 4 | Partly | 12 |
| I | TIPP3-fast + BSCAMPP(p) | 2 | 2 | 2 | 5 | Reference package and reads only | 11 |
| E | Quartet amalgamation local search | 3 | 2 | 2 | 4 | Partly (ASTRAL/ASTRID trees) | 11 |
| H | BSCAMPP subtree selection | 4 | 2 | 3 | 2 | No | 11 |
| F | Missing-data ASTRID | 4 | 1 | 3 | 3 | Yes | 11 |
| J | Query-specific HMMs | 4 | 2 | 2 | 3 | No | 11 |

**Top 5:** A, C, D, B, K.

A, C and D share the property you asked for: a small, well-defined component that can be swapped while everything else in the pipeline stays fixed, so the experiment is paired by replicate.

---

## 4. Recommended plan: GTM-Blend in 4 weeks

1. **Week 1: harness.**
   - Download RNASim1000, 1000M1-HF and Cox1-HET (both the data and the `_Analysis` tarballs; under 1 GB in total), and the TreeMerge species-tree data.
   - Rescore every published tree (GTM, TreeMerge, Constrained-INC, IQ-TREE, RAxML-NG, FastTree) against the true trees, and check them against the paper's figures. This is the same validation layer you built for MAGUS.
   - Install GTM from github.com/vlasmirnov/GTM (Python and dendropy) and confirm that it reproduces the published GTM trees.
2. **Week 2: method.**
   - Implement the constrained SPR search with bitset-based checks that every T_i is still induced.
   - Implement a Fitch-parsimony filter and ML re-scoring of the top candidates.
   - Unit-test on the 8-leaf caterpillar counterexample from the GTM paper.
3. **Week 3: experiments.**
   - Run the paired comparison against GTM on every replicate and setting (subset sizes 120 and 500; FastTree and IQ-TREE guide trees).
   - Run a small caterpillar or unbalanced simulation where blending is needed.
   - Run the species-tree setting on the TreeMerge data.
   - Report FN/RF, ML score and runtime.
4. **Week 4: write-up.**
   - A negative result is still a valid paper ("blending does not help when subsets come from the guide tree").
   - Pair it with the caterpillar case where it does help.

Fallback if week 2 stalls: switch to C (ASTRID-Pro). The same tree-handling code carries over.

---

## 5. Things not verified

- The true tree for 1000M1-HF inside the DTM tarballs, and the replicate counts for RNASim1000 and Cox1-HET.
- The contents of CAMUS dataset V2.
- How ASTRID-multi aggregates copies (I believe it averages all copy pairs); the exact STAG distance rule; and ASTRAL-Pro3's rooting and tagging internals.
- Whether IQ-TREE or RAxML-NG can accept several partial constraint trees.
- The Algorithms 2021 DTM paper's numeric results. MDPI blocked the full text, so the GTM accuracy claims above come from the lecture slides and the abstract.

## Sources

- Course: https://tandy.cs.illinois.edu/CS581-Fall2026.html, https://tandy.cs.illinois.edu/CS581-Fall2026-lectures.html, https://tandy.cs.illinois.edu/CS581-Fall2026-hw.html, https://tandy.cs.illinois.edu/CS581-Sp2025-projects.html, https://tandy.cs.illinois.edu/CS581-Sp2025-lectures.html
- Lecture text: `notes/lectures/CS581-Fall2026-firstday.txt`, `CS581-phylogenomics-2026-part2.txt`, `581-Divide-and-Conquer-trees-2023.txt`, `CS581-intro-to-phylogenetic-networks.txt`
- Illinois Data Bank records cited inline (metadata via the `.json` endpoints)
- GTM paper full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC7161100/
- TIPP3: https://pmc.ncbi.nlm.nih.gov/articles/PMC11970662/
- BSCAMPP preprint: https://doi.org/10.1101/2022.10.26.513936
- SDSR: https://arxiv.org/abs/2603.10215
- Kim et al. 2026: https://pmc.ncbi.nlm.nih.gov/articles/PMC13483592/
