# Prior art check: quartet amalgamation / summary species-tree methods

Scope: a short literature check (about 30 minutes of web search, 2026-10-09) for a CS581 pilot. It tests three candidate ideas:
- **(a)** Local search (SPR/NNI hill-climbing) on the (weighted) quartet score, starting from a TREE-QMC/wQMC/ASTRAL tree.
- **(b)** Harvesting bipartitions from several heuristic runs (ASTRAL, TREE-QMC, wQFM-TREE, ASTRID, bootstrap/jackknife gene-tree subsets) to enlarge ASTRAL's restricted search space X, then running the exact DP.
- **(c)** Downweighting quartets or gene trees likely affected by HGT, or other HGT-robust quartet methods.

Legend: **[verified]** means checked against the paper, PMC, or official docs during this session. **[memory]** means the citation or DOI comes from prior knowledge and should be double-checked before citing.

---

## 1. Quartet-amalgamation / summary methods

### TREE-QMC
Han Y, Molloy EK. *Improving quartet graph construction for scalable and accurate species tree estimation from gene trees.* Genome Research 33(7):1042-1052 (2023). doi:10.1101/gr.277629.122 [verified]
- Builds on Snir & Rao's Quartet Max Cut. It constructs the quartet graph directly from gene trees (no explicit quartet list), using normalization schemes for quartets on artificial taxa, in O(n^3 k) time.
- (a) No local search; it is a single divide-and-conquer run. (b) No. (c) No.

### Weighted TREE-QMC ("wTREE-QMC", TREE-QMC v2/v3)
Han Y, Molloy EK. *Improved robustness to gene tree incompleteness, estimation errors, and systematic homology errors with weighted TREE-QMC.* Systematic Biology (2025) syaf009. doi:10.1093/sysbio/syaf009. Preprint doi:10.1101/2024.09.27.615467. [verified title/venue; DOI pattern from dissertation citation]
- Adds wASTRAL-style weighting to TREE-QMC: quartets are weighted by gene-tree branch support, branch length, or a "hybrid" of the two. Later versions also add a mode for unrooted quartet input and multi-copy input.
- (a) No local search, as far as the paper and README show. (b) No. (c) Weighting targets estimation error, missing data, and paralogy or homology errors, **not HGT** specifically.

### wQFM and wQFM-TREE
Mahbub M, Wahab Z, Reaz R, Rahman MS, Bayzid MS. *wQFM: highly accurate genome-scale species tree estimation from weighted quartets.* Bioinformatics 37(21):3734-3743 (2021). doi:10.1093/bioinformatics/btab428 [memory]
Rafi A, Rumi AMS, Hakim SA, Tahmid MT, Momin RJI, Zaman TA, Reaz R, Bayzid MS. *wQFM-TREE: highly accurate and scalable quartet-based species tree inference from gene trees.* Bioinformatics Advances 5:vbaf053 (2025). doi:10.1093/bioadv/vbaf053 [verified]
- wQFM is a divide-and-conquer bipartition method with a Fiduccia-Mattheyses-style refinement **of each bipartition** against an explicit weighted quartet set. wQFM-TREE avoids generating Θ(n^4) quartets and runs in about O(n^3 k log n). It beats ASTRAL in 25/27 conditions (200-1000 taxa).
- (a) **Partial.** The FM-style moves are local search over taxon bipartitions *inside* the divide step. They are not SPR/NNI moves on a complete tree, and there is no post-hoc refinement of a finished tree. (b) No. (c) No.
- Related: Hasan NB, Biswas A, Wahab Z, Mahbub M, Reaz R, Bayzid MS. *Leveraging weighted quartet distributions for enhanced species tree inference from genome-wide data.* Genome Biol Evol (2025), PMC12401674. It weights quartets by their frequency across Bayesian or bootstrap gene-tree distributions. That is weighting for uncertainty, not HGT. DOI not retrieved.

### wQMC / QMC
Avni E, Cohen R, Snir S. *Weighted quartets phylogenetics.* Systematic Biology 64(2):233-242 (2015). doi:10.1093/sysbio/syu087 [memory]
- Runs Quartet Max Cut on an explicit set of weighted quartets.
- (a) No. (b) No. (c) No. It did, however, prove highly accurate under high HGT in Davidson et al. 2015 (see Section 3).

### wASTRAL
Zhang C, Mirarab S. *Weighting by gene tree uncertainty improves accuracy of quartet-based species trees.* Mol Biol Evol 39(12):msac215 (2022). doi:10.1093/molbev/msac215 [verified]
- Weights each gene-tree quartet by branch support ("-u"), by terminal branch length ("-l"), or by both ("hybrid", the default). This gives threshold-free contraction. It is implemented in ASTER as `wastral` and uses ASTER's placement-based search, not ASTRAL-III's DP.
- (a) Uses ASTER's search (placement + NNI + DP; see below). (b) Same as ASTER. (c) **No.** Weights model gene-tree estimation error and are motivated by ILS. Nothing is MSC-deviation or HGT-based.

### ASTRAL-III (and ASTRAL-II search space)
Zhang C, Rabiee M, Sayyari E, Mirarab S. *ASTRAL-III: polynomial time species tree reconstruction from partially resolved gene trees.* BMC Bioinformatics 19(Suppl 6):153 (2018). doi:10.1186/s12859-018-2129-y [memory]
Mirarab S, Warnow T. *ASTRAL-II.* Bioinformatics 31(12):i44-i52 (2015). doi:10.1093/bioinformatics/btv234 [memory]
- Runs an exact DP over a restricted cluster set X. By default X contains the gene-tree bipartitions plus heuristic additions: greedy-consensus polytomy resolutions (capped by `--polylimit`) and similarity- or UPGMA-based completions from ASTRAL-II. ASTRAL-III bounds the size of X so that runtime is polynomial.
- Options [verified from the tutorial]:
  - `-e` / `-f` input "extra trees … to expand its search space"; these do **not** contribute to the score.
  - `-x` runs the exact (unrestricted) DP, feasible up to about 18 taxa.
  - `-q` scores a given tree.
  - In ASTRAL 5.x there is also `--extraLevel` / `-p` to control how many extra bipartitions are added [memory, unverified].
- (a) No. (b) **The mechanism exists (`-e`)**, and adding ASTRID or other method trees via `-e` is an old, documented recommendation (the ASTRID paper by Vachaspati & Warnow 2015 suggests it). No paper systematically harvests from *multiple methods* (TREE-QMC + wQFM-TREE + ASTER + ...). (c) No.

### ASTRAL-MP
Yin J, Zhang C, Mirarab S. *ASTRAL-MP: scaling ASTRAL to very large datasets using randomization and parallelization.* Bioinformatics 35(20):3961-3969 (2019). doi:10.1093/bioinformatics/btz211 [memory]
- Parallelizes ASTRAL-III on GPU and vector units. It also expands and limits X by repeatedly running randomized ASTRAL-style heuristics on **subsamples of taxa**: a "randomized sampling" of clusters added to X.
- (a) No. (b) **Partially.** X is enlarged from many randomized subsample-based runs of ASTRAL's own heuristics, but not from runs of other methods. (c) No.

### FASTRAL
Dibaeinia P, Tabe-Bordbar S, Warnow T. *FASTRAL: improving scalability of phylogenomic analysis.* Bioinformatics 37(16):2317-2324 (2021). doi:10.1093/bioinformatics/btab093 [verified]
- Builds X as the union of bipartitions of **ASTRID trees computed on random gene-tree subsamples**. The setup uses m = 51 subsamples: all genes, 10×50%, 20×25% and 20×10%. Polytomies are resolved with UPGMA, and then the ASTRAL DP runs restricted to that X. ASTRAL was modified so that gene-tree bipartitions are **not** added.
- The goal is *speed*: a smaller X with accuracy comparable to ASTRAL.
- (a) No. (b) **Essentially yes for a single method (ASTRID)** with jackknife-style gene subsets. The difference is that FASTRAL *shrinks* X, while (b) would *enlarge* the default X with bipartitions from several different methods. (c) No.
- Liu B, Warnow T. *Scalable species tree inference with external constraints* (J Comput Biol 2023; preprint doi:10.1101/2021.11.05.467436). FASTRAL-J / NJst-J use user constraint trees. Not directly relevant.

### ASTER package: ASTRAL-IV, ASTRAL-Pro 2/3, CASTER, wASTRAL
Zhang C, Nielsen R, Mirarab S. *ASTER: a package for large-scale phylogenomic reconstructions.* Mol Biol Evol 42(8):msaf172 (2025). doi:10.1093/molbev/msaf172 [verified]
Zhang C, Nielsen R, Mirarab S. *CASTER: direct species tree inference from whole-genome alignments.* Science (2025). doi:10.1126/science.adk9688 [memory]
Zhang C, Mirarab S. *ASTRAL-Pro 2: ultrafast species tree reconstruction from multi-copy gene family trees.* Bioinformatics 38(21):4949-4950 (2022). doi:10.1093/bioinformatics/btac620 [memory]
Zhang C, Scornavacca C, Molloy EK, Mirarab S. *ASTRAL-Pro.* Mol Biol Evol 37(11):3292-3307 (2020). doi:10.1093/molbev/msaa139 [memory]

**Algorithm, confirmed from the MBE 2025 text and the ASTRAL-IV tutorial (`tutorial/astral4.md`):**
1. Greedy initial trees by **sequential placement of species in random orders** into a growing tree. Each placement is optimal for the objective.
2. **NNI moves to improve each initial tree.** The paper says, verbatim: "nearest-neighbor interchange (NNI) moves to improve each initial tree".
3. **DP "akin to ASTRAL" over tripartitions harvested from the initial (and refined) trees.** This finds the best tree within the union of the harvested tripartitions.
4. **Rounds:** "R (4 by default) rounds of search", then "repeatedly performs S (4 by default) rounds of subsampling and exploration until no improvement found". `-R` runs more rounds. Divide-and-conquer is used for speed.
5. Option `-g` takes **user-supplied candidate or guide trees as hints** that are added to the search; this is the ASTER analogue of ASTRAL's `-e`. Option `-c` places taxa onto a given tree, and `-C` / `-r 1 -s 0` scores a fixed tree.

- So ASTER **already** combines local search (NNI, not SPR) with harvesting tripartitions from many randomized runs and then an exact DP. Your recollection is correct, with one addition: NNI refinement is part of each round. The SPR-style "placement" happens only during stepwise addition, not as prune-and-regraft of existing subtrees. That last point is my reading of the paper; it does not mention SPR.
- (a) **Largely done (NNI) inside ASTER.** (b) **Done internally** across its own randomized runs. `-g` lets the user inject trees from other methods, but no published study evaluates multi-method harvesting. (c) No HGT weighting. ASTRAL-Pro is argued to be robust to HGT only because quartet methods are, with no explicit modelling.

### Q-SPR (most direct prior art for (a))
Arasti S, Mirarab S. *Median quartet tree search algorithms using optimal subtree prune and regraft.* Algorithms Mol Biol 19:12 (2024). doi:10.1186/s13015-024-00257-3. WABI 2023 version: LIPIcs WABI 2023 paper 4. [verified]
- Places a pruned subtree optimally (for quartet score) in quasi-linear time using a hierarchical decomposition tree. This drives **SPR hill-climbing on the median quartet tree score** against k reference (gene) trees. Nodes are visited in a random or heuristically weighted order, and the search stops after a round with no improvement.
- It was **started from ASTRAL-III and ASTER trees** (experiment E4) and from stepwise-addition trees (E2).
- Results: the quartet score improves in a minority of cases. In E2 it beat ASTRAL-III in 66 of 294 cases (+0.012% on average); in E4 it improved 129 of 600 replicates, mostly by less than 0.5%. Accuracy (nRF) **did not improve consistently**: in E4, 45 better versus 84 worse. Q-SPR is about 40× slower than ASTRAL-III.
- **Unweighted** quartets only. No TREE-QMC or wQFM starting trees, and no wASTRAL-style weights. HGT appears only in its 10k-taxon dataset (ILS+HGT).
- (a) **Done for unweighted gene-tree quartets from ASTRAL/ASTER starts. The authors found little accuracy gain.**

---

## 2. Site-based, network and other quartet methods
- **SVDquartets.** Chifman J, Kubatko LS. Bioinformatics 30(23):3317-3324 (2014). doi:10.1093/bioinformatics/btu530 [memory]. Infers site-based quartets, which are then assembled with QFM (PAUP*) or QMC. Not (a), (b) or (c).
- **Deep-learning quartet classifier.** Zou Z, Zhang H, Guan Y, Zhang J. *Deep residual neural networks resolve quartet molecular phylogenies.* Mol Biol Evol 37(5):1495-1507 (2020). doi:10.1093/molbev/msz307 [memory]. Per-quartet topology from an alignment. Could supply quartet weights for wQMC or TREE-QMC; it has no HGT-aware component.
- **QuCo.** Rabiee M, Mirarab S. *QuCo: quartet-based co-estimation of species trees and gene trees.* Bioinformatics 38(Suppl 1):i413-i421 (2022). doi:10.1093/bioinformatics/btac265 [verified]. Per-quartet MSC likelihood using gene-tree posteriors, then updates the gene-tree posteriors. It targets gene-tree estimation error (for example, long-branch attraction), **not HGT**. (c) No, although the reweighting structure is analogous.
- **MSCquartets: quartet concordance tests and the T3/star tests, plus QDC, WQDC and NANUQ.** Rhodes JA, Baños H, Mitchell JD, Allman ES. *MSCquartets 1.0.* Bioinformatics 37(12):1766-1768 (2021). doi:10.1093/bioinformatics/btaa868 [memory]. Allman ES, Mitchell JD, Rhodes JA. *Gene tree discord, simplex plots, and statistical tests under the coalescent.* Syst Biol 71(4):929-942 (2022). doi:10.1093/sysbio/syab008 [memory]. These test each 4-taxon set's concordance factors against MSC expectations (tree-like versus star versus non-tree-like). That provides the **statistic** for (c), but it is used for network inference and diagnostics, not for downweighting inputs to a species-tree search.
- **NANUQ.** Allman ES, Baños H, Rhodes JA. Algorithms Mol Biol 14:24 (2019). doi:10.1186/s13015-019-0159-2 [verified DOI]. Builds a level-1 network from quartets that fail tree-likeness tests.
- **TINNiK.** Allman ES, Baños H, Mitchell JD, Rhodes JA. *TINNiK: inference of the tree of blobs of a species network under the coalescent model.* Algorithms Mol Biol (2024). doi:10.1186/s13015-024-00266-2 [verified]. Infers a tree of blobs, using quartet tests to separate tree-like structure from reticulate "blobs". It is HGT/hybridization-aware at the topology level, but it does not reweight quartets for a species-tree search.
- **Quintet / QuNet:** not found as distinct published tools in this search. Possibly confused with MSCquartets routines or Quintet Rooting (Tabatabaee et al. 2022, which is unrelated). Treat as unverified.
- **Stenz NWM, Larget B, Baum DA, Ané C.** *Exploring tree-like and non-tree-like patterns using genome sequences: an example using the inbreeding plant species Arabidopsis thaliana.* Syst Biol 64(5):809-823 (2015). doi:10.1093/sysbio/syv039 [memory]. Uses quartet concordance factors (from BUCKy) and goodness-of-fit to the MSC (the "Quartet MaxCut + CF tests" idea). Diagnostic only.
- **Mallo D, Posada D.** *Multilocus inference of species trees and DNA barcoding.* Phil Trans R Soc B 371:20150335 (2016). doi:10.1098/rstb.2015.0335 [memory]. Mallo also co-authored the SimPhy simulator (Syst Biol 2016, doi:10.1093/sysbio/syv082), which simulates ILS+HGT+GDL. That makes SimPhy useful for (c) experiments.

## 3. HGT theory and robustness of quartet methods
- **Roch S, Snir S.** *Recovering the treelike trend of evolution despite extensive lateral genetic transfer: a probabilistic analysis.* J Comput Biol 20(2):93-112 (2013). doi:10.1089/cmb.2012.0234 [memory DOI]. Under random LGT, with a bounded (almost linear) number of events per gene, the dominant quartet topology equals the species-tree quartet. This makes quartet-majority methods consistent under bounded HGT.
- **Daskalakis C, Roch S.** *Species trees from gene trees despite a high rate of lateral genetic transfer: a tight bound.* SODA 2016. arXiv:1508.01962 [memory: title and venue; arXiv ID verified as the constant-rate HGT result]. Shows the species tree is identifiable from unrooted gene-tree topologies under a constant per-edge HGT rate, with a tight bound on the rate.
- **Davidson R, Vachaspati P, Mirarab S, Warnow T.** *Phylogenomic species tree estimation in the presence of incomplete lineage sorting and horizontal gene transfer.* BMC Genomics 16(Suppl 10):S1 (2015). doi:10.1186/1471-2164-16-S10-S1 [verified]. A simulation of ILS+HGT (50-taxon SimPhy-style). ASTRAL-2 and wQMC stay highly accurate even at high HGT; concatenation and NJst degrade. The authors conjecture that quartet methods are consistent under ILS+HGT.
- **Solís-Lemus C, Yang M, Ané C.** *Inconsistency of species tree methods under gene flow.* Syst Biol 65(5):843-851 (2016). doi:10.1093/sysbio/syw030 [memory]. Shows that quartet methods can be inconsistent under (structured) gene flow or introgression. This is the main motivation for (c).
- **Zhang C, Mirarab S (2022), wASTRAL** (above): not HGT-specific.
- **Other 2023-2026 work:** this search did not find any paper that downweights quartets or gene trees by **deviation from MSC quartet-frequency expectations, or by an HGT likelihood**, and then feeds them into ASTRAL, TREE-QMC or wQFM. HGT-aware work found is either detection (gene-tree/species-tree reconciliation tools such as RANGER-DTL and ALE; quartet-test diagnostics in MSCquartets) or networks (NANUQ, TINNiK, PhyloNet, SNaQ). Quartet-based "highway" HGT detection exists (Bansal, Banay, Gogarten, Shamir 2011, "Detecting highways of horizontal gene transfer", J Comput Biol 18(9):1087, doi:10.1089/cmb.2011.0098 [memory]), but it detects HGT and does not reweight a species-tree search. Coverage of 2025-2026 preprints is limited by the 30-minute budget.

---

## 4. Novelty assessment

**(a) SPR/NNI local search on the (weighted) quartet score from a TREE-QMC/wQMC/ASTRAL start: PARTIALLY DONE, leaning "already done".**
- Q-SPR (Arasti & Mirarab 2024) already does optimal-SPR hill-climbing on the quartet score starting from ASTRAL-III and ASTER trees. It found tiny score gains and no consistent nRF improvement.
- ASTER already runs NNI refinement inside every search round. wQFM uses FM-style local moves within bipartitions.
- **Remaining gaps:**
  - Weighted scores (wASTRAL/wTREE-QMC hybrid weights) as the search objective.
  - Starting from TREE-QMC or wQFM-TREE trees. These are not DP-optimal over any X, so they may leave more room for improvement than ASTRAL/ASTER starts.
  - Explicit-quartet-set input, for example SVDquartets or deep-learning quartet weights, where ASTER-style DP over gene trees is unavailable.
- Expect a small effect size. Evaluate both the score and the RF error, because Q-SPR shows that a better score does not imply a more accurate tree.

**(b) Harvesting bipartitions from multiple heuristics to enlarge X, then exact DP: PARTIALLY DONE.**
- The mechanism exists: ASTRAL `-e` / `-f` and ASTER `-g`.
- FASTRAL builds X from ASTRID runs on jackknife gene-tree subsets, but it *replaces and shrinks* X.
- ASTRAL-MP enlarges X with randomized runs of its own heuristics.
- ASTER internally harvests tripartitions from multiple randomized placement+NNI runs before its DP.
- **Open piece:** a systematic study of *multi-method* harvesting (TREE-QMC + wQFM-TREE + ASTRID + ASTER + bootstrap or jackknife runs) added to the *default* X. This would measure the gain in the score and accuracy of the X-restricted optimum, and the size of X. It is a modest, mostly empirical contribution. The strongest version would show cases where ASTRAL-III's X misses the bipartitions found by TREE-QMC or wQFM-TREE.

**(c) HGT-aware downweighting of quartets or gene trees: LARGELY OPEN.**
- The theory says quartet methods are robust to bounded HGT (Roch & Snir; Daskalakis & Roch; Davidson et al.).
- The MSCquartets tests supply a per-quartet MSC-deviation statistic.
- Existing weighting (wASTRAL, wTREE-QMC, wQFM, Hasan et al. 2025) targets gene-tree estimation error, support and branch length, not HGT.
- **No method was found that uses an MSC-deviation or HGT signal to downweight inputs to a quartet amalgamation search.** This is the most novel of the three ideas.
- **Caveats:**
  - Under bounded HGT, quartet majority is already consistent, so gains may appear only in high-HGT or highway regimes.
  - Use SimPhy with HGT (as in Davidson et al. 2015) to construct those regimes.
  - Compare against plain wASTRAL and wTREE-QMC, and against simply filtering outlier gene trees.
