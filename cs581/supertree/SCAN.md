# Supertree scan (methods and datasets), 2026-10-09

Items marked **(unverified)** were not checked against the paper itself.

## Methods

| method | objective / idea | largest published scale I found | scaling notes | code |
|---|---|---|---|---|
| MRP (PAUP*/TNT) | max parsimony on the Baum-Ragan matrix | ~5,500 taxa (STK supertree); MRP did not finish in 14 days on SMIDGenOG-5500 (cited in the SCS paper) | NP-hard; slow beyond a few thousand taxa | PAUP*, TNT |
| MRL (Nguyen, Mirarab & Warnow, *AMB* 2012, doi:10.1186/1748-7188-7-3) | max likelihood (RAxML) on the MRP matrix | CPL 2,228 | ML heuristic on a binary matrix. **Here I use FastTree on the MRP matrix ("MRL-FT")**, which takes about 2 min at 10k taxa | RAxML / FastTree |
| FastMRP **(unverified citation)** | MRP via a fast parsimony ratchet | — | — | — |
| SuperFine (Swenson et al., *Syst Biol* 2012, doi:10.1093/sysbio/syr092) | Strict Consensus Merger + polytomy refinement (MRP/MRL/QMC) | CPL 2,228; SMIDGen 1,000 | the SCM can collapse to a star at large n (D&C lecture slide) | python2, PAUP* |
| Bad Clade Deletion (Fleischauer & Böcker, *MBE* 2017, doi:10.1093/molbev/msx191; beam search *PeerJ* 2018, doi:10.7717/peerj.4987) | min-cut on the MRP matrix with a weighted clade-deletion objective | SMIDGenOG-5500; SCS-DCM 10,000 (about 2 h per instance with the GSCM preprocessing) | polynomial; ~2 h at 10k | Java (BCD) |
| FastRFS (Vachaspati & Warnow, *Bioinformatics* 2017, doi:10.1093/bioinformatics/btw600) | exact Robinson-Foulds supertree within an ASTRAL-style constrained clade set X | CPL 2,228 (FastRFS-enhanced 3,282 s) | O(n^2 \|X\| k); needs a Bazel build and modified ASTRAL (not built here) | github pranjalv123/FastRFS |
| Exact-RFS-2 / GreedyRFS (Yu et al., *AMB* 2021, doi:10.1186/s13015-021-00189-2) | exact RFS of **two** trees in O(n^2 \|X\|); GreedyRFS merges pairs in SCM order | SMIDGen 500, 20% scaffold, 9 replicates; 501-species ILS D&C | the paper reports criterion scores only, no tree error and no runtimes; beats FastRFS on RFS score when there are few source trees | github yuxilin51/GreedyRFS |
| ASTRAL-II/III as a supertree method (Mirarab & Warnow 2015; Zhang et al. 2018) | max quartet support within the clade set X | SMIDGen 1,000 (RF 11.6–16.9%, FastRFS paper Table 2) | memory and time grow quickly with n and \|X\| | Java 5.7.8 |
| ASTER / ASTRAL-IV (`astral4`, Zhang et al.) | quartet placement plus subsampling search | — | fast. **Finding here: about twice ASTRAL-III's error on SMIDGen supertree input** (see REPORT) | C++ |
| TREE-QMC (Han & Molloy, *Genome Res* 2023, doi:10.1101/gr.277629.122) | weighted quartet Max-Cut, normalized for artificial taxa | the paper calls it "promising as a supertree method" | O(n^3 k); 340 s on one SMIDGen-1000 replicate here | C++ |
| wQFM-TREE (*Bioinf Adv* 2025, doi:10.1093/bioadv/vbaf053) | weighted quartet FM partitioning, run directly on the trees | species-tree data, ≤ ~1,000 taxa | I found no supertree evaluation | Java |
| Asteroid (Morel, Williams & Stamatakis, *Bioinformatics* 2023, doi:10.1093/bioinformatics/btac832) | minimum balanced evolution on internode distances with missing data | thousands of genes; robust at >80% missing | fast | C++ |
| Spectral Cluster Supertree (McArthur et al., *Front Mol Biosci* 2024, doi:10.3389/fmolb.2024.1432495) | spectral clustering of a proper-cluster graph (rooted inputs) | 10,000 taxa, DCM source trees: ~20 s versus ~2 h for BCD | BCD has better RF; SCS has better matching-cluster distance | PyPI `sc-supertree` |
| SDSR (Reshef et al., arXiv:2603.10215, 2026, preprint) | spectral D&C for *species* trees, with MSC recovery guarantees | up to 10× faster than CA-ML/ASTRAL at comparable accuracy (abstract) | — | — |
| DTMs: NJMerge (*AMB* 2019), TreeMerge (*Bioinformatics* 2019, doi:10.1093/bioinformatics/btz344), GTM (Smirnov & Warnow, 2020/21) | merge **disjoint** subset trees using a guide tree or distance matrix | GTM: RNASim 50k (ML trees); GTM+ASTRAL on species trees | GTM runs in O(n) on top of the guide tree; no blending | github vlasmirnov/GTM, ekmolloy/treemerge |

## Datasets with true trees

| dataset | DOI / location | size | contents |
|---|---|---|---|
| SMIDGen (Swenson et al., *AMB* 2010, doi:10.1186/1748-7188-5-8) | Illinois Data Bank doi:10.13012/B2IDB-2952208_V1 (`SuperFine.zip`, 13.6 MB) | 100/500 taxa: 30 replicates per scaffold density; 1,000 taxa: 10 per density; densities 20/50/75/100% | `sm_data.R.model_tree` (true tree) and `sm_data.R.source_trees`: about 25 ML source trees, one scaffold plus clade-based trees. Same archive: CPL, THPL, mammals, marsupials and seabirds, with published SuperFine/MRP/RFS/SCM trees |
| SMIDGenOG / SMIDGenOG-5500 (Fleischauer & Böcker 2016/2017) | Zenodo doi:10.5281/zenodo.11118022 (`scs_analysis_data.zip`, 1.1 GB; members fetched by HTTP range) | 100/500/1,000 taxa × 4 scaffold densities × 30 replicates; the "10000" folder holds 10 replicates of ~5,500 taxa with ~505 RAxML source trees | outgroup-rooted model trees (pruned), RAxML source trees and alignments |
| SCS birth-death DCM-IQ (McArthur et al. 2024) | same Zenodo archive, `birth_death/` | 500 / 1,000 / 2,000 / 5,000 / 10,000 taxa × 10 replicates; Rec-I-DCM3 subsets ≤ 50 or ≤ 100 taxa | model trees, IQ-TREE 2 source trees (10k sites), "exact" DCM source trees (true trees restricted), sequences |
| SuperTriplets benchmark | same Zenodo archive | 101 taxa; 25/50/75% deletion; 100 problems | rooted model and source trees |
| TreeMerge species-tree data (Molloy & Warnow 2019) | doi:10.13012/B2IDB-9570561_V1 (0.74 GB) | 1,000 ingroup taxa, 2 tree heights (10M, 500K), 20 replicates; exons and introns | true species trees; subset decompositions (≤120 taxa) from NJ-AGID, NJ-logdet and RP starting trees; ASTRAL and RAxML subset trees; NJMerge, NJMerge-2 and TreeMerge outputs. A DTM benchmark (disjoint subsets), not a supertree benchmark |
| NJMerge data | doi:10.13012/B2IDB-1424746_V1 (25 GB) | as above, with gene trees and alignments | inputs to the TreeMerge data |
| Exact-RFS-2 / GreedyRFS data | no Data Bank record found; the paper uses SMIDGen-500 (20% scaffold, 9 replicates) from the SuperFine data and a 501-species ASTRAL-II ILS simulation | — | code at github yuxilin51/GreedyRFS |
| DTMs for ML trees (Park et al.) | doi:10.13012/B2IDB-7008049 (34 GB) | RNASim up to 50k | — |
| INC in D&C | doi:10.13012/B2IDB-8518809 (1.6 GB) | — | — |

The only true-tree supertree data I found at **10,000+ taxa** is the SCS birth-death 10k set (and SMIDGenOG-5500 at ~5.5k). No public 100k-taxon supertree benchmark turned up.
