# Prior art: rogue taxa, divergent homologs, and their effect on MSA and tree estimation

Compiled 2026-10-10 for CS581 project scoping. Each DOI was checked against the Crossref API
(api.crossref.org/works/<DOI>), and the title, authors and venue it returned matched the citation.
"verified" means Crossref resolved it. Notes come from abstracts and what I already knew about these
papers; I did not reread the full texts.

## 1. Rogue taxon identification

- **Aberer AJ, Krompass D, Stamatakis A.** Pruning rogue taxa improves phylogenetic accuracy: an efficient algorithm and webservice. *Syst Biol* 62(1):162-166 (2013; online 2012). DOI 10.1093/sysbio/sys078 (verified).
  RogueNaRok. A graph-based algorithm that finds the leaf sets whose pruning most increases the support or resolution of a consensus of bootstrap trees (the RBIC criterion). In simulation, pruning rogues brought consensus and ML trees closer to the true tree. It works on trees only, not on alignments.
- **Pattengale ND, Aberer AJ, Swenson KM, Stamatakis A, Moret BME.** Uncovering hidden phylogenetic consensus in large data sets. *IEEE/ACM TCBB* 8(4):902-911 (2011). DOI 10.1109/TCBB.2011.28 (verified).
  The direct predecessor of RogueNaRok. It defines the RBIC objective and a slower heuristic for it.
- **Aberer AJ, Stamatakis A.** A simple and accurate method for rogue taxon identification. *IEEE BIBM* 2011, pp. 118-122. DOI unverified (I could not confirm it on Crossref).
  Another RogueNaRok predecessor.
- **Thorley JL, Wilkinson M.** Testing the phylogenetic stability of early tetrapods. *J Theor Biol* 200(3):343-344 (1999). DOI 10.1006/jtbi.1999.0999 (verified).
  The usual citation for the leaf stability index, which measures how consistently a leaf's quartets resolve across bootstrap trees. It is in RadCon (Thorley & Page 2000, *Bioinformatics* 16:486, DOI 10.1093/bioinformatics/16.5.486, verified) and in Phyutility.
- **Wilkinson M.** Common cladistic information and its consensus representation: reduced Adams and reduced cladistic consensus trees and profiles. *Syst Biol* 43(3):343-368 (1994). DOI 10.1093/sysbio/43.3.343 (verified).
  Reduced consensus. It finds the agreement that remains once unstable ("wildcard") taxa are left out.
- **Wilkinson M.** Coping with abundant missing entries in phylogenetic inference using parsimony. *Syst Biol* 44(4):501-514 (1995). DOI 10.2307/2413657 (verified; the 10.1093/sysbio/44.4.501 form did not resolve).
  Introduces safe taxonomic reduction (STR): taxa that can be removed a priori without changing the relationships inferred among the rest. Defined for parsimony and missing data.
- **Wilkinson M.** Majority-rule reduced consensus trees and their use in bootstrapping. *Mol Biol Evol* 13(3):437-444 (1996). DOI 10.1093/oxfordjournals.molbev.a025604 (verified).
  Extends reduced consensus to bootstrap support, which is what most rogue-pruning work optimizes.
- **Siu-Ting K, Pisani D, Creevey CJ, Wilkinson M.** Concatabominations: identifying unstable taxa in morphological phylogenetics using a heuristic extension to safe taxonomic reduction. *Syst Biol* 64(1):137-143 (2015). DOI 10.1093/sysbio/syu066 (verified).
  A heuristic STR for unstable taxa. There is follow-up work from the same group: Serra Silva et al. 2025 *Syst Biol*, "Coping with ineffective overlap in multilocus phylogenetics", DOI 10.1093/sysbio/syaf044 (verified).
- **Smith MR.** Using information theory to detect rogue taxa and improve consensus trees. *Syst Biol* 71(5):1088-1094 (2022). DOI 10.1093/sysbio/syab099 (verified).
  Scores taxa by the phylogenetic information content of the consensus. Implemented in the R package `Rogue` (with TreeTools). Its reduced consensus trees are better resolved and more accurate than those from RogueNaRok or LSI.
- **Mai U, Mirarab S.** TreeShrink: fast and accurate detection of outlier long branches in collections of phylogenetic trees. *BMC Genomics* 19(Suppl 5):272 (2018). DOI 10.1186/s12864-018-4620-2 (verified).
  Mirarab lab. Finds leaves that inflate the diameter of each gene tree, then removes them from that gene only. It is more conservative than rogue removal and cuts gene tree discordance more at a matched amount of filtering. The usual pipeline is align, build a tree, TreeShrink, then realign or retree. That loop makes it a natural starting point for "pre-alignment filtering".
- **Comte A, Tricou T, Tannier E, Joseph J, Siberchicot A, Penel S, et al.** PhylteR: efficient identification of outlier sequences in phylogenomic datasets. *Mol Biol Evol* 40(11):msad234 (2023). DOI 10.1093/molbev/msad234 (verified).
  Builds a multi-gene distance-matrix analysis (DISTATIS) that flags gene x species outliers, including outliers caused by topology and not only by long branches. Compared against TreeShrink.
- **de Vienne DM, Ollier S, Aguileta G.** Phylo-MCOA: a fast and efficient method to detect outlier genes and species in phylogenomics using multiple co-inertia analysis. *Mol Biol Evol* 29(6):1587-1598 (2012). DOI 10.1093/molbev/msr317 (verified; note that 10.1093/sysbio/syr130 is a different paper).
  The precursor to PhylteR. It finds complete outlier genes and species, and cell-wise (gene x species) outliers.
- **Struck TH.** TreSpEx: detection of misleading signal in phylogenetic reconstructions based on tree information. *Evol Bioinform* 10:51-67 (2014). DOI 10.4137/EBO.S14239 (verified).
  Uses patristic distances and LB-scores to flag long-branch taxa and genes.
- **Kück P, Meid SA, Groß C, Wägele JW, Misof B.** AliGROOVE: visualization of heterogeneous sequence divergence within multiple sequence alignments and detection of inflated branch support. *BMC Bioinformatics* 15:294 (2014). DOI 10.1186/1471-2105-15-294 (verified).
  Works at the alignment level. It flags single taxa with heterogeneous divergence, which can cause long-branch or rogue behaviour, and that site masking misses.

## 2. Long-branch attraction, taxon sampling, fast-evolving taxa

- **Felsenstein J.** Cases in which parsimony or compatibility methods will be positively misleading. *Syst Zool* 27(4):401-410 (1978). DOI 10.2307/2412923 (verified).
  The original "Felsenstein zone" (LBA) result. Parsimony is statistically inconsistent when two long branches are separated by a short internal edge.
- **Bergsten J.** A review of long-branch attraction. *Cladistics* 21(2):163-193 (2005). DOI 10.1111/j.1096-0031.2005.00059.x (verified).
  The standard review. Remedies it covers include adding taxa that break up long branches, removing long-branch taxa, removing fast sites, better models, and outgroup choice.
- **Hillis DM.** Inferring complex phylogenies. *Nature* 383:130-131 (1996). DOI 10.1038/383130a0 (verified).
  Simulated a 228-taxon tree and found that dense taxon sampling makes large trees tractable.
- **Zwickl DJ, Hillis DM.** Increased taxon sampling greatly reduces phylogenetic error. *Syst Biol* 51(4):588-598 (2002). DOI 10.1080/10635150290102339 (verified).
  Adding taxa that subdivide long branches reduces error more than adding characters does.
- **Pollock DD, Zwickl DJ, McGuire JA, Hillis DM.** Increased taxon sampling is advantageous for phylogenetic inference. *Syst Biol* 51(4):664-671 (2002). DOI 10.1080/10635150290102357 (verified).
  A rebuttal to the view that adding taxa can hurt. This is the counterpoint to "remove rogues".
- **Hedtke SM, Townsend TM, Hillis DM.** Resolution of phylogenetic conflict in large data sets by increased taxon sampling. *Syst Biol* 55(3):522-529 (2006). DOI 10.1080/10635150600697358 (verified).
- **Kück P, Mayer C, Wägele JW, Misof B.** Long branch effects distort maximum likelihood phylogenies in simulations despite selection of the correct model. *PLoS ONE* 7(5):e36593 (2012). DOI 10.1371/journal.pone.0036593 (verified).
  Shows that LBA also affects ML when the model is correct, for finite sequence lengths.
- **Philippe H, Brinkmann H, Lavrov DV, Littlewood DTJ, Manuel M, Wörheide G, et al.** Resolving difficult phylogenetic questions: why more sequences are not enough. *PLoS Biol* 9(3):e1000602 (2011). DOI 10.1371/journal.pbio.1000602 (verified).
  Argues that systematic error from fast-evolving taxa and genes, and from contamination, dominates in phylogenomics. Recommends removing fast-evolving taxa and better models.
- **Sanderson MJ, Shaffer HB.** Troubleshooting molecular phylogenetic analyses. *Annu Rev Ecol Syst* 33:49-72 (2002). DOI 10.1146/annurev.ecolsys.33.010802.150509 (verified).
  A review that covers rogue taxa and taxon deletion.

## 3. Effect of distant, fast-evolving or fragmentary homologs on MSA (Warnow lab and others)

- **Nguyen N-pD, Mirarab S, Kumar K, Warnow T.** Ultra-large alignments using phylogeny-aware profiles. *Genome Biol* 16:124 (2015). DOI 10.1186/s13059-015-0688-z (verified).
  UPP. It picks a backbone of "full-length" sequences (a random subset of length-filtered sequences), aligns it with PASTA, and adds the remaining sequences with an ensemble of HMMs. Other methods lose much more accuracy than UPP when fragments are present. The design idea is to keep problematic sequences out of the backbone, a form of filtering by length.
- **Park M, Ivanovic S, Chu G, Shen C, Warnow T.** UPP2: fast and accurate alignment of datasets with fragmentary sequences. *Bioinformatics* 39(1):btad007 (2023). DOI 10.1093/bioinformatics/btad007 (verified).
  Selects the HMM from the ensemble for each query sequence. It is more accurate than UPP and other leading MSA methods under sequence length heterogeneity.
- **Smirnov V, Warnow T.** MAGUS: Multiple sequence Alignment using Graph clUStering. *Bioinformatics* 37(12):1666-1672 (2021). DOI 10.1093/bioinformatics/btaa992 (verified).
  Divide and conquer: it aligns subsets and merges them with the Graph Clustering Merger. On its own it does not handle fragments well.
- **Smirnov V, Warnow T.** Phylogeny estimation given sequence length heterogeneity. *Syst Biol* 70(2):268-282 (2021). DOI 10.1093/sysbio/syaa058 (verified).
  Introduces MAGUS+UPP (MAGUS backbone of full-length sequences, fragments added with UPP). Fragments degrade both alignment and tree accuracy, and FastTree does poorly on alignments with fragments.
- **Shen C, Zaharias P, Warnow T.** MAGUS+eHMMs: improved multiple sequence alignment accuracy for fragmentary sequences. *Bioinformatics* 38(4):918-924 (2022). DOI 10.1093/bioinformatics/btab788 (verified).
- **Shen C, Park M, Warnow T.** WITCH: improved multiple sequence alignment through weighted consensus hidden Markov model alignment. *J Comput Biol* 29(8):782-801 (2022). DOI 10.1089/cmb.2021.0585 (verified).
  Adds each query sequence using a weighted set of HMMs (a bitscore-based weighting), on top of the UPP framework. WITCH-NG (Liu & Warnow, bioRxiv 2022, DOI 10.1101/2022.08.08.503232, verified) is a faster version.
- **Shen C, Liu B, Williams KM, Warnow T.** EMMA: a new method for computing multiple sequence alignments given a constraint subset alignment. *Algorithms Mol Biol* 18:21 (2023). DOI 10.1186/s13015-023-00247-x (verified).
  Adds sequences to a fixed constraint alignment while keeping it intact, using MAFFT --add on subproblems. This is the natural tool for an experiment that holds the clean taxa's backbone fixed and then adds rogue taxa.
- **Mirarab S, Nguyen N, Guo S, Wang L-S, Kim J, Warnow T.** PASTA: ultra-large multiple sequence alignment for nucleotide and amino-acid sequences. *J Comput Biol* 22(5):377-386 (2015). DOI 10.1089/cmb.2014.0156 (verified).
- **Mirarab S, Warnow T.** FastSP: linear time calculation of alignment accuracy. *Bioinformatics* 27(23):3250-3258 (2011). DOI 10.1093/bioinformatics/btr553 (verified).
  Computes SP-score, modeler/developer scores, SPFN and SPFP. Restricting both the estimated and the reference MSA to the retained taxa before scoring answers the core project question.
- **Nute M, Warnow T.** Scaling statistical multiple sequence alignment to large datasets. *BMC Genomics* 17(Suppl 10):764 (2016). DOI 10.1186/s12864-016-3101-8 (verified).
- **Katoh K, Frith MC.** Adding unaligned sequences into an existing alignment using MAFFT and LAST. *Bioinformatics* 28(23):3144-3146 (2012). DOI 10.1093/bioinformatics/bts578 (verified).
  MAFFT --add, which keeps the existing alignment fixed.
- **Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, et al.** Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. *Mol Syst Biol* 7:539 (2011). DOI 10.1038/msb.2011.75 (verified).
- **Sievers F, Dineen D, Wilm A, Higgins DG.** Making automated multiple alignments of very large numbers of protein sequences. *Bioinformatics* 29(8):989-995 (2013). DOI 10.1093/bioinformatics/btt093 (verified).
  Progressive aligners lose accuracy on a fixed reference subset (the HomFam design) as more homologs are added. This is the closest existing measurement of "accuracy of the core set versus how many extra sequences are added". The variable there is the number of extra sequences, not how divergent they are.
- **Boyce K, Sievers F, Higgins DG.** Simple chained guide trees give high-quality protein multiple sequence alignments. *PNAS* 111(29):10556-10561 (2014). DOI 10.1073/pnas.1405628111 (verified).
  Rebutted by Tan G, Gil M, Löytynoja A, Goldman N, Dessimoz C, *PNAS* 112(2):E99 (2015), DOI 10.1073/pnas.1417526112 (verified). Both use the core-subset scoring design, and they disagree about how guide trees should be judged.
- **Eddy SR.** Accelerated profile HMM searches. *PLoS Comput Biol* 7(10):e1002195 (2011). DOI 10.1371/journal.pcbi.1002195 (verified).
  HMMER3. Bitscores or E-values against a profile of the backbone can flag non-homologous or contaminant sequences. UPP and WITCH rely on this implicitly.
- **Castresana J.** Selection of conserved blocks from multiple alignments for their use in phylogenetic analysis. *Mol Biol Evol* 17(4):540-552 (2000). DOI 10.1093/oxfordjournals.molbev.a026334 (verified). This is Gblocks.
- **Capella-Gutiérrez S, Silla-Martínez JM, Gabaldón T.** trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. *Bioinformatics* 25(15):1972-1973 (2009). DOI 10.1093/bioinformatics/btp348 (verified). trimAl also has a `-resoverlap/-seqoverlap` option that removes whole sequences, not just columns.
- **Criscuolo A, Gribaldo S.** BMGE: a new software for selection of phylogenetic informative regions from multiple sequence alignments. *BMC Evol Biol* 10:210 (2010). DOI 10.1186/1471-2148-10-210 (verified).
- **Tan G, Muffato M, Ledergerber C, Herrero J, Goldman N, Gil M, Dessimoz C.** Current methods for automated filtering of multiple sequence alignments frequently worsen single-gene phylogenetic inference. *Syst Biol* 64(5):778-791 (2015). DOI 10.1093/sysbio/syv033 (verified).
  Column filtering usually makes single-gene trees worse. This is a caution for any filtering-helps hypothesis, though it concerns columns, not taxa.
- **Steenwyk JL, Buida TJ, Li Y, Shen X-X, Rokas A.** ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference. *PLoS Biol* 18(12):e3001007 (2020). DOI 10.1371/journal.pbio.3001007 (verified).

## 4. Sequence- and taxon-level outlier detection in alignments

- **Jehl P, Sievers F, Higgins DG.** OD-seq: outlier detection in multiple sequence alignments. *BMC Bioinformatics* 16:269 (2015). DOI 10.1186/s12859-015-0702-1 (verified).
  Flags sequences whose mean gap-based distance to the others is anomalous. On MSA benchmarks, removing them did not significantly change accuracy. That neutral result is directly relevant to the project.
- **Di Franco A, Poujol R, Baurain D, Philippe H.** Evaluating the usefulness of alignment filtering methods to reduce the impact of errors on evolutionary inferences. *BMC Evol Biol* 19:21 (2019). DOI 10.1186/s12862-019-1350-2 (verified).
  Introduces HmmCleaner, which masks sequence segments that fit a profile HMM poorly. Compares segment filtering with block filtering.
- **Zhang C, Zhao Y, Braun EL, Mirarab S.** TAPER: pinpointing errors in multiple sequence alignments despite varying rates of evolution. *Methods Ecol Evol* 12(11):2145-2158 (2021). DOI 10.1111/2041-210X.13696 (verified).
  Mirarab lab. Finds outlier stretches in each species while separating real fast evolution and long branches from error. The paper's warning that divergence is not error applies to rogue removal as well.
- **Whelan S, Irisarri I, Burki F.** PREQUAL: detecting non-homologous characters in sets of unaligned homologous sequences. *Bioinformatics* 34(22):3929-3930 (2018). DOI 10.1093/bioinformatics/bty448 (verified).
  Filters residues before alignment, the counterpart to removing taxa before alignment. A related tool is Divvier (Ali, Bogusz, Whelan, *MBE* 36:2340, 2019, DOI 10.1093/molbev/msz142, verified).
- **Borowiec ML.** Spruceup: fast and flexible identification, visualization, and removal of outliers from large multiple sequence alignments. *JOSS* 4(42):1635 (2019). DOI 10.21105/joss.01635 (verified).
  Removes per-taxon outlier windows by distance.
- **Tumescheit C, Firth AE, Brown K.** CIAlign: a highly customisable command line tool to clean, interpret and visualise multiple sequence alignments. *PeerJ* 10:e12983 (2022). DOI 10.7717/peerj.12983 (verified).
  Includes removal of divergent sequences (a percent-identity-to-consensus threshold) as well as column cleaning.
- **Dunn CD.** SequenceBouncer: a method to remove outlier entries from a multiple sequence alignment. *bioRxiv* (2020). DOI 10.1101/2020.11.24.395459 (verified; a preprint, and I found no journal version).
  An entropy-based, IQR-thresholded outlier removal. After realignment there are fewer heavily gapped columns, but it reports no accuracy-versus-reference measurement.
- **Chiner-Oms A, González-Candelas F.** EvalMSA: a program to evaluate multiple sequence alignments and detect outliers. *Evol Bioinform* 12:EBO.S40583 (2016). DOI 10.4137/EBO.S40583 (verified).
  Estimates how much MSA quality would improve if the flagged sequences were dropped.

## 5. Work from 2020 to 2026 on taxon removal, adding homologs, and alignment accuracy

- **Kim D, Gil M, Katoh K, Dessimoz C.** AmpliPhy improves gene trees by adding homologous sequences without affecting alignments. *bioRxiv* (2026). DOI 10.64898/2026.01.26.701724 (verified via Crossref; a preprint).
  The closest recent study. It uses an enrichment-impoverishment design: add homologs, then score the original sequences. Enrichment consistently improves the gene trees, but the effect on alignment quality of the original set is only marginal. Adding sequences with MAFFT --add keeps most of the benefit. This is the opposite direction to rogue removal, so it is a key comparison point.
- **Sayyari E, Whitfield JB, Mirarab S.** Fragmentary gene sequences negatively impact gene tree and species tree reconstruction. *Mol Biol Evol* 34(12):3279-3291 (2017). DOI 10.1093/molbev/msx261 (verified).
  Mirarab lab. Removing fragmentary sequences (filtering by taxon or sequence) improves gene trees and species trees. This predates 2020 but is the most direct evidence that removing sequences helps tree accuracy.
- In the Warnow lab, the 2020-2023 line of work is UPP2, WITCH, WITCH-NG, MAGUS+eHMMs, EMMA, and Smirnov & Warnow 2021 (all above). It handles problem sequences by keeping them out of the backbone, not by deleting them. Its benchmarks score all sequences, or report fragment and full-length accuracy together. I found no Warnow-lab paper that scores the accuracy of the remaining taxa as divergent or rogue taxa are added or removed.
- Several searches ("rogue taxa alignment accuracy", "divergent sequences MSA accuracy remaining", "outlier sequence removal multiple sequence alignment", "TreeShrink alignment") found no other dedicated 2020-2026 study. The search was limited to web search, so Google Scholar coverage was not exhaustive.

## Gap

I found no study that directly measures how including rogue or long-branch homologs changes the
SP-score (or SPFN/SPFP) of the MSA *restricted to the remaining taxa*, or that tests whether
pre-filtering such taxa before alignment recovers that accuracy. Three lines of work come close.
Sievers et al. 2013 and Boyce et al. 2014 track accuracy on a fixed reference subset as more
homologs are added, but the variable is the number of added sequences, not how divergent or rogue
they are. AmpliPhy (2026 preprint) finds almost no effect on alignment of the original set when
homologs are added. OD-seq reports a neutral benchmark effect when outliers are removed. The
Warnow-lab methods (UPP, UPP2, MAGUS+UPP, WITCH, EMMA) put length outliers in the "query" set
rather than the backbone, but they do this for fragments, not for long-branch or rogue sequences,
and they do not report accuracy on the retained subset. Tree-level filters (RogueNaRok, Smith's
Rogue, TreeShrink, PhylteR) are judged on tree accuracy only. Tan et al. 2015 warns that filtering
can hurt, though their filtering was by column.

A controlled experiment would fill this gap: simulate with long branches inserted (for example,
with INDELible, RAxML-NG simulation or Alisim). Then compare three settings:
(a) align all taxa and restrict to the clean taxa;
(b) remove rogues (identified by TreeShrink, RogueNaRok or OD-seq), realign, and score;
(c) build a backbone from the clean taxa and add the rogues with EMMA or UPP.
Score each with FastSP and with tree error on the shared taxa. As far as I can tell this
has not been published, but my search was not exhaustive (no Google Scholar, no full-text search).
