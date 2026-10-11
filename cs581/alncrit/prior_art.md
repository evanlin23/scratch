# Prior art: which MSA accuracy criteria predict tree accuracy; effect of over-alignment / compression

DOIs were checked against Crossref and PubMed/Europe PMC (Oct 2026). Claims were checked against the full text where it was
open access (OA) or against the abstract otherwise. These are flagged as [FT] = read full text, [Abs] = abstract only, [Mem] = recalled from memory and not re-verified.
Terminology: SP-score = 1-SPFN (recall); Modeler score = 1-SPFP (precision); TC = total-column score; tree error = FN/RF.

## Warnow-group simulation studies (nucleotide, ROSE/Indelible/RNASim; RAxML/FastTree)

- **Nelesen S, Liu K, Zhao D, Linder CR, Warnow T (2008).** The effect of the guide tree on multiple sequence alignments
  and subsequent phylogenetic analyses. *Pac Symp Biocomput* 13:25-36. DOI 10.1142/9789812776136_0004 (verified, Crossref). [FT]
  Here "SP-error" means SPFN only. Better guide trees barely changed SP-error but made RAxML trees much more accurate (ProbCons, FTA).
  Quote: "not all errors are of equal importance ... although FTA alignments are often 'worse' with respect to the SP-error,
  trees estimated from FTA alignments can be more accurate than trees estimated from other alignments with lower SP-error";
  alignments with similar SP-error "can still be very different from each other." This is the earliest clear evidence that SPFN alone does not predict tree error.
- **Liu K, Raghavan S, Nelesen S, Linder CR, Warnow T (2009).** Rapid and accurate large-scale coestimation of sequence
  alignments and phylogenetic trees. *Science* 324(5934):1561-1564. DOI 10.1126/science.1171243 (verified). [Abs; details Mem]
  SATé improved both alignment SPFN and ML tree FN rate over two-phase methods on 1000-taxon simulations. The paper reports
  SPFN as its main alignment criterion. It does not quantify the alignment-tree correlation, and I could not confirm whether its supplement compares SPFN with SPFP.
- **Liu K, Linder CR, Warnow T (2010).** Multiple sequence alignment: a major challenge to large-scale phylogenetics.
  *PLoS Currents* 2:RRN1198. DOI 10.1371/currents.RRN1198 (verified). [FT]
  Quote: "alignment error, measured using SP-FN, is not particularly predictive of tree error". Opal and Prank+GT swap rank positions between
  the two criteria. "ML trees on Prank+GT alignments have reasonably low missing branch rates although Prank+GT alignments have high SP-FN error",
  and Opal is the reverse. This is direct evidence that a gappy, under-aligned (low-FP) alignment can give good trees despite high SPFN.
- **Liu K, Warnow T, Holder MT, Nelesen SM, Yu J, Stamatakis AP, Linder CR (2012).** SATé-II: very fast and accurate
  simultaneous estimation of multiple sequence alignments and phylogenetic trees. *Syst Biol* 61(1):90-106.
  DOI 10.1093/sysbio/syr095 (verified). [Abs]
  SATé-boosted alignments are more accurate, and so are the trees built on them. The authors found a correlation between tree/alignment pair quality and ML score
  with gaps treated as missing data. However, directly optimizing that ML score gives very poor alignments and trees, and every tree is optimal under JC.
  So the likelihood (gaps treated as missing) is not a usable alignment criterion. Relevant because over-alignment is cheap under gaps-as-missing ML.
- **Liu K, Warnow T (2012).** Treelength optimization for phylogeny estimation. *PLoS ONE* 7(3):e33104.
  DOI 10.1371/journal.pone.0033104 (verified). [FT]
  BeeTLe finds shorter treelength than POY under three gap penalties (Simple-1, Simple-2, Affine), but its trees are worse than ML on standard alignments.
  Treelength-implied alignments have higher SPFN than standard methods and much higher than SATé. So minimizing treelength, an explicit
  optimization criterion, improves neither alignment nor tree. This echoes Ogden & Rosenberg 2007 (*Syst Biol* 56:182, DOI 10.1080/10635150701281102), where
  POY-implied alignments were less accurate than ClustalW in 99.95% of cases.
- **Mirarab S, Nguyen N, Guo S, Wang L-S, Kim J, Warnow T (2015).** PASTA: ultra-large multiple sequence alignment for
  nucleotide and amino-acid sequences. *J Comput Biol* 22(5):377-386. DOI 10.1089/cmb.2014.0156 (verified). [FT via summary]
  The paper reports alignment (SP-error, TC) and tree FN side by side, but makes no explicit correlation claim. On RNASim 50K, SATé-II recovered 30 fully correct columns
  versus 311 for PASTA, and had tree FN of 12.6% versus 8.2%. On the parameters they varied, pairs score and tree error moved together only slightly.
- **Nguyen N-pD, Mirarab S, Kumar K, Warnow T (2015).** Ultra-large alignments using phylogeny-aware profiles (UPP).
  *Genome Biol* 16:124. DOI 10.1186/s13059-015-0688-z (verified via Europe PMC record). [FT]
  Key quote: "UPP alignments tend to have lower SP-error rates than PASTA alignments but also lower TC scores, indicating that
  these two criteria are not that well correlated. However, ML trees based on PASTA alignments ... are typically more accurate."
  So here mean(SPFN,SPFP) ranks methods the wrong way for tree accuracy, while TC ranks them the right way.
- **Nute M, Warnow T (2016).** Scaling statistical multiple sequence alignment to large datasets. *BMC Genomics*
  17(Suppl 10):764. DOI 10.1186/s12864-016-3101-8 (verified). [FT]
  Reports Modeler score, SP-score, TC and Delta-RF together. PASTA+BAli-Phy gains most on TC, and its trees are nearly always better than default PASTA, "although default PASTA
  has slightly better SP-scores ... on several of the Rose S1 replicates". MAFFT L-INS-i is worse on all alignment criteria but has the best trees on Rose M1.
  Once again, SP-score and Modeler score do not fully predict tree error. Modeler results were "nearly identical to SP-score" here.
- **Smirnov V, Warnow T (2021).** MAGUS: Multiple sequence Alignment using Graph clUStering. *Bioinformatics* 37(12):1666-1672.
  DOI 10.1093/bioinformatics/btaa992 (verified). [FT] Reports only the mean of SPFN and SPFP (both are in the supplement) and does **not** evaluate trees.
  The paper makes no statement about which criterion predicts tree error.
- **Smirnov V, Warnow T (2021).** Phylogeny estimation given sequence length heterogeneity. *Syst Biol* 70(2):268-282.
  DOI 10.1093/sysbio/syaa058 (verified). [FT] Reports SPFN and SPFP, and tree FN/FP, for UPP, PASTA and placement pipelines on fragmentary data.
  UPP+RAxML gives the best trees, and PASTA alignments have the highest SPFN/SPFP under fragmentation. There is no explicit criterion-vs-tree correlation.
  (Other Smirnov work checked and found not relevant to alignment-vs-tree criteria: Smirnov & Warnow 2020 GTM, *BMC Genomics* 21:235, DOI 10.1186/s12864-020-6605-1;
  Shen, Zaharias & Warnow 2022 MAGUS+eHMMs, *Bioinformatics* 38(4):918, DOI 10.1093/bioinformatics/btab788, which reports alignment error only.)
- **Wang L-S, Leebens-Mack J, Wall PK, Beckmann K, dePamphilis CW, Warnow T (2011).** The impact of multiple protein
  sequence alignment on phylogenetic estimation. *IEEE/ACM TCBB* 8(4):1108-1119. DOI 10.1109/TCBB.2009.68 (verified, PubMed). [Abs]
  Alignment accuracy is positively correlated with tree accuracy. The correlation is strongest when sequences are hard to align and is small when
  alignment error is generally low. The tree gain from a better alignment ranges "from quite small to substantial".
- **Nute M, Saleh E, Warnow T (2019).** Evaluating statistical multiple sequence alignment in comparison to other alignment
  methods on protein data sets. *Syst Biol* 68(3):396-411. DOI 10.1093/sysbio/syy068 (verified). [FT]
  Introduces the **expansion ratio** (estimated length / true length): below 1 means over-aligned or compressed, above 1 means under-aligned. Under-alignment shows as
  low SPFP with high SPFN, i.e. Modeler score above SP-score. BAli-Phy under-aligns on biological benchmarks but not on simulated data. On simulated proteins, Delta-RF under the hardest condition was
  28% Clustal, 20% PRANK, 9% Probalign, 7% ProbCons, 4% Muscle, 1% BAli-Phy/PRIME and 0% MAFFT-G-INS-i. PRANK (long alignments) gave poor trees here.
  The paper does not regress tree error on SPFN versus SPFP.

## Phylogeny-aware gaps, over-alignment, empirical tree-based tests

- **Löytynoja A, Goldman N (2008).** Phylogeny-aware gap placement prevents errors in sequence alignment and evolutionary
  analysis. *Science* 320(5883):1632-1635. DOI 10.1126/science.1158395 (verified). [Abs]
  Standard progressive aligners systematically infer "excess deletions and substitutions, too few insertions". They over-align by
  compressing independent insertions into shared columns. PRANK yields longer, gappier, less compressed alignments that are better for downstream analyses.
  Their downstream evaluation focused on indel and substitution inference rather than tree topology (Mem).
- **Dessimoz C, Gil M (2010).** Phylogenetic assessment of alignments reveals neglected tree signal in gaps. *Genome Biol*
  11(4):R37. DOI 10.1186/gb-2010-11-4-r37 (verified). [FT]
  Tree-based tests on real orthologs: (i) consistency-based aligners do not beat matrix-based ones; (ii) gaps carry phylogenetic signal,
  and PRANK is best at gap placement on amino acids when judged by gap-only parsimony trees; (iii) removing gaps or variable regions worsens trees;
  (iv) "alignment variability poorly predicts tree accuracy", with -r_s < 0.16 for amino acids. Bootstrap support predicts tree accuracy better
  than alignment variability (P < 0.006). Cites Wong et al. 2008 (*Science* 319:473, DOI 10.1126/science.1151532): alignment variability vs tree variability, r_s = 0.53.
- **Ogden TH, Rosenberg MS (2006).** Multiple sequence alignment accuracy and phylogenetic inference. *Syst Biol* 55(2):314-328.
  DOI 10.1080/10635150500541730 (verified). [Abs] On average, topological accuracy falls as alignment error rises. The effect is strong on pectinate trees
  and nearly absent for balanced ultrametric trees. Individual analyses can be accurate or inaccurate regardless of alignment accuracy.
- **Tan G, Muffato M, Ledergerber C, Herrero J, Goldman N, Gil M, Dessimoz C (2015).** Current methods for automated filtering
  of multiple sequence alignments frequently worsen single-gene phylogenetic inference. *Syst Biol* 64(5):778-791.
  DOI 10.1093/sysbio/syv033 (verified). [FT] Filtering (Gblocks, trimAl, etc.) on average gives worse trees and more strongly supported wrong branches.
  Removing up to roughly 20-30% of sites has little effect, and trees deteriorate rapidly beyond that. Shorter alignments that lose sites, even unreliable ones, hurt trees.
  This is about column removal, which is different from compression.
- **Talavera G, Castresana J (2007).** Improvement of phylogenies after removing divergent and ambiguously aligned blocks.
  *Syst Biol* 56(4):564-577. DOI 10.1080/10635150701472164 (verified). [Mem] Reports the opposite of Tan et al.: Gblocks-type trimming improved trees in their simulations.
- **Morrison DA (2006).** Multiple sequence alignment for phylogenetic purposes. *Aust Syst Bot* 19:479-539. DOI 10.1071/SB06020
  (verified). [Mem] A review arguing that alignment for phylogenetics should target homology (positional) rather than structural or score
  optima, and that benchmark SP/TC scores are not phylogenetic criteria. See also Morrison, Morgan & Kelchner 2015, *Aust Syst Bot* 28:46,
  DOI 10.1071/SB15001 (verified). No quantitative SPFN vs SPFP result.
- **Chowdhury B, Garai G (2017).** A review on multiple sequence alignment from the perspective of genetic algorithm.
  *Genomics* 109(5-6):419-431. DOI 10.1016/j.ygeno.2017.06.007 (verified). [Abs] Surveys GA and multi-objective GA objective functions
  (SP-type scores, gap penalties, etc.). It does not evaluate trees, so it is relevant only as a catalogue of alignment criteria. A follow-up is Chowdhury & Garai 2020,
  *Soft Computing*, a bi-objective GA, DOI 10.1007/s00500-020-04917-5 (verified, abstract unavailable).
- Not located in the time available: a Zaharias/Park paper that specifically correlates alignment criteria with tree error. (Zaharias, Smirnov & Warnow 2023 *TCBB* 20:1700,
  DOI 10.1109/TCBB.2022.3191848, is about the alignment-merging problem. Park & Warnow 2023 HMMerge, *Bioinf Adv* vbad052, DOI 10.1093/bioadv/vbad052,
  reports alignment error only.)

## Synthesis (5 lines)
1. Known: alignment error is positively but loosely correlated with ML tree error. The link is strongest on hard, divergent data and weak when alignment error is low (Wang 2011, Ogden & Rosenberg 2006).
2. Known: SPFN/SP-score alone often mis-ranks methods for tree accuracy (Nelesen 2008 FTA; Liu 2010 Prank+GT vs Opal). Mean(SPFN,SPFP) can disagree with TC, and TC tracked tree accuracy better in UPP vs PASTA (Nguyen 2015).
3. Known: PRANK-style long, under-compressed alignments sometimes give good trees with high SPFN (Liu 2010) and sometimes bad trees (Nute 2019 protein simulations, 20% Delta-RF). Removing columns usually hurts (Tan 2015, Dessimoz & Gil 2010).
4. Open: no study we found directly tests whether SPFP (over-alignment) or SPFN (under-alignment) better predicts ML tree error, e.g. by partial correlation or regression across many methods and conditions. The Warnow-group papers mostly report SPFN or the mean of the two.
5. Open: whether longer, less compressed alignments (expansion ratio of 1 or more) systematically give better ML trees. The evidence is anecdotal and conflicting; the expansion-ratio idea is in Nute 2019 but has not been tied quantitatively to RF.
