# Prior art: how much do MSA method rankings depend on the simulator?

CS581 project pilot. Literature review done 2026-10-10, timeboxed to about 40 minutes.

**How DOIs were checked.** Each DOI was resolved through the Crossref API (title, first author, year and venue matched). Most abstracts were read through Europe PMC. Where the text says "per abstract", the claim comes from the abstract only, not the full text. Anything that could not be checked against a primary source is marked **unverified**.

**Background.** The project's starting fact is that MAGUS improved on PASTA by about 28% on ROSE simulations but only about 2% on RNASim (Smirnov & Warnow 2021). That fact comes from the project framing. I did not re-derive it here. The MAGUS paper reports its RNASim numbers only in figures (Figs. 5–6).

---

## 1. Simulators and how realistic MSA benchmarks are

**Simulators**

- **ROSE.** Stoye J, Evers D, Meyer F (1998). "Rose: generating sequence families." *Bioinformatics* 14(2):157–163. DOI: 10.1093/bioinformatics/14.2.157 (verified).
  - ROSE evolves sequences down a guide tree with substitutions, insertions and deletions, and logs the true alignment.
  - Indel lengths come from a user-specified distribution. In Warnow-lab use this gives the S/M/L ("short/medium/long gap") model conditions.
  - It is the basis of the SATé 1000-taxon conditions (1000S1–S3, 1000M1–M4, 1000L1–L3), which are the datasets where MAGUS gains most over PASTA.

- **INDELible.** Fletcher W, Yang Z (2009). "INDELible: a flexible simulator of biological sequence evolution." *Mol Biol Evol* 26(8):1895–1898. DOI: 10.1093/molbev/msp098 (verified).
  - Supports several indel-length models (geometric, negative binomial, Zipf/power-law, Lavalette), heterogeneous substitution models, and codon models.
  - The Warnow lab's "Indelible 10K" (10000M2–M4) datasets come from it.

- **RNASim.** Guo S, Wang L-S, Kim J (2009). "Large-scale simulation of RNA macroevolution by an energy-dependent fitness model." arXiv:0912.2326. No journal DOI. The arXiv record exists (verified); the arXiv-assigned DOI 10.48550/arXiv.0912.2326 follows the standard pattern but I did not check it separately.
  - Uses a mutation–selection population-genetics model with fitness given by RNA secondary-structure folding energy.
  - The authors report that its sequences have more "statistical complexity" for phylogeny reconstruction than ROSE or Seq-Gen output.
  - Used in PASTA, UPP, MAGUS and TWILIGHT, with up to 1M sequences.
  - Nute et al. (2016, Table 1) quantify how it differs from ROSE:

    | Dataset | Avg gap length | Max gap length | p-distance |
    |---|---|---|---|
    | RNASim | 3.1 | 68 | lower |
    | ROSE L1 | 13.2 | — | — |
    | ROSE M1 | 9.9 | — | — |
    | ROSE S1 | 3.9 | — | — |
    | INDELible M2 | 5.6 | — | — |

    RNASim also has longer sequences and many more, shorter indels. This is a candidate explanation for the 28% vs 2% MAGUS gap.

- **AliSim.** Ly-Trong N, Naser-Khdour S, Lanfear R, Minh BQ (2022). "AliSim: a fast and versatile phylogenetic sequence simulator for the genomic era." *Mol Biol Evol* 39(5):msac092. DOI: 10.1093/molbev/msac092 (verified).
  - Built-in indel-length distributions are Geometric, Negative Binomial, Zipfian and Lavalette.
  - The **default is Zipfian with exponent 1.7**, citing Benner 1993 and Cartwright 2009.
  - The benchmark used insertion and deletion rates of 0.03 and 0.09, with truncated Zipf lengths (a=1.7, max=50).
  - Simulates 1M sequences in about 1.4 h and 1.3 GB.
  - This makes AliSim the natural tool for this project's indel-length factor.

- **AliSim-HPC.** Ly-Trong N et al. (2023). "AliSim-HPC: parallel sequence simulator for phylogenetics." *Bioinformatics* 39(9):btad540. DOI: 10.1093/bioinformatics/btad540 (verified).

- **Dawg.** Cartwright RA (2005). "DNA assembly with gaps (Dawg): simulating sequence evolution." *Bioinformatics* 21(Suppl 3):iii31–iii38. DOI: 10.1093/bioinformatics/bti1200 (verified).
  - The first widely used simulator with explicitly **power-law** indel lengths, which the author calls "biologically realistic".

- **indel-Seq-Gen.**
  - Strope CL, Scott SD, Moriyama EN (2007). "indel-Seq-Gen: a new protein family simulator incorporating domains, motifs, and indels." *Mol Biol Evol* 24(3):640–649. DOI: 10.1093/molbev/msl195 (verified).
  - v2.0: Strope CL, Abel K, Scott SD, Moriyama EN (2009). *Mol Biol Evol* 26(11):2581–2593. DOI: 10.1093/molbev/msp174 (verified).

- **SimPhy.** Mallo D, De Oliveira Martins L, Posada D (2016). "SimPhy: phylogenomic simulation of gene, locus, and species trees." *Syst Biol* 65(2):334–344. DOI: 10.1093/sysbio/syv082 (verified).
  - Models lineage, gene, and gene-by-lineage rate heterogeneity (uncorrelated relaxed clocks) and calls INDELible for sequences. See Section 5.

- **Other simulators** (all verified; not central here):
  - Seq-Gen: Rambaut & Grassly 1997, *CABIOS/Bioinformatics* 13(3):235. DOI: 10.1093/bioinformatics/13.3.235.
  - PhyloSim: Sipos et al. 2011, *BMC Bioinf* 12:104. DOI: 10.1186/1471-2105-12-104.
  - TE-sequence simulator: Hubley et al. 2022, *NAR Genom Bioinf* 4(2):lqac040. DOI: 10.1093/nargab/lqac040. It simulates neutral transposable-element families with power-law indels. MAFFT and Refiner were best at low-to-medium divergence; Refiner was uniquely good on high-divergence, fragmented data.
  - Fast indel simulation: Wygoda E et al. (2025/2026). "Efficient algorithms for simulating sequences along a phylogenetic tree." *Bioinformatics* 42(1):btaf686. DOI: 10.1093/bioinformatics/btaf686 (verified). Now the simulation engine inside SpartaABC.

**Estimating indel parameters from real data (Pupko lab)**

- **Simulation-based fitting.** Levy Karin E, Rabin A, Ashkenazy H, Shkedy D, Avram O, Cartwright RA, Pupko T (2015). "Inferring indel parameters using a simulation-based approach." *Genome Biol Evol* 7(12):3226–3238. DOI: 10.1093/gbe/evv212 (verified).
  - Fits indel parameters to an MSA by parametric simulation.
  - Shows that indel parameters "substantially vary" between mammals, bacteria and retroviruses (per abstract).

- **SpartaABC.** Levy Karin E, Shkedy D, Ashkenazy H, Cartwright RA, Pupko T (2017). "Inferring rates and length-distributions of indels using approximate Bayesian computation." *Genome Biol Evol* 9(5):1280–1294. DOI: 10.1093/gbe/evx084 (verified).
  - Web server: Ashkenazy H et al. (2017). *Nucleic Acids Res* 45(W1):W453. DOI: 10.1093/nar/gkx322 (verified).

- **Separate insertion and deletion models.** Loewenthal G, Rapoport D, Avram O, Moshe A, Wygoda E, Itzkovitch A, Israeli O, Azouri D, Cartwright RA, Mayrose I, Pupko T (2021). "A probabilistic model for indel evolution: differentiating insertions from deletions." *Mol Biol Evol* 38(12):5769–5781. DOI: 10.1093/molbev/msab266 (verified).
  - A model with separate insertion and deletion rates and length distributions fits most empirical datasets better than the symmetric model.
  - Deletion rate exceeds insertion rate in the majority of datasets. ROSE and INDELible defaults do not capture this asymmetry.

- **Which length distribution fits.** Wygoda E, Loewenthal G, Moshe A, Alburquerque M, Mayrose I, Pupko T (2024). "Statistical framework to determine indel-length distribution." *Bioinformatics* 40(2):btae043. DOI: 10.1093/bioinformatics/btae043 (verified).
  - Uses ABC model selection with alignment-bias correction.
  - **Zipf fits the vast majority of empirical alignments**, though the best model varies by alignment.
  - Notes that gap-counting is biased because aligners assume geometric (affine) lengths. This is directly relevant to the geometric-vs-Zipf factor.

**Empirical indel-length distributions**

- **Benner, Cohen & Gonnet (1993).** "Empirical and structural models for insertions and deletions in the divergent evolution of proteins." *J Mol Biol* 229(4):1065–1082. DOI: 10.1006/jmbi.1993.1105 (verified).
  - Gap probability is proportional to length^-1.7 (Zipfian) across the whole range examined, and roughly independent of evolutionary distance.

- **Qian & Goldstein (2001).** "Distribution of indel lengths." *Proteins* 45(1):102–104. DOI: 10.1002/prot.1129 (verified).
  - Structure-based alignments of distant proteins fit a four-component multi-exponential.

- **Chang & Benner (2004).** "Empirical analysis of protein insertions and deletions determining parameters for the correct placement of gaps in protein sequence alignments." *J Mol Biol* 341(2):617–631. DOI: 10.1016/j.jmb.2004.05.045 (verified).
  - Gap lengths are approximately Zipfian, proportional to L^-1.8.

- **Gu & Li (1995).** *J Mol Evol* 40:464–473. DOI: 10.1007/BF00164032 (verified).
  - Pseudogene indel sizes support a logarithmic (power-law-like) gap penalty.

- **Zhang & Gerstein (2003).** *Nucleic Acids Res* 31(18):5338. DOI: 10.1093/nar/gkg745 (verified).
  - Pseudogene-based substitution and indel patterns. Cited for context only.

- **Cartwright (2009).** "Problems and solutions for estimating indel rates and length distributions." *Mol Biol Evol* 26(2):473–480. DOI: 10.1093/molbev/msn275 (verified).
  - In primate and rodent introns, a zeta power-law fits much better than a geometric model, with exponent about 1.6–1.7 and about 12–16 indels per 100 substitutions.
  - **Using the geometric/affine model "introduces artifacts into evolutionary analysis."**

**Realism of simulated alignments**

- **Trost J, Haag J, Höhler D, Jacob L, Stamatakis A, Boussau B (2024).** "Simulations of sequence evolution: how (un)realistic they are and why." *Mol Biol Evol* 41(1):msad277. DOI: 10.1093/molbev/msad277 (verified; Crossref lists the online date as 2023).
  - Gradient-boosted-tree and CNN classifiers separate AliSim-simulated MSAs from empirical ones with balanced accuracy ≥0.93 (CNN) under every model tested, from JC up to LG+C60, with and without indels.
  - Adding indels did **not** make simulations harder to detect: CNN balanced accuracy was ≥0.97 for DNA and about 0.996 for protein.
  - Even SpartaABC-fitted indels were detectable (GBT 0.94). Only copying real gap patterns ("mimick") brought GBT down to 0.77.
  - The main signals were site-rate heterogeneity (simulated rates are too uniform) and composition.
  - The paper does not test whether method rankings change.

- **Höhler D, Pfeiffer W, Ioannidis V, Stockinger H, Stamatakis A (2022).** "RAxML Grove: an empirical phylogenetic tree database." *Bioinformatics* 38(6):1741–1742. DOI: 10.1093/bioinformatics/btab863 (verified).
  - More than 60,000 empirical trees and model parameters, published to make simulation trees, tree shapes and parameters more realistic.

- **Abadi S, Azouri D, Pupko T, Mayrose I (2019).** "Model selection may not be a mandatory step for phylogeny reconstruction." *Nat Commun* 10:934. DOI: 10.1038/s41467-019-08822-w (verified).
  - Uses empirical-like simulations. Context only: conclusions about model choice depend on how the data were simulated.

---

## 2. Clock violation and deep divergence in MSA and co-estimation studies

- **SATé-I.** Liu K, Raghavan S, Nelesen S, Linder CR, Warnow T (2009). "Rapid and accurate large-scale coestimation of sequence alignments and phylogenetic trees." *Science* 324(5934):1561–1564. DOI: 10.1126/science.1171243 (verified).
  - Introduces the ROSE 100- and 1000-taxon model conditions.
  - Model trees were generated with r8s (birth–death) and then made non-ultrametric, with branch lengths scaled by a "tree height" factor. Sequences are 1000 bp under GTR+Γ with S/M/L gaps.
  - This description comes from secondary sources that reuse the data:
    - Brown & Owen, arXiv:1708.00294.
    - Brown & Truszkowski, arXiv:1111.0379. They say they "multiplied the length of each branch by a factor chosen uniformly at random from [0.5, 2] to deviate the trees from ultrametricity. This methodology follows … Liu et al."
  - The exact SATé multiplier range is **unverified against the SATé SOM**; I could not retrieve the supplement. The clock deviation is a single fixed setting, not a varied factor.

- **SATé-II.** Liu K, Warnow TJ, Holder MT, Nelesen SM, Yu J, Stamatakis AP, Linder CR (2012). "SATé-II: very fast and accurate simultaneous estimation of multiple sequence alignments and phylogenetic trees." *Syst Biol* 61(1):90–106. DOI: 10.1093/sysbio/syr095 (verified).
  - Uses the same ROSE 1000-taxon conditions (1000M1 is hardest, 1000M4 easiest).

- **Guide-tree effects.** Nelesen S, Liu K, Zhao D, Linder CR, Warnow T (2008). "The effect of the guide tree on multiple sequence alignments and subsequent phylogenetic analyses." *Pac Symp Biocomput* 2008:25–36. DOI: 10.1142/9789812776136_0004 (verified).
  - On ROSE-type simulations, better guide trees changed SP accuracy little but substantially improved RAxML trees, especially for ProbCons and FTA.

- **Gap penalties and treelength.** Liu K, Nelesen S, Raghavan S, Linder CR, Warnow T (2009). "Barking up the wrong treelength: the impact of gap penalty on alignment and tree accuracy." *IEEE/ACM TCBB* 6(1):7–21. DOI: 10.1109/TCBB.2008.63 (verified).

- **PRANK.** Löytynoja A, Goldman N (2008). "Phylogeny-aware gap placement prevents errors in sequence alignment and evolutionary analysis." *Science* 320(5883):1632–1635. DOI: 10.1126/science.1158395 (verified).
  - Progressive aligners systematically under-infer insertions and over-infer deletions and substitutions on simulated data.
  - Which method looks best depends on whether the benchmark's indel process separates insertions from deletions. This links to Loewenthal 2021.

- **BeeTLe and treelength.** Liu K, Warnow T (2012). "Treelength optimization for phylogeny estimation." *PLoS ONE* 7(3):e33104. DOI: 10.1371/journal.pone.0033104 (verified).
  - Treelength (POY/BeeTLe) trees are less accurate than ML on good alignments across simulated and biological data.

- **PASTA.** Mirarab S, Nguyen N, Guo S, Wang L-S, Kim J, Warnow T (2015). "PASTA: ultra-large multiple sequence alignment for nucleotide and amino-acid sequences." *J Comput Biol* 22(5):377–386. DOI: 10.1089/cmb.2014.0156 (verified).
  - Conference version: RECOMB 2014, LNCS 8394. DOI: 10.1007/978-3-319-05269-4_15 (verified).
  - Introduces RNASim datasets of up to 200K sequences in the Warnow-lab benchmark suite. Guo, Wang and Kim are co-authors.

- **UPP.** Nguyen N, Mirarab S, Kumar K, Warnow T (2015). "Ultra-large alignments using phylogeny-aware profiles." *Genome Biol* 16:124. DOI: 10.1186/s13059-015-0688-z (verified).
  - Benchmarks on ROSE NT, INDELible 10K, RNASim (10K–1M) and ROSE AA. Relative performance of UPP, PASTA and MAFFT differs by simulator.

- **MAGUS.** Smirnov V, Warnow T (2021). "MAGUS: Multiple sequence Alignment using Graph clUStering." *Bioinformatics* 37(12):1666–1672. DOI: 10.1093/bioinformatics/btaa992 (verified).
  - On the ROSE 1000-sequence conditions MAGUS beats PASTA on 9 of 10 and ties on 1. The largest gap is on 1000L3, about 10–11% vs 17–18% SP-error.
  - On large CRW 16S data the gains are about 1–2% or a tie. The paper describes RNASim as a non-standard model "with positive selection".

- **Recursive MAGUS.** Smirnov V (2021). "Recursive MAGUS: scalable and accurate multiple sequence alignment." *PLoS Comput Biol* 17(10):e1008950. DOI: 10.1371/journal.pcbi.1008950 (verified; Crossref lists Smirnov as sole author).

- **PASTA+BAli-Phy.** Nute M, Warnow T (2016). "Scaling statistical multiple sequence alignment to large datasets." *BMC Genomics* 17(Suppl 10):764. DOI: 10.1186/s12864-016-3101-8 (verified).
  - Documents that RNASim, ROSE and INDELible occupy different parts of "data space" (the gap-length table in Section 1).
  - Short gaps combined with high substitution rates are the hardest conditions.

- **BAli-Phy on proteins.** Nute M, Saleh E, Warnow T (2019). "Evaluating statistical multiple sequence alignment in comparison to other alignment methods on protein data sets." *Syst Biol* 68(3):396–411. DOI: 10.1093/sysbio/syy068 (verified).
  - **This is a clear ranking reversal.** BAli-Phy has the best precision and recall on 120 simulated datasets, but consistently lower recall than many methods on 1192 biological benchmark datasets. It "systematically underaligns" on biological data, with no sign of this on simulated data.

- **Tree shape and clock regime in an MSA study.** Ogden TH, Rosenberg MS (2006). "Multiple sequence alignment accuracy and phylogenetic inference." *Syst Biol* 55(2):314–328. DOI: 10.1080/10635150500541730 (verified).
  - **The closest prior study to this project's design.** Simulated with indels on pectinate (caterpillar), balanced and random topologies, under ultrametric-equal, ultrametric-random and non-ultrametric-random branch lengths.
  - The drop in topological accuracy as alignment error rises was "much more pronounced" on pectinate trees. On balanced, ultrametric, equal-length trees, alignment error had little average effect on tree accuracy.
  - Alignment accuracy falls as a branch and its neighbouring branches get longer.
  - Only ClustalW-era aligners and small trees were used, and it did not rank modern aligners.

- **POY vs ClustalW.** Ogden TH, Rosenberg MS (2007). "Alignment and topological accuracy of the direct optimization approach via POY and traditional phylogenetics via ClustalW + PAUP*." *Syst Biol* 56(2):182–193. DOI: 10.1080/10635150701281102 (verified).
  - Same factorial design: pectinate, balanced and random trees under clocklike, non-clocklike and ultrametric branch lengths.
  - ClustalW alignments beat POY-implied alignments in 99.95% of cases.

- **Gap placement bias.** Golubchik T, Wise MJ, Easteal S, Jermiin LS (2007). "Mind the gaps: evidence of bias in estimates of multiple sequence alignments." *Mol Biol Evol* 24(11):2433–2442. DOI: 10.1093/molbev/msm176 (verified).
  - Aligners misplace gaps even in otherwise identical sequences.
  - Rankings depend on the indel configuration (overlapping vs non-overlapping deletions). MAFFT G-INS-i and DIALIGN-T were best for non-overlapping deletions.

- **Genomic alignment uncertainty.** Lunter G et al. (2008). "Uncertainty in homology inferences: assessing and improving genomic sequence alignment." *Genome Res* 18(2):298–309. DOI: 10.1101/gr.6725608 (verified).
  - More than 15% of aligned bases are wrong at human–mouse divergence.
  - Modelling the evolutionary process correctly gives only modest gains.

**Clock violation in tree-estimation studies (no MSA)**

- **Nakhleh L, Roshan U, Vawter L, Warnow T (2002).** "Estimating the deviation from a molecular clock." WABI 2002, LNCS 2452:287–299. DOI: 10.1007/3-540-45784-4_22 (verified).
  - Defines a "stretch" measure. Parsimony and ML estimates of clock deviation tend to underestimate it. This is the finding as reported in a search summary, not from my own reading of the paper.

- **Moret BME, Roshan U, Warnow T (2002).** "Sequence-length requirements for phylogenetic methods." WABI 2002, LNCS 2452:343–356. DOI: 10.1007/3-540-45784-4_26 (verified).
  - Birth–death trees with controlled deviation from ultrametricity. Required sequence length grows with the deviation and with tree height. The claim is from a search summary; Crossref has no abstract for this paper.

- **Nakhleh L, Moret BME, Roshan U, St. John K, Sun J, Warnow T (2002).** "The accuracy of fast phylogenetic methods for large datasets." *Pac Symp Biocomput* 2002:211–222. DOI: 10.1142/9789812799623_0020 (verified).
  - "Random birth-death trees, with controlled deviations from ultrametricity." The method ranking (Weighbor vs DCM-NJ+MP) depends on sequence length.

- **Felsenstein (1978)** on long-branch attraction. *Syst Zool* 27(4):401–410. DOI: 10.2307/2412923 (verified).

- **Kuhner & Felsenstein (1994).** "A simulation comparison of phylogeny algorithms under equal and unequal evolutionary rates." *Mol Biol Evol* 11(3):459–468. DOI: 10.1093/oxfordjournals.molbev.a040126 (verified).
  - **Classic evidence that method rankings change when the clock is violated.**

- **Huelsenbeck (1995).** *Syst Biol* 44(1):17–48. DOI: 10.2307/2413481 (verified).

- **Kolaczkowski & Thornton (2004).** "Performance of maximum parsimony and likelihood phylogenetics when evolution is heterogeneous." *Nature* 431:980–984. DOI: 10.1038/nature02917 (verified).
  - Under heterotachy, ML and Bayesian MCMC can become inconsistent, and parsimony outperforms them. **A simulator-dependent ranking reversal.**

- **Lopez, Casane & Philippe (2002).** "Heterotachy, an important process of protein evolution." *Mol Biol Evol* 19(1):1–7. DOI: 10.1093/oxfordjournals.molbev.a003973 (verified).

- **Mai U, Sayyari E, Mirarab S (2017).** "Minimum variance rooting…" *PLoS ONE* 12(8):e0182238. DOI: 10.1371/journal.pone.0182238 (verified).
  - Clock deviation hurts gene-tree rooting. Relevant to guide-tree and rooting steps.

---

## 3. Benchmark-dependent rankings

- **Iantorno S, Gori K, Goldman N, Gil M, Dessimoz C (2014).** "Who watches the watchmen? An appraisal of benchmarks for multiple sequence alignment." *Methods Mol Biol* 1079:59–73. DOI: 10.1007/978-1-62703-646-7_4 (verified).
  - Reviews simulation, consistency, structural and phylogenetic benchmarks. Concludes there is "no universally applicable means of benchmarking MSA."

- **Blackshields G, Wallace IM, Larkin M, Higgins DG (2006).** "Analysis and comparison of benchmarks for multiple sequence alignment." *In Silico Biol* 6(4):321–339. DOI: 10.3233/ISB-00245 (verified).
  - Across 6 benchmarks, "method performance is dependent on the input sequences." Recommends combining benchmarks to detect over-optimisation for one benchmark.

- **Edgar RC (2010).** "Quality measures for protein alignment benchmarks." *Nucleic Acids Res* 38(7):2145–2153. DOI: 10.1093/nar/gkp1196 (verified).
  - In BAliBASE, 87% of sequences have unknown structure, 20% of columns mix folds, and 30% of core-block columns have conflicting secondary structure.
  - "Calls into question their ability to determine reliable algorithm rankings."

- **Blackburne BP, Whelan S (2012).** "Measuring the distance between multiple sequence alignments." *Bioinformatics* 28(4):495–502. DOI: 10.1093/bioinformatics/btr701 (verified).
  - Introduces MetAl metrics, including ones based on gap position and indel-event position on the tree.

- **Blackburne BP, Whelan S (2013).** "Class of multiple sequence alignment algorithm affects genomic analysis." *Mol Biol Evol* 30(3):642–653. DOI: 10.1093/molbev/mss256 (verified).
  - On 200 chordate families, aligners split into similarity-based (progressive and consistency) and evolution-based (PRANK, statistical) classes.
  - Tree estimates, branch lengths and inferred positive selection depend on the class.

- **Chatzou M, Magis C, Chang J-M, Kemena C, Bussotti G, Erb I, Notredame C (2016).** "Multiple sequence alignment modeling: methods and applications." *Brief Bioinform* 17(6):1009–1023. DOI: 10.1093/bib/bbv099 (verified).
  - Includes a section on how empirical and simulated benchmarks relate and how that has shaped method development.

- **Tan G, Muffato M, Ledergerber C, Herrero J, Goldman N, Gil M, Dessimoz C (2015).** "Current methods for automated filtering of multiple sequence alignments frequently worsen single-gene phylogenetic inference." *Syst Biol* 64(5):778–791. DOI: 10.1093/sysbio/syv033 (verified).
  - Uses phylogeny-based tests on empirical and simulated data. Shows the conclusion holds across both.

- **Thompson JD, Linard B, Lecompte O, Poch O (2011).** "A comprehensive benchmark study of multiple sequence alignment methods: current challenges and future perspectives." *PLoS ONE* 6(3):e18093. DOI: 10.1371/journal.pone.0018093 (verified).
  - BAliBASE-style evaluation. Errors concentrate in local motifs, disordered regions and fragmentary sequences.

- **Nuin PAS, Wang Z, Tillier ERM (2006).** "The accuracy of several multiple sequence alignment programs for proteins." *BMC Bioinformatics* 7:471. DOI: 10.1186/1471-2105-7-471 (verified).
  - More than 30,000 Simprot simulations. Accuracy is "extremely dependent on the number of indels", while **indel size has a weaker effect**.
  - Rankings on Simprot and BAliBASE were "consistent in most cases". MAFFT L-INS-i and ProbCons were best.
  - This is evidence *against* large indel-length effects, but only for small protein alignments.

- **Dessimoz C, Gil M (2010).** "Phylogenetic assessment of alignments reveals neglected tree signal in gaps." *Genome Biol* 11:R37. DOI: 10.1186/gb-2010-11-4-r37 (verified).
  - A phylogeny-based benchmark that gives rankings different from structural benchmarks (per the Iantorno review).

- **Mirarab S, Warnow T (2011).** "FastSP: linear time calculation of alignment accuracy." *Bioinformatics* 27(23):3250–3258. DOI: 10.1093/bioinformatics/btr553 (verified).
  - The standard SP-error, TC and modeler/developer score tool. Use it for this project.

- **Linder CR, Suri R, Liu K, Warnow T (2010).** "Benchmark datasets and software for developing and testing methods for large-scale multiple sequence alignment and phylogenetic inference." *PLoS Currents* RRN1195. DOI: 10.1371/currents.RRN1195 (verified).
  - Source of the ROSE-based benchmarks.

- **Ezawa K (2016).** "Characterization of multiple sequence alignment errors using complete-likelihood score and position-shift map." *BMC Bioinformatics* 17:133. DOI: 10.1186/s12859-016-0945-5 (verified).
  - Uses power-law indel simulations. Many MAFFT errors reflect a mismatch between the aligner's score and the true likelihood, not search failure.

- **Hubley et al. (2022).** NAR Genom Bioinf; see Section 1. Rankings shift with divergence and fragmentation.

**Simulated vs empirical rankings for tree methods (context)**

- **Zhou X, Shen X-X, Hittinger CT, Rokas A (2018).** "Evaluating fast maximum likelihood-based phylogenetic programs using empirical phylogenomic data sets." *Mol Biol Evol* 35(2):486–503. DOI: 10.1093/molbev/msx302 (verified).
  - Empirical rankings (IQ-TREE > RAxML/ExaML > PhyML/FastTree by likelihood) were partly sensitive to data properties.

- **Zaharias P, Grosshauser M, Warnow T (2022).** "Re-evaluating deep neural networks for phylogeny estimation: the issue of taxon sampling." *J Comput Biol* 29(1):74–89. DOI: 10.1089/cmb.2021.0383 (verified).
  - DNN superiority on narrow simulations does not hold under other simulation regimes.

- **FAMSA2.** Gudyś A et al. (2026). "Fast and accurate multiple-protein-sequence alignment at scale with FAMSA2." *Nat Biotechnol*. DOI: 10.1038/s41587-026-03095-3 (Crossref-verified).
  - Uses simulated reference MSAs plus extHomFam. I could not access the full text, so the simulator settings and any simulated-vs-empirical ranking differences are **unverified**.

- **TWILIGHT.** Tseng Y-H, Walia S, Turakhia Y (2025). "Ultrafast and ultralarge multiple sequence alignments using TWILIGHT." *Bioinformatics* 41(Suppl 1). DOI: 10.1093/bioinformatics/btaf212 (verified).
  - Benchmarked on RNASim (aligned 1M sequences in under 30 min) and other sets.

---

## 4. Tree shape effects

- **Ogden & Rosenberg (2006, 2007).** See Section 2. These are the only MSA studies found that cross **pectinate/balanced/random topology with clock regime**. Alignment error hurt tree accuracy far more on pectinate trees.

- **Blum MGB, François O (2006).** "Which random processes describe the tree of life? A large-scale study of phylogenetic tree imbalance." *Syst Biol* 55(4):685–691. DOI: 10.1080/10635150600889625 (verified).
  - TreeBASE trees are on average much more imbalanced than Yule predicts, and are consistent with Aldous' beta-splitting model with β around −1 (between Yule, β=0, and PDA, β=−1.5).
  - The exact estimate is **unverified from the full text**. Secondary sources confirm the "more imbalanced than Yule; fits Aldous' beta-splitting model" finding.

- **Aldous DJ (2001).** "Stochastic models and descriptive statistics for phylogenetic trees, from Yule to today." *Stat Sci* 16(1):23–34. DOI: 10.1214/ss/998929474 (verified).
  - The beta-splitting family, which spans caterpillar-like to balanced shapes. Gives a single knob for the tree-shape factor.

- **Mooers AØ, Heard SB (1997).** "Inferring evolutionary process from phylogenetic tree shape." *Q Rev Biol* 72(1):31–54. DOI: 10.1086/419657 (verified).

- **Kirkpatrick M, Slatkin M (1993).** "Searching for evolutionary patterns in the shape of a phylogenetic tree." *Evolution* 47(4):1171–1181. DOI: 10.1111/j.1558-5646.1993.tb02144.x (verified; JSTOR duplicate 10.2307/2409983).
  - Colless and Sackin indices. Sackin (1972) DOI: 10.2307/2412292 (verified).

- **Stadler T (2011).** "Simulating trees with a fixed number of extant species." *Syst Biol* 60(5):676–684. DOI: 10.1093/sysbio/syr029 (verified).
  - TreeSim, for birth–death trees.

- **Duchêne DA, Duchêne S, Ho SYW (2015).** "Tree imbalance causes a bias in phylogenetic estimation of evolutionary timescales using heterochronous sequences." *Mol Ecol Resour* 15(4):785–794. DOI: 10.1111/1755-0998.12352 (verified).
  - Imbalance biases clock-based dating.

- **Höhler et al. (2022), RAxML Grove** (Section 1). An empirical source of tree shapes.

- **Divide-and-conquer methods.** In the DTM/GTM papers (Smirnov & Warnow 2020, *BMC Genomics* 21(Suppl 2):235, DOI 10.1186/s12864-020-6605-1, verified), in NJMerge (Molloy & Warnow 2019, *AMB* 14:14, DOI 10.1186/s13015-019-0151-x, verified), and in DACTAL (Nelesen et al. 2012, *Bioinformatics* 28(12):i274, DOI 10.1093/bioinformatics/bts218, verified):
  - I found **no explicit caterpillar-vs-balanced experiment**. Model trees come from SimPhy or ROSE/r8s birth–death processes.
  - This matters because centroid-edge decomposition (PASTA, MAGUS) behaves differently on caterpillar trees, where subsets span long paths, than on balanced trees. That effect appears untested.

---

## 5. How to simulate clock violation

- **Drummond AJ, Ho SYW, Phillips MJ, Rambaut A (2006).** "Relaxed phylogenetics and dating with confidence." *PLoS Biol* 4(5):e88. DOI: 10.1371/journal.pbio.0040088 (verified).
  - The uncorrelated lognormal (UCLN) and exponential relaxed clocks. The coefficient of variation of rates is the standard "clocklikeness" measure.
  - Real datasets range from near-clock to strongly non-clock. No significant autocorrelation was found in 3 large datasets.

- **Thorne JL, Kishino H, Painter IS (1998).** "Estimating the rate of evolution of the rate of molecular evolution." *Mol Biol Evol* 15(12):1647–1657. DOI: 10.1093/oxfordjournals.molbev.a025892 (verified).
  - The autocorrelated lognormal clock.
  - Also Kishino, Thorne & Bruno (2001), *Mol Biol Evol* 18(3):352. DOI: 10.1093/oxfordjournals.molbev.a003811 (verified).

- **Lepage T, Bryant D, Philippe H, Lartillot N (2007).** "A general comparison of relaxed molecular clock models." *Mol Biol Evol* 24(12):2669–2680. DOI: 10.1093/molbev/msm193 (verified).

- **Ho SYW, Duchêne S, Duchêne D (2015).** "Simulating and detecting autocorrelation of molecular evolutionary rates among lineages." *Mol Ecol Resour* 15(4):688–696. DOI: 10.1111/1755-0998.12320 (verified).
  - **NELSI** R package, which simulates branch rates under several clock models (strict, UCLN, autocorrelated) on a chronogram.
  - Rate autocorrelation is hard to detect.

- **Lanfear R, Welch JJ, Bromham L (2010).** "Watching the clock: studying variation in rates of molecular evolution between species." *Trends Ecol Evol* 25(9):495–503. DOI: 10.1016/j.tree.2010.06.007 (verified).
  - A review of real lineage rate variation.

- **SimPhy** (Mallo et al. 2016; Section 1). Lineage-specific, gene-specific and gene-by-lineage rate heterogeneity multipliers, drawn from gamma distributions (uncorrelated).

- **Mirarab S, Warnow T (2015).** "ASTRAL-II: coalescent-based species tree estimation with many hundreds of taxa and thousands of genes." *Bioinformatics* 31(12):i44–i52. DOI: 10.1093/bioinformatics/btv234 (verified).
  - Uses SimPhy "rate heterogeneity modifiers… introducing deviations from molecular clock and rate heterogeneity between genes". The parameters are in Supplementary Table S1.
  - Has become the standard Warnow/Mirarab-lab recipe for "deviation from ultrametricity".
  - The exact gamma parameters are **unverified**; they are in the supplement, which I did not open.

- **Older Warnow-lab "deviation from ultrametricity" recipe.** Multiply each branch of a birth–death (r8s) tree by an i.i.d. random factor. Used in:
  - Nakhleh et al. 2002 and Moret et al. 2002 (Section 2).
  - SATé/ROSE, reportedly with a factor in [0.5, 2] according to Brown & Truszkowski, arXiv:1111.0379.
  - This is uncorrelated, multiplicative clock violation and gives the project a simple one-parameter factor.

- **Sanderson MJ (2003).** "r8s: inferring absolute rates of molecular evolution and divergence times in the absence of a molecular clock." *Bioinformatics* 19(2):301–302. DOI: 10.1093/bioinformatics/19.2.301 (verified).
  - The tool used to generate the SATé ROSE model trees.

---

## Gap

No published study that I could find systematically varies **(a) degree of clock violation, (b) indel-length distribution (geometric vs Zipf/power-law), and (c) tree shape (caterpillar / balanced / birth–death or Aldous β)** as controlled factors, and then measures how the **ranking and relative gains** of modern large-scale aligners change. "Modern" here means MAGUS, PASTA, MAFFT (auto, L-INS-i, PartTree), FAMSA/FAMSA2, TWILIGHT, Clustal Omega, UPP/WITCH and EMMA.

The closest work has clear limits:

- **Ogden & Rosenberg (2006, 2007)** cross tree shape with clock regime, but only for ClustalW/POY on small trees. They measure the effect on downstream trees, not aligner rankings.
- **Nuin et al. (2006)** vary indel size for small protein alignments and find a weak effect.
- **Nute et al. (2016, 2019)** and the MAGUS/UPP/PASTA papers show that gains differ between ROSE, INDELible, RNASim and biological data. Each simulator, though, changes many things at once: substitution model, selection, gap-length distribution, rate and sequence length.
  - For example, RNASim's average gap length is 3.1 vs 9.9–13.2 for ROSE M1/L1, and it has a lower p-distance.
  - So the 28%-vs-2% MAGUS-over-PASTA difference has never been attributed to any one factor.
- **Trost et al. (2024)** and **Wygoda et al. (2024)** show that standard simulations are detectably unrealistic, and that empirical indels are mostly Zipfian rather than ROSE-like geometric. Neither tests consequences for aligner rankings.
- The SATé/ROSE conditions use a single fixed, uncorrelated clock-deviation recipe. They never vary it.
- In the DTM, PASTA and MAGUS papers I found no caterpillar-vs-balanced experiment, even though centroid-edge decomposition should be sensitive to tree shape.

A factorial AliSim (or INDELible) experiment would fill this gap. It would vary a branch-rate multiplier (UCLN σ or a [1/k, k] multiplier), indel-length distribution (geometric vs Zipf a≈1.7, matched on mean length and gap rate), and Aldous β or caterpillar/balanced/Yule trees, while holding diameter and indel rate fixed. Results would be scored with FastSP and reported as relative MAGUS-vs-PASTA gains and Kendall-τ of method ranks.
