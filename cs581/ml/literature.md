# Literature review: maximum likelihood (ML) phylogeny estimation, looking for tractable new-method projects

CS581 (Algorithmic Computational Genomics), UIUC. Compiled 2026-10-09.

**Verification policy.** I checked every citation below against Crossref, DataCite (for Illinois Data Bank DOIs), Europe PMC full text, or the publisher or preprint page. Where I took a number from a paper, I read it from that paper's full text, and the table is named. Items I could not verify are marked **[UNVERIFIED]**. A placeholder means I deliberately did not invent a number.

**Error metric used throughout.** FN rate (missing-branch rate) = |bipartitions(true) \ bipartitions(estimated)| / (internal edges of the true tree). When both trees are binary, FN = FP = RF rate. Delta-FN = FN(tree on estimated alignment) − FN(tree on true alignment). This is the definition used in the UPP and PASTA papers.

---

## 0. Benchmark datasets most relevant to apples-to-apples comparison

| Dataset | Origin | Where to get it | Notes |
|---|---|---|---|
| ROSE 1000-taxon (1000M1–M4, 1000S1–S3, 1000L1–L3, plus M5/S5/L5 easy conditions), 20 replicates each | SATé paper: Liu et al. 2009, *Science* 324:1561, doi:[10.1126/science.1171243](https://doi.org/10.1126/science.1171243) | Re-packaged with the MAGUS data: Illinois Data Bank [10.13012/B2IDB-2643961_V1](https://doi.org/10.13012/B2IDB-2643961_V1) (verified on DataCite; the description says it contains "10 simulated nucleotide model conditions from the SATe paper, each with 20 replicates"). Original location: https://sites.google.com/eng.ucsd.edu/datasets/alignment/sate-i. Also on Dryad [10.5061/dryad.n9r3h](https://doi.org/10.5061/dryad.n9r3h), cited in Liu et al. 2011 | GTR+indel. Reference trees are almost fully resolved: 99.4–99.6% for M1–M3, 97.3% for M4 (Smirnov & Warnow 2021, Table 3). M1 is the hardest condition and M4 the easiest. |
| RNASim (one replicate of 1M sequences; 1K, 10K and 50K subsamples) | Used in PASTA: Mirarab et al. 2015, *J Comput Biol* 22:377, doi:[10.1089/cmb.2014.0156](https://doi.org/10.1089/cmb.2014.0156) | 1K (20 subsets) and 10K (10 subsets) are in the MAGUS IDB entry above. RNASim1000 (5 subsets), Cox1-HET and 1000M1-HF, together with the analyses (baseline trees), are in IDB [10.13012/B2IDB-7008049_V1](https://doi.org/10.13012/B2IDB-7008049_V1) (Park, Zaharias, Warnow 2021). Recursive-MAGUS data: [10.13012/B2IDB-1048258_V1](https://doi.org/10.13012/B2IDB-1048258_V1) | Non-i.i.d. model with RNA-structure selection, so ML is misspecified. Reference trees are binary. |
| 1000M1-HF (half the sequences fragmented to 25% of median length) | Smirnov & Warnow 2021, *Syst Biol* 70:268, doi:[10.1093/sysbio/syaa058](https://doi.org/10.1093/sysbio/syaa058) | IDB-7008049 (above). Smirnov & Warnow's own low- and high-fragmentation versions of 1000M1–M4, RNASim and 16S.M/23S.M are on Dryad [10.5061/dryad.8pk0p2nj8](https://doi.org/10.5061/dryad.8pk0p2nj8) (link given in the paper; I could not open the Dryad landing page through the API) | Fragmentary-data benchmark. |
| Placement benchmarks: nt78, 16S.B.ALL, PEWO LTP_s128_SSU, PEWO green85 | Wedell, Cai, Warnow (SCAMPP) | IDB [10.13012/B2IDB-9257957_V1](https://doi.org/10.13012/B2IDB-9257957_V1) (verified) | Contains backbone trees, re-estimated branch lengths and chosen query sequences. |
| RNASim variable-size placement (5K–200K backbones) | Balaban et al. 2020 (APPLES) | Dryad [10.5061/dryad.78nf7dq](https://doi.org/10.5061/dryad.78nf7dq). Fragmentary queries: IDB [10.13012/B2IDB-8788479_V1](https://doi.org/10.13012/B2IDB-8788479_V1) (verified) | |

---

## 1. Fast ML heuristics and the accuracy/runtime trade-off

**Methods (all verified):**
- FastTree 2. Price, Dehal, Arkin 2010, *PLoS ONE* 5:e9490, doi:[10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490). Builds a minimum-evolution/NJ-like start, runs a bounded number of NNI/SPR rounds and ML NNIs, then stops. It does not search to convergence.
- RAxML-NG. Kozlov, Darriba, Flouri, Morel, Stamatakis 2019, *Bioinformatics* 35:4453, doi:[10.1093/bioinformatics/btz305](https://doi.org/10.1093/bioinformatics/btz305).
- IQ-TREE. Nguyen et al. 2015, *MBE* 32:268, doi:[10.1093/molbev/msu300](https://doi.org/10.1093/molbev/msu300).
- IQ-TREE 2. Minh et al. 2020, *MBE* 37:1530, doi:[10.1093/molbev/msaa015](https://doi.org/10.1093/molbev/msaa015).
- IQ-TREE 3. Wong et al., preprint EcoEvoRxiv 2025, doi:[10.32942/X2P62N](https://doi.org/10.32942/X2P62N); published in *MBE* 2026, doi:[10.1093/molbev/msag117](https://doi.org/10.1093/molbev/msag117) (Crossref-verified). Its new features are mainly about models (mixtures, concordance factors, a simulator), not tree-search speed.
- IQ-TREE `-fast`. From the IQ-TREE manual (https://iqtree.github.io/doc/Command-Reference): it "will just construct two starting trees: maximum parsimony and BIONJ, which are then optimized by NNI". This is a FastTree-like mode. The default start set is 100 parsimony trees plus a BIONJ tree.
- VeryFastTree. Piñeiro, Abuín, Pichel 2020, *Bioinformatics* 36:4658, doi:[10.1093/bioinformatics/btaa582](https://doi.org/10.1093/bioinformatics/btaa582). A parallel and vectorised re-implementation of FastTree-2.

**Comparisons:**
- Liu, Linder, Warnow 2011, *PLoS ONE* 6:e27731, doi:[10.1371/journal.pone.0027731](https://doi.org/10.1371/journal.pone.0027731). They ran RAxML, FastTree and a time-limited RAxML on the ROSE 1000-taxon conditions with six alignments each (true, SATé, MAFFT, ClustalW, QuickTree, PartTree), plus 16S datasets of up to 27,634 sequences.
  - With equal time, FastTree is much more accurate than time-limited RAxML.
  - Run to completion, RAxML gets better ML scores but only slightly better topologies. On the harder conditions the average advantage was 0.5% on true alignments, 0.5% on SATé and 1.1% on MAFFT.
  - On poor alignments (ClustalW, PartTree) FastTree is often more accurate.
  - On 16S.B.ALL with the reference alignment, FastTree had a 3.9% missing-branch rate and time-limited RAxML 13.8%. The reference tree here is itself a RAxML tree, so RAxML scores 0%.
  - Per-condition FN values for 1000M2, M3 and L1 appear only as bar charts (Fig. 1), not tables. **I could not extract exact numbers.**
- Zhou, Shen, Hittinger, Rokas 2018, *MBE* 35:486, doi:[10.1093/molbev/msx302](https://doi.org/10.1093/molbev/msx302). They tested 19 empirical phylogenomic matrices with up to 200 taxa. For concatenation, IQ-TREE got the best likelihoods with RAxML/ExaML a close second. FastTree was fastest but had lower likelihoods and more divergent topologies.
- **Park, Zaharias, Warnow 2021, *Algorithms* 14:148, doi:[10.3390/a14050148](https://doi.org/10.3390/a14050148).** This is the most useful modern head-to-head comparison: RAxML-NG 1.0.1, IQ-TREE 2.0.6 and FastTree 2.1.10, all under GTR+G on the true alignment, 16 cores and 64 GB. FN rates:

| Condition | FastTree | IQ-TREE 2 | RAxML-NG | Best DTM pipeline (GTM) |
|---|---|---|---|---|
| RNASim1000 (Table 3) | 14.9% | 15.1% | 15.1% (24 h cap hit) | 14.4% (IQ-TREE start) |
| Cox1-HET, 2341 seqs (Table 4) | 23.9% | 19.6% | 18.2% (7.1 h) | 18.7% |
| 1000M1-HF (Table 5) | 50.9% | 30.2% | 24.9% (24 h cap hit) | ~28.4–28.5% |
| RNASim10K (text) | 10.8% (2 h) | 10.9% (~75 h) | 12.3% (168 h cap hit) | 10.1% (3.5 h) |
| RNASim50K (text) | 8.0% (8.3 h) | failed (memory) | 100% error (168 h) | 7.5% (14.1 h) |

  Main points from this paper:
  - FastTree is competitive on "easy" RNASim data but loses badly on fragmentary data (1000M1-HF) and on heterotachy (Cox1-HET).
  - RAxML-NG and IQ-TREE are not reliably usable at 50K sequences under these resource limits.

**Open problem.** Can we get FastTree-like runtime with RAxML-NG-like accuracy on hard conditions (fragmentary, heterotachous) at 1K–50K taxa?

**Student-scale idea.** Run FastTree (or IQ-TREE `-fast`), then a *short* RAxML-NG search (`--tree FT.nwk`, limited SPR radius or a time budget), on RNASim1000, Cox1-HET and 1000M1-HF from IDB-7008049. Use the stored trees in that entry as baselines.
- Criteria: FN, final log-likelihood, wall-clock time, peak memory.
- Already done? Partly. Park et al. used IQ-TREE start trees inside DTM pipelines, not as "FastTree → RAxML-NG polish". I did not find a published systematic "FT start + budgeted RAxML-NG" curve on these exact datasets. The idea is not deeply novel, though, because practitioners routinely use FastTree start trees.

---

## 2. Starting-tree effects and search difficulty

- Morel, Barbera, Czech, Bettisworth, Hübner, Lutteropp, Serdari, Kostaki, Mamais, Kozlov, Pavlopoulos, Paraskevis, Stamatakis 2021, "Phylogenetic Analysis of SARS-CoV-2 Data Is Difficult", *MBE* 38:1777, doi:[10.1093/molbev/msaa314](https://doi.org/10.1093/molbev/msaa314). They ran RAxML-NG from 50 random and 50 parsimony starts. "Tree searches initiated on parsimony starting trees yielded phylogenies with better log-likelihood scores (consistently >400 log-likelihood units)." Average pairwise relative RF between the ML trees was about 0.78, which shows a rugged likelihood surface for many-taxa, few-sites data.
- Haag, Höhler, Bettisworth, Stamatakis 2022, "From Easy to Hopeless — Predicting the Difficulty of Phylogenetic Analyses" (Pythia), *MBE* 39:msac254, doi:[10.1093/molbev/msac254](https://doi.org/10.1093/molbev/msac254).
  - A random-forest regressor predicts difficulty in [0,1] from MSA features and parsimony trees. It was trained on 3,250 TreeBASE MSAs and reports MAPE 2.9%.
  - Computing the prediction takes about 1/5 the time of one ML inference.
  - Code: https://github.com/tschuelia/PyPythia
- Togkousidis, Kozlov, Haag, Höhler, Stamatakis 2023, "Adaptive RAxML-NG", *MBE* 40:msad227, doi:[10.1093/molbev/msad227](https://doi.org/10.1093/molbev/msad227). Tested on 9,515 empirical and 5,000 simulated MSAs. Speedups were more than 10× on easy and hard datasets (53% of MSAs), and about 94% of the trees were statistically indistinguishable from standard RAxML-NG trees. RAxML-NG v2 now defaults to `--tree auto`; v1.2 used 10 random plus 10 parsimony starts (RAxML-NG wiki, Input-data page).
- In Park et al. 2021 (Section 1), switching the DTM starting tree from FastTree to IQ-TREE changed FN by **14.1 points** on 1000M1-HF (42.4% vs 28.4%). The difference was ≤0.3 points on RNASim1000 and Cox1-HET. So starting trees matter most on fragmentary data.

**Student idea.** Run a starting-tree study for RAxML-NG and IQ-TREE on ROSE 1000M1–M4 and RNASim1K, comparing random, parsimony, BIONJ, FastTree and GTM-pipeline starts.
- Criteria: FN, log-likelihood, time to reach X% of the best log-likelihood.
- Novelty: low to medium. Starting-tree effects are known (Morel 2021; Pythia/adaptive work), but the Warnow-lab simulated benchmarks with true trees have not been used for this, as far as I found. The project is very tractable.

---

## 3. Divide-and-conquer (DCM and disjoint tree mergers) for ML

**History (all verified):**
- DCM: Huson, Nettles, Warnow 1999, *J Comput Biol* 6:369, doi:[10.1089/106652799318337](https://doi.org/10.1089/106652799318337).
- Rec-I-DCM3: Roshan, Moret, Warnow, Williams 2004, *Proc. IEEE CSB 2004*, doi:[10.1109/CSB.2004.1332422](https://doi.org/10.1109/CSB.2004.1332422).
- Constrained-INC: Zhang, Rao, Warnow 2019, *Algorithms Mol Biol* 14:2, doi:[10.1186/s13015-019-0136-9](https://doi.org/10.1186/s13015-019-0136-9).
- INC-ML: Le, Sy, Molloy, Zhang, Rao, Warnow 2021, **"Using Constrained-INC for Large-Scale Gene Tree and Species Tree Estimation"**, *IEEE/ACM TCBB* 18:2, doi:[10.1109/TCBB.2020.2990867](https://doi.org/10.1109/TCBB.2020.2990867). The title in the request ("Using INC within divide-and-conquer…") is not the published title. Park et al. summarise the result: INC-ML "was much less accurate than RAxML and often not as accurate as FastTree 2".
- NJMerge: Molloy & Warnow, RECOMB-CG 2018 LNCS, doi:[10.1007/978-3-030-00834-5_15](https://doi.org/10.1007/978-3-030-00834-5_15). Journal version: *Algorithms Mol Biol* 14:14 (2019), doi:[10.1186/s13015-019-0151-x](https://doi.org/10.1186/s13015-019-0151-x).
- TreeMerge: Molloy & Warnow 2019, *Bioinformatics* 35:i417, doi:[10.1093/bioinformatics/btz344](https://doi.org/10.1093/bioinformatics/btz344).
- GTM: Smirnov & Warnow 2020, "Unblended disjoint tree merging using GTM improves species tree estimation", *BMC Genomics* 21(Suppl 2):235, doi:[10.1186/s12864-020-6605-1](https://doi.org/10.1186/s12864-020-6605-1). Polynomial time; adds edges between subset trees so as to minimise RF to a guide tree. In the species-tree setting it generally matched or beat NJMerge and TreeMerge and was much faster.
- RF-supertree merger: Yu, Le, Christensen, Molloy, Warnow, bioRxiv 2020, doi:[10.1101/2020.05.16.099895](https://doi.org/10.1101/2020.05.16.099895).

**Does GTM+ML match RAxML on large data?** (Park et al. 2021; numbers in Section 1)
- The GTM pipeline uses an IQ-TREE or FastTree start, centroid decomposition to at most 500 taxa per subset, and IQ-TREE subset trees.
- Accuracy: it matches or beats IQ-TREE and FastTree everywhere. It beats RAxML-NG on RNASim1000, 10K and 50K, but is 3.5 points worse than RAxML-NG on 1000M1-HF and 0.5 points worse on Cox1-HET.
- Runtime: 1.2–3.1 h on the 1K–2.3K datasets, versus RAxML-NG's 7.1 h to more than 24 h.
- The authors' own conclusion: DTM pipelines "do not reliably match or improve the accuracy of RAxML-NG" when RAxML-NG is feasible. **This is the gap a project can attack.**

**Student ideas.**
1. *GTM + RAxML-NG polish.* Use the GTM tree as `--tree` for a short RAxML-NG search, or do SPR moves restricted to edges near subset boundaries. That targets the "no blending" weakness, since GTM can only add edges between subsets.
2. Use RAxML-NG instead of IQ-TREE for the subset trees, with subsets of 250–1000 taxa.

Data and baselines: IDB-7008049, which includes the analyses. Criteria: FN, log-likelihood, time, memory. I found no paper doing (1). It looks **novel, tractable and low-risk**, because GTM, RAxML-NG and the data are all public. The RNASim10K/50K subsets might not be in IDB-7008049 (its description lists only RNASim1000 and Cox1-Het). They can be resampled from the 1M RNASim set or taken from the MAGUS IDB 10K subsets.

---

## 4. ML on estimated or gappy alignments: trimming, masking, co-estimation

**Trimming and filtering:**
- Gblocks: Castresana 2000, *MBE* 17:540, doi:[10.1093/oxfordjournals.molbev.a026334](https://doi.org/10.1093/oxfordjournals.molbev.a026334).
- trimAl: Capella-Gutiérrez et al. 2009, *Bioinformatics* 25:1972, doi:[10.1093/bioinformatics/btp348](https://doi.org/10.1093/bioinformatics/btp348).
- BMGE: Criscuolo & Gribaldo 2010, *BMC Evol Biol* 10:210, doi:[10.1186/1471-2148-10-210](https://doi.org/10.1186/1471-2148-10-210).
- GUIDANCE2: Sela, Ashkenazy, Katoh, Pupko 2015, *NAR* 43:W7, doi:[10.1093/nar/gkv318](https://doi.org/10.1093/nar/gkv318).
- TAPER: Zhang, Zhao, Braun, Mirarab 2021, *Methods Ecol Evol*, doi:10.1111/2041-210X.13696. **[Partially verified: the DOI appeared only through Crossref review records.]**
- **Tan et al. 2015**, *Syst Biol* 64:778, doi:[10.1093/sysbio/syv033](https://doi.org/10.1093/sysbio/syv033): "trees obtained from filtered MSAs are on average worse than those obtained from unfiltered MSAs". Light filtering (≤20% of positions) has little impact.
- **ClipKIT**: Steenwyk, Buida, Li, Shen, Rokas 2020, *PLoS Biol* 18:e3001007, doi:[10.1371/journal.pbio.3001007](https://doi.org/10.1371/journal.pbio.3001007). It keeps parsimony-informative sites instead of removing "bad" ones, was tested on about 140K alignments, and outperformed the other trimmers.

**Co-estimation and alignment:**
- SATé: Liu et al. 2009 (above).
- SATé-II: Liu et al. 2012, *Syst Biol* 61:90, doi:[10.1093/sysbio/syr095](https://doi.org/10.1093/sysbio/syr095).
- PASTA: Mirarab et al. 2015 (above).
- MAGUS: Smirnov & Warnow 2021, *Bioinformatics* 37:1666, doi:[10.1093/bioinformatics/btaa992](https://doi.org/10.1093/bioinformatics/btaa992).
- UPP: Nguyen, Mirarab, Kumar, Warnow 2015, *Genome Biol* 16:124, doi:[10.1186/s13059-015-0688-z](https://doi.org/10.1186/s13059-015-0688-z).

**Published tree-error numbers I found.**
- **MAGUS main text reports only alignment error (SPFN/SPFP), not tree error.** I could not find MAGUS-alignment FN rates for 1000M2, M3, L1 or RNASim1K in the main text. They may be in its supplement **[NOT FOUND]**.
- UPP paper, Table 2 (full-length data; FastTree on each alignment; average Delta-FN):

| Collection | UPP | PASTA | MAFFT | Muscle | Clustal |
|---|---|---|---|---|---|
| ROSE NT (1000-taxon) | 1.3 | 1.3 | 5.8 | 8.4 | 24.3 |
| RNASim 10K | 0.8 | 0.4 | 3.5 | 7.3 | 10.4 |

- UPP paper, Table 3 (fragmentary data; average Delta-FN):

| Collection | UPP | PASTA | MAFFT |
|---|---|---|---|
| ROSE NT | 1.9 | 25.2 | 18.0 |
| RNASim 10K | 3.1 | 21.9 | 6.2 |

- So on full-length ROSE and RNASim, PASTA-alignment trees are within about 0.4–1.3 FN points of true-alignment trees. Little headroom is left for "better alignment → better ML tree" on full-length data. The headroom is in fragmentary data and in the ML search itself.

**Student idea.** Gap-aware column weighting or masking before ML: drop columns with more than X% gaps, or keep only parsimony-informative sites. Apply it to MAGUS or PASTA alignments of ROSE 1000M1–M4 and RNASim1K, then run FastTree and RAxML-NG.
- Criteria: FN, Delta-FN versus the true alignment, runtime.
- Novelty: low. Tan 2015 and ClipKIT already show filtering rarely helps. Still, this is a clean negative or positive replication on the Warnow-lab benchmarks, and very tractable.

---

## 5. ML with fragmentary sequences

- Sayyari, Whitfield, Mirarab 2017, *MBE* 34:3279, doi:[10.1093/molbev/msx261](https://doi.org/10.1093/molbev/msx261). Fragmentary sequences hurt gene-tree and species-tree accuracy. FastTree was poor even on true alignments (101-taxon simulations, as summarised by Smirnov & Warnow 2021). They also proposed filtering fragments.
- Smirnov & Warnow 2021, *Syst Biol* 70:268, doi:[10.1093/sysbio/syaa058](https://doi.org/10.1093/sysbio/syaa058). They compared "align all, then ML" with "tree on full-length sequences, then place the fragments".
  - Best pipeline: UPP alignment followed by RAxML. Placement pipelines were clearly worse, and FastTree was poor on fragmentary alignments.
  - FN rates from Tables 4 and 5 (main text):

| Condition | Method | 1000M1 | 1000M2 | 1000M3 | 1000M4 | RNASim (1K) |
|---|---|---|---|---|---|---|
| Low fragmentation (25% fragments, ~50% length) | PASTA-RAxML | 24.6 | 18.1 | 9.5 | 6.1 | 18.6 |
| | UPP(R)-RAxML | 15.7 | 12.8 | 9.4 | 6.1 | 18.5 |
| | UPP(R)-pplacer | 21.5 | 18.3 | 15.0 | 11.2 | 24.3 |
| High fragmentation (50% fragments, ~25% length) | PASTA-RAxML | 76.5 | 61.6 | 35.5 | 16.4 | 43.6 |
| | UPP(R)-RAxML | 37.0 | 30.4 | 23.7 | 16.7 | 37.7 |
| | UPP(R)-pplacer | 48.8 | 43.7 | 38.0 | 32.0 | 50.7 |

  - FastTree variants are in their Dryad supplement (Tables 8–11) **[not extracted]**.

**Student idea.** Two-phase ML:
1. Run RAxML-NG on the full-length sequences only.
2. Use that tree as a *starting tree* (not a constraint) for RAxML-NG on all sequences. Alternatively, use it as a GTM guide tree, with fragments placed by EPA-ng/pplacer before a final short polish.

Evaluate on 1000M1-HF (IDB-7008049) and on the Smirnov–Warnow low/high fragmentation data, where the baselines are the table above. Novelty: I found no paper that does "placement-initialised ML search" on these data. **Medium novelty, tractable.** Risk: the gain may be small compared with UPP+RAxML.

---

## 6. ML-based phylogenetic placement

**Methods (all verified):**
- pplacer: Matsen, Kodner, Armbrust 2010, *BMC Bioinformatics* 11:538, doi:[10.1186/1471-2105-11-538](https://doi.org/10.1186/1471-2105-11-538).
- EPA-ng: Barbera et al. 2019, *Syst Biol* 68:365, doi:[10.1093/sysbio/syy054](https://doi.org/10.1093/sysbio/syy054).
- APPLES: Balaban, Sarmashghi, Mirarab 2020, *Syst Biol* 69:566, doi:[10.1093/sysbio/syz063](https://doi.org/10.1093/sysbio/syz063).
- APPLES-2: Balaban, Jiang, Roush, Zhu, Mirarab 2022, *Mol Ecol Resour* 22:1213, doi:[10.1111/1755-0998.13527](https://doi.org/10.1111/1755-0998.13527).
- SCAMPP: Wedell, Cai, Warnow 2023, *IEEE/ACM TCBB* 20:1417, doi:[10.1109/TCBB.2022.3170386](https://doi.org/10.1109/TCBB.2022.3170386). Each query picks a small placement subtree, for example B = 2000 leaves, and pplacer or EPA-ng then places it within that subtree. This scales to 200K-leaf backbones with accuracy at or above APPLES-2.
- SCAMPP+FastTree: Chu & Warnow 2023, *Bioinformatics Advances* 3:vbad008, doi:[10.1093/bioadv/vbad008](https://doi.org/10.1093/bioadv/vbad008). Uses FastTree instead of RAxML for the backbone numeric parameters.
- BSCAMPP (batch): Wedell, Shen, Warnow 2025, *IEEE TCBB*, doi:[10.1109/TCBBIO.2025.3562281](https://doi.org/10.1109/TCBBIO.2025.3562281). Preprint: doi:[10.1101/2022.10.26.513936](https://doi.org/10.1101/2022.10.26.513936).

**Criterion.** Delta error (leave-one-out): the increase in FN when the query is added to the backbone, compared with the true tree restricted to backbone + query. Runtime and peak memory are also reported.

**Datasets.** nt78, ROSE 1000M1 and M5, RNASim, 16S.B.ALL, PEWO LTP_s128_SSU and green85 (IDB-9257957), and RNASim 5K–200K (APPLES Dryad).

**Student idea.** Alternative placement-subtree selection for SCAMPP, for example by k-mer/Jaccard similarity or by alignment-free distance instead of Hamming distance to the backbone leaves, or several small subtrees with likelihood voting. Novelty: medium. Many SCAMPP variants already exist, so check the 2026 *Methods Mol Biol* chapter "Phylogenetic Placement Using SCAMPP and Batch-SCAMPP" (doi:[10.1007/978-1-0716-4836-0_3](https://doi.org/10.1007/978-1-0716-4836-0_3)) first. The project is very tractable because the code and data are public. This is placement rather than "ML tree estimation" proper.

---

## 7. Deep learning and machine learning for ML phylogenetics

- Suvorov, Hochuli, Schrider 2020, *Syst Biol* 69:221, doi:[10.1093/sysbio/syz060](https://doi.org/10.1093/sysbio/syz060). CNN classification of **quartets** (4 taxa).
- Zou, Zhang, Guan, Zhang 2020, *MBE* 37:1495, doi:[10.1093/molbev/msz307](https://doi.org/10.1093/molbev/msz307). Residual networks for **quartets**.
- Azouri, Abadi, Mansour, Mayrose, Pupko 2021, *Nat Commun* 12:1983, doi:[10.1038/s41467-021-22073-8](https://doi.org/10.1038/s41467-021-22073-8). Machine learning ranks SPR moves to guide the search.
- Azouri, Granit, Alburquerque, Mansour, Pupko, Mayrose 2024, *MBE* 41:msae105, doi:[10.1093/molbev/msae105](https://doi.org/10.1093/molbev/msae105). Reinforcement learning ("tree reconstruction game"), tested on "dozens of sequences". It is about 3× faster than standard software on 15-sequence × 18 kbp data. Code: https://github.com/michaelalb/ThePhylogeneticGame
- Smith & Hahn 2023, *Bioinformatics* 39:btad543, doi:[10.1093/bioinformatics/btad543](https://doi.org/10.1093/bioinformatics/btad543). A GAN for phylogenetic inference (small trees).
- PhyloGFN: Zhou et al., ICLR 2024 (https://proceedings.iclr.cc/paper_files/paper/2024/hash/9837dc00ff67d176373268ed48042d49-Abstract-Conference.html; arXiv:2310.08774). Parsimony and Bayesian sampling on standard DS1–DS8 benchmarks (tens of taxa); this is not ML point estimation.
- Phyloformer: Nesterenko, Blassel, Veber, Boussau, Jacob 2025, *MBE* 42:msaf051, doi:[10.1093/molbev/msaf051](https://doi.org/10.1093/molbev/msaf051).
  - It learns distances or trees from the MSA. Topologically it "falls behind maximum likelihood approaches, especially as the number of sequences increases".
  - It used up to 24.2 GB of GPU RAM at 180 sequences, so it cannot handle larger trees on 32 GB GPUs.
  - Code and data: https://github.com/lucanest/Phyloformer
- Pythia and Adaptive RAxML-NG are described in Section 2.

**Honest assessment.** Deep-learning tree inference has so far been shown on quartets up to roughly 100–200 taxa. None of it is competitive with RAxML-NG or IQ-TREE topologically on 1000+ taxa. The realistic way to use it is as a *component*, for example Phyloformer on ≤100-taxon DTM subsets, or a learned difficulty or move-ranking model.

**Student idea.** Use Phyloformer (pretrained, GPU) to build subset trees inside a GTM pipeline on RNASim1000. Novelty: high. Risk: also high, because the pretrained models target protein or specific nucleotide models, a GPU is needed, and accuracy may be poor.

---

## Candidate project ideas (ranked by 4-week tractability)

| Rank | Idea | Tractability (4 wk) | Risk | Novelty / already done? | Data and baselines available |
|---|---|---|---|---|---|
| 1 | **GTM pipeline + short RAxML-NG polish** (GTM tree as `--tree`, budgeted SPR), aimed at closing the gap to RAxML-NG on 1000M1-HF (28.4% vs 24.9%) and Cox1-HET (18.7% vs 18.2%) while staying faster | High | Low–Med | Not found in the literature. Park et al. 2021 state the gap explicitly. | IDB-7008049 (data and analyses); RNASim from the MAGUS IDB |
| 2 | **Starting-tree study / FastTree→RAxML-NG time-budget curve** on ROSE 1000M1–M4, RNASim1K and 1000M1-HF | Very high | Low | Partly done (Morel 2021 empirical; adaptive RAxML-NG; Park DTM starts). Not done systematically on these simulated benchmarks. | MAGUS IDB, IDB-7008049 |
| 3 | **Fragment-aware two-phase ML** (full-length tree → starting/guide tree → place fragments → polish) | Med–High | Med | Smirnov & Warnow 2021 compared placement with align+ML, but not placement-initialised ML search. Appears novel. | Smirnov–Warnow Dryad; IDB-7008049 1000M1-HF; FN baselines in Section 5 |
| 4 | **Gap/fragment column masking before ML** on MAGUS/PASTA alignments | Very high | Low (but likely a null result) | Largely done (Tan 2015; ClipKIT 2020). New only on these benchmarks. | MAGUS IDB (includes alignments) |
| 5 | **Pythia-guided method choice** (FastTree vs IQ-TREE `-fast` vs RAxML-NG) by predicted difficulty | High | Med | Adaptive RAxML-NG (2023) already does the within-RAxML version. Cross-tool selection appears untested on these data. | ROSE/RNASim; PyPythia |
| 6 | **SCAMPP subtree-selection variant** | High | Low–Med | Many SCAMPP/BSCAMPP variants exist, so check before starting. | IDB-9257957; APPLES Dryad |
| 7 | **Phyloformer subset trees inside GTM** | Low–Med | High (GPU, model mismatch, memory) | Not found. | RNASim1000; Phyloformer GitHub |

**Flags.**
- The following figures are in bar charts or supplements, which I could not extract:
  - per-condition FN for RAxML and FastTree on the *true* alignment of 1000M2, M3 and L1 (Liu et al. 2011 Fig. 1; SATé 2009)
  - MAGUS-alignment tree error
  - Smirnov–Warnow FastTree tables (Dryad)
- Sayyari 2017, PASTA and SATé full texts were not open-access to my tools. Claims about them come from the abstracts or from secondary summaries in Smirnov & Warnow 2021 and Park et al. 2021.
