AI-assisted (Claude), literature notes for CS581 project

# Novelty audit: backbone-vote filtering and re-weighting of MAGUS/GCM edges ("gcmvote")

Audit date: 2026-10-11. This was a literature-only audit using web search, PubMed E-utilities, Crossref, and Semantic Scholar
citation lists, plus a read of the local MAGUS 0.2 source (`code/MAGUS/magus/align/merge/graph_build/graph_builder.py`).

Ideas under test:
- **(1) Support filter:** keep a GCM edge only if at least k backbones vote for it.
- **(2) Posterior from a mixture:** fit an unsupervised zero-truncated binomial or beta-binomial mixture on (votes k, exposure n)
  and take the posterior that an edge is true.
- **(3) Strength-weighted votes:** weight each backbone's vote by its fraction of realised residue pairs.
- **(4) Per-backbone reliability:** estimate each backbone's reliability Dawid–Skene style (EM over edges).
- **(5) Effective number of votes:** count overlap-aware effective independent votes.
- **(6) Dataset-level gate:** a reference-free rule that decides whether to filter at all.

Closeness codes: **I** = identical, **C** = closely related, **L** = loosely related, **–** = no relation.

## Key facts about the baseline (verified)

- **MAGUS (Smirnov & Warnow 2021).** An edge exists if "at least one pair of letters" from the two columns shares a backbone
  column, and "the weight of each such edge is the number of such pairs". Backbones are "not constrained to be consistent
  ... with each other". MCL is described as "culling the most weakly supported edges". This is the only edge pruning, and it
  has no explicit cutoff. Figure 2 varies only the number and size of the backbones; defaults are 10 backbones of 200
  sequences.
- **MAGUS 0.2 code.** `addAlignmentFileToGraph` adds `avalue * bvalue` per backbone column, which is the raw pair count. There
  is no vote counting and no threshold. `--graphbuildrestrict` only drops edges that conflict within a subset.
- **MWT-AM (Zaharias, Smirnov & Warnow, TCBB 2023).** The weight is defined as "if x and y are two letters that are aligned in q
  backbone alignments, then this particular pair (x, y) contributes q to the total weight". In the future-work section:
  "future research could explore modifying how GCM defines the weights on the pairs of columns". The GCM variants they test
  (FM, MWT-greedy/search, opt) change clustering and trace, not weights.
- **Later MAGUS-family work.** Recursive MAGUS, MAGUS+eHMMs, EMMA and UPP2 leave GCM's backbone weights unchanged. On
  2026-10-11 we scanned Semantic Scholar citations of MWT-AM (4 citing papers) and of MAGUS (69 citing papers; titles from
  2023 on). None modifies GCM weights: TWILIGHT, FAMSA2, learnMSA2 and MuSAlS are different methods.

## Paper-by-paper

### A. MAGUS / GCM family

| Paper | What it does | (1) | (2) | (3) | (4) | (5) | (6) |
|---|---|---|---|---|---|---|---|
| Smirnov & Warnow, MAGUS, *Bioinformatics* 37:1666 (2021), doi:10.1093/bioinformatics/btaa992 | Baseline GCM: count-weighted alignment graph, then MCL | L (MCL culls weak edges implicitly) | – | L (raw count already measures strength) | – | – | – |
| Zaharias, Smirnov & Warnow, MWT-AM, *IEEE/ACM TCBB* 20:1700 (2023), doi:10.1109/TCBB.2022.3191848 | Poses MWT-AM with backbone-count weights; GCM as heuristic; changing weights is listed as future work | – | – | – | – | – | – |
| Smirnov & Warnow, Recursive MAGUS, *PLoS CB* 17:e1008950 (2021), doi:10.1371/journal.pcbi.1008950 | Recursion on large subsets; GCM unchanged | – | – | – | – | – | – |
| Shen, Park & Warnow, WITCH, *JCB* 29:782 (2022), doi:10.1089/cmb.2021.0585 | **Weighted GCM.** Each query's k extended alignments get weight w = P(H_i given q) from adjusted bitscores; "The weight of the edge is the weight of the HMM". Ablation: weighted beats unweighted. The prior for HMMs on overlapping sequence sets is s_i / S#. Future work: "consensus alignments, potentially equipped with statistically defined weights" | – | L (posterior over sources, not over edges) | **C** (per-source weighted votes inside GCM) | C/L (per-source weights, but from model likelihood, not from agreement or EM) | L (overlap-aware prior) | – |
| Liu & Warnow, WITCH-NG, *Bioinf Adv* 3:vbad024 (2023), doi:10.1093/bioadv/vbad024 | Faster WITCH; same weighted GCM ("weighted by the support defined by these extended alignments") | – | – | C | L | – | – |
| Park & Warnow, HMMerge, *Bioinf Adv* 3:vbad052 (2023), doi:10.1093/bioadv/vbad052 | Combines the HMMs with WITCH weights into one HMM | – | – | L | L | – | – |
| Shen et al., MAGUS+eHMMs (*Bioinformatics* 2022, doi:10.1093/bioinformatics/btab788); EMMA (*AMB* 18:21, 2023, doi:10.1186/s13015-023-00247-x) | Add sequences to a backbone; GCM weights untouched | – | – | – | – | – | – |

### B. Consensus and ensemble alignment

| Paper | What it does | (1) | (2) | (3) | (4) | (5) | (6) |
|---|---|---|---|---|---|---|---|
| Bucka-Lassen, Caprani & Hein, ComAlign, *Bioinformatics* 15:122 (1999), doi:10.1093/bioinformatics/15.2.122 | Combines good sub-alignments from several MSAs by dynamic programming | L | – | – | – | – | – |
| Prasad et al., Consensus, *Bioinformatics* 19:1682 (2003), doi:10.1093/bioinformatics/btg211 | Consensus of 5 aligners (pairwise target–template); removes segments likely to be wrong | C/L | – | – | – | – | – |
| Wallace et al., M-Coffee, *NAR* 34:1692 (2006), doi:10.1093/nar/gkl091 | Library built from 8 aligners' MSAs; "the original T-Coffee weight is multiplied by the method weight". Tests 4 method-weight schemes: VarCov (inverse covariance), ACL and THG (tree-based), ACC (accuracy). Results: weights "fail to significantly outperform a simple combination ... (No weights)" and "do not appear to properly address the problem of method redundancy". Duplicated methods lower accuracy. "Consistency is only useful as an accuracy indicator when methods are unlikely to commit exactly the same error" | L (no threshold) | – | **C** | **C** (per-method weights, supervised or heuristic; negative result) | **C** (redundancy weights that discount correlated voters) | – |
| Lassmann & Sonnhammer, MUMSA, *NAR* 33:7120 (2005), doi:10.1093/nar/gki1020 | Overlap score between alternative alignments. **AOS** is the mean pairwise overlap, a per-case reference-free difficulty measure (r = 0.81 vs accuracy on BAliBASE; a 0.8 cutoff flags 96–100% of cases below 80% accuracy). **MOS** is a per-alignment score that weights residue pairs by "the number of m − 1 alignments that contain σ". No correction for non-independent aligners | L | – | – | C/L (unsupervised per-alignment reliability from agreement; one step, no EM) | – | **C** |
| Lassmann & Sonnhammer, *BMC Bioinf* 8(S5):S9 (2007), doi:10.1186/1471-2105-8-S5-S9 | Consensus built from pairs of aligned residues (POARs) with stringency f: "f = 3 requires all POARs in the final alignment occur in at least three input alignments". Fraction correct rises 25–100% | **C, near-identical rule** (vote threshold on residue pairs before merging) | – | – | – | – | – |
| Collingridge & Kelly, MergeAlign, *BMC Bioinf* 13:117 (2012), doi:10.1186/1471-2105-13-117 | DAG of columns; edge weight = "number of MSAs containing that transition"; inputs unweighted. Column score = proportion of MSAs, fitted to precision as f(x) = x^m with m = 0.124 (supervised calibration). No threshold | L | L (supervised calibration of support to precision) | L | – | – | – |
| Modzelewski & Dojer, MSARC, *AMB* 9:12 (2014), doi:10.1186/1748-7188-9-12 | Clusters residues on a weighted residue graph to form columns (a structural analogue of GCM) | L | – | – | – | – | – |
| Muller et al., AQUA, *Bioinformatics* 26:263 (2010), doi:10.1093/bioinformatics/btp651 | Runs MUSCLE and MAFFT, refines with RASCAL, then picks per family by NorMD | – | – | – | – | – | C/L (per-dataset choice by a reference-free score) |
| Kececioglu & DeBlasio, Facet / parameter advising, *JCB* 20:259 (2013), doi:10.1089/cmb.2013.0007; DeBlasio & Kececioglu, "Ensemble multiple sequence alignment via advising", ACM-BCB 2015, pp. 452–461, doi:10.1145/2808719.2808766 (also book chapter 2017, doi:10.1007/978-3-319-64918-4_7) | Picks the aligner or parameters per input with a learned reference-free accuracy estimator | – | – | – | – | – | **C** |
| Edgar, Muscle5, *Nat Commun* 13:6968 (2022), doi:10.1038/s41467-022-34630-w | Ensembles from HMM perturbation and guide-tree permutation. CC = fraction of replicates containing a column. The manual says CC "over-estimates if you try to interpret it as a probability", that dispersion (mean column-difference between replicates) "predicts error rate", and that `-maxcc` picks a replicate | C (support fraction) | L (motivates calibration) | – | – | – | **C** |

### C. Column reliability and probabilistic consistency

| Paper | What it does | (1) | (2) | (3) | (4) | (5) | (6) |
|---|---|---|---|---|---|---|---|
| Notredame et al., T-Coffee, *JMB* 302:205 (2000), doi:10.1006/jmbi.2000.4042 | Weighted residue-pair library (percent identity) with consistency extension. Zero-weight pairs are not stored; no documented weight-threshold pruning was found | L | – | L | – | – | – |
| Landan & Graur, HoT, *MBE* 24:1380 (2007), doi:10.1093/molbev/msm060 | Agreement between alignments of the sequences in original (heads) and reversed (tails) order | L | – | – | – | – | L |
| Penn et al., GUIDANCE, *MBE* 27:1759 (2010), doi:10.1093/molbev/msq066 | Residue-pair score is "the proportion of MSAs where this pair is aligned together" over bootstrap-guide-tree MSAs. Residue and column scores are averages. Server default removes columns below 0.93 (about 12% FPR and 78% TPR) | **C** | – | L (aggregates pairs to columns) | – | – | – |
| Sela et al., GUIDANCE2, *NAR* 43:W7 (2015), doi:10.1093/nar/gkv318 | Adds gap-penalty sampling and HoT co-optimals (400 alternatives by default). Pairs are not weighted | **C** | – | – | – | – | – |
| Chang, Di Tommaso & Notredame, TCS, *MBE* 31:1625 (2014), doi:10.1093/molbev/msu117 | Transitive consistency of each aligned pair with a library; filters or up-weights reliable columns | **C** | – | L | – | – | L |
| Do et al., ProbCons, *Genome Res* 15:330 (2005), doi:10.1101/gr.2821705; Roshan & Livesay, Probalign, *Bioinformatics* 22:2715 (2006), doi:10.1093/bioinformatics/btl472 | Pair-HMM or partition-function posterior that two residues are aligned; probabilistic consistency | – | L (a posterior for each pair, but from a sequence model, not from votes) | – | – | – | – |
| Bradley et al., FSA, *PLoS CB* 5:e1000392 (2009), doi:10.1371/journal.pcbi.1000392 | Posterior-based sequence annealing | – | L | – | – | – | – |
| Kim & Ma, PSAR, *NAR* 39:6359 (2011), doi:10.1093/nar/gkr334 | Reliability = agreement with sampled suboptimal alignments | L | – | – | – | – | – |
| Wu, Chatterji & Eisen, ZORRO, *PLoS ONE* 7:e30288 (2012), doi:10.1371/journal.pone.0030288 | Column confidence from pair-HMM posteriors; masking | L | L | – | – | – | – |
| Ali, Bogusz & Whelan, Divvier, *MBE* 36:2340 (2019), doi:10.1093/molbev/msz142 | Pair-HMM homology posteriors per column, then UPGMA clustering with a cutoff calibrated at 1% FDR on BAliBASE; splits or filters | C/L (threshold the evidence, then cluster) | L | – | – | – | – |
| Capella-Gutiérrez et al., trimAl, *Bioinformatics* 25:1972 (2009), doi:10.1093/bioinformatics/btp348 | Consistency score across several alignments. `automated1` "implements a heuristic to decide the most appropriate mode" for each alignment | L | – | – | – | – | C/L |
| Tan et al., *Syst Biol* 64:778 (2015), doi:10.1093/sysbio/syv033 | Automated filtering often worsens single-gene trees | – | – | – | – | – | L (motivates a gate) |

### D. Unsupervised aggregation outside MSA (method-level precedents)

| Paper | What it does | (1) | (2) | (3) | (4) | (5) | (6) |
|---|---|---|---|---|---|---|---|
| Dawid & Skene, *Applied Statistics* 28:20 (1979), doi:10.2307/2346806 | EM estimate of per-annotator error rates with no ground truth | – | C | – | **I (as a method)** | – | – |
| Chen, Mackey, Vermunt & Roos, *PLoS ONE* 2:e383 (2007), doi:10.1371/journal.pone.0000383 | Latent class analysis estimates the sensitivity and specificity of orthology (homology) methods with no gold standard, and assesses dependence between methods | – | C | – | **C** | L | – |
| Cantarel et al., BAYSIC, *BMC Bioinf* 15:104 (2014), doi:10.1186/1471-2105-15-104 | Unsupervised Bayesian latent class model over variant callers; user sets a threshold on the posterior | C | **C** | – | **C** | – | – |
| Parisi et al., SML, *PNAS* 111:1253 (2014), doi:10.1073/pnas.1219097111 | Ranks and combines predictors without labels (spectral method) | – | – | – | C | – | – |
| Jaffe et al., *AISTATS*, PMLR 51:351 (2016) | Unsupervised ensembles with conditionally dependent classifiers (DREAM somatic mutation calling) | – | – | – | C | **C** | – |
| Kallus et al., ROPE, arXiv:1702.07685 (2017) | Over-dispersed beta-binomial mixture (null and non-null) on edge-selection counts across resamples, giving edge q-values | L | **C** (same model family, "votes out of B" for graph edges) | – | – | – | – |
| Pledger, *Biometrics* 56:434 (2000), doi:10.1111/j.0006-341x.2000.00434.x | Finite binomial mixtures for closed capture–recapture, where all-zero histories are unobserved (zero truncation) | – | **C** (the zero-truncated binomial mixture machinery) | – | – | – | – |
| Henikoff & Henikoff, *JMB* 243:574 (1994), doi:10.1016/0022-2836(94)90032-9; Kish, *Survey Sampling* (1965) | Sequence redundancy weights; effective sample size | – | – | – | – | L | – |

Search found no MSA paper that applies Dawid–Skene, latent class, or a mixture model to agreement counts. Several queries were
tried: Dawid–Skene plus MSA, latent class plus residue pairs, mixture model plus alignment agreement, crowdsourced alignment
aggregation (Phylo), and word-alignment combination.

## Verdicts

1. **(1) Support filter: KNOWN as a technique; no application to GCM was found.** Lassmann & Sonnhammer (2007) use essentially
   the same rule (keep residue pairs found in at least f input alignments before merging). GUIDANCE/GUIDANCE2, TCS and Muscle5
   CC threshold support fractions. Search found the rule in neither the MAGUS papers nor the MAGUS code. Claim it only as "the
   first application to GCM's alignment graph".
2. **(2) Zero-truncated binomial / beta-binomial mixture posterior: PARTIALLY NOVEL.** This is the strongest claim. Search found
   no MSA or merger paper with an unsupervised mixture on vote counts. The machinery is standard elsewhere: Pledger 2000
   (zero-truncated binomial mixtures), ROPE (beta-binomial mixture on resampling counts of graph edges), BAYSIC (latent-class
   posterior plus threshold). Muscle5's note that CC over-estimates as a probability is useful motivation.
3. **(3) Strength-weighted votes: PARTIALLY NOVEL, minor.** Weighted votes inside GCM already exist (WITCH and WITCH-NG,
   per-alignment weights), as do method or pair weights in consensus (M-Coffee, T-Coffee). Search did not find this exact
   normalisation (each backbone's fraction of realised residue pairs). It is a re-normalisation of GCM's existing count weight.
4. **(4) Dawid–Skene per-backbone reliability: PARTIALLY NOVEL.** Search found no application to MSA. The method is classical
   (Dawid & Skene) and has been used for homology and variant calls (Chen 2007, BAYSIC). M-Coffee found that per-method weights,
   including accuracy-based ones, gave no significant gain. MUMSA's MOS is a one-step analogue. MAGUS backbones are
   exchangeable (same aligner, random subsets), so expect reliabilities near uniform.
5. **(5) Overlap-aware effective number of votes: PARTIALLY NOVEL.** The principle of discounting correlated voters is known:
   M-Coffee VarCov and tree weights, Jaffe 2016, WITCH's overlap prior, Henikoff weights, Kish's effective sample size. Search
   found no overlap-derived effective vote count for GCM backbones.
6. **(6) Reference-free dataset-level gate: PARTIALLY NOVEL, weak.** The concept is close to MUMSA AOS (difficulty from
   inter-alignment agreement), Muscle5 dispersion, Facet aligner advising, AQUA and trimAl `automated1`. Search found no gate
   that decides whether to filter GCM's evidence. Reviewers will ask to compare against an AOS or dispersion statistic computed
   on the 10 backbones over shared sequences.

## Must-cite (6)

1. Zaharias, Smirnov & Warnow, TCBB 2023, doi:10.1109/TCBB.2022.3191848. It defines the weights and states the future-work hook.
2. Lassmann & Sonnhammer 2005, doi:10.1093/nar/gki1020, together with 2007, doi:10.1186/1471-2105-8-S5-S9. They give the AOS
   gate, the MOS reliability score and the f-of-m support rule.
3. Wallace et al. 2006, M-Coffee, doi:10.1093/nar/gkl091. It tested method weighting and redundancy weighting, with a negative
   result.
4. Shen, Park & Warnow 2022, WITCH, doi:10.1089/cmb.2021.0585. It is the existing weighted GCM with statistically defined
   per-source weights.
5. Penn et al. 2010, GUIDANCE, doi:10.1093/molbev/msq066, or Sela et al. 2015, GUIDANCE2. It defines residue-pair support as a
   fraction of perturbed alignments, plus a threshold.
6. Dawid & Skene 1979, doi:10.2307/2346806. Pair it with a mixture precedent: Kallus et al. ROPE (arXiv:1702.07685) or Pledger
   2000.

Strongly recommended as well: MAGUS (the base method), Muscle5 (CC and dispersion), TCS, MergeAlign, BAYSIC, Jaffe et al. 2016,
and Kececioglu & DeBlasio 2013.

## Safe and unsafe wording

- **Safe:** "To our knowledge, backbone-agreement filtering and unsupervised vote-count models have not been applied to GCM's
  alignment graph." Also: "Following WITCH's weighted GCM and M-Coffee's method weighting, we estimate reliability from
  agreement alone."
- **Avoid:** "first agreement-based filter" (Lassmann 2007, GUIDANCE, TCS). Also avoid "first weighted GCM" (WITCH), "new
  consistency score", and "first reference-free difficulty predictor" (MUMSA AOS, Muscle5 dispersion, Facet).

## Not accessed or unverified

- Smirnov's PhD thesis and the MAGUS supplement.
- Full text of DeBlasio & Kececioglu ACM-BCB 2015.
- The 2026 bioRxiv "A deep learning–based score to evaluate multiple sequence alignments". It cites MWT-AM and may bear on (6);
  only the title was seen.
- Whether ROPE was published in a journal.
- T-Coffee source code, checked for any undocumented library-pruning flag.
- The trimAl `automated1` decision rule. The paper defers it to the online documentation.
