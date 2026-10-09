# Is MAGUS the state of the art? A literature review for a CS581 MSA project

*Prepared 2026-10-09. Starting point: Smirnov & Warnow, "MAGUS: Multiple sequence Alignment using Graph clUStering", Bioinformatics 37(12):1666–1672 (2021; online 2020-11-30), doi:10.1093/bioinformatics/btaa992.*

**How the citations were checked.** I resolved every journal DOI below against Crossref and every dataset DOI against DataCite on 2026-10-09. I read MAGUS's own source code (github.com/vlasmirnov/MAGUS) to confirm defaults. Numbers come from the papers' own text. Where a paper reports a number only in a figure, I say that rather than guess. Anything I could not verify is marked **(unverified)**.

---

## TL;DR

- **For the job it was built for, MAGUS is still near the top.** That job is aligning large datasets (about 1k–50k sequences) of full-length nucleotide sequences that evolve at high rates. On that data its accuracy is matched by UPP2 and TWILIGHT but not clearly beaten. Its merger, GCM, is still the best published way to merge many disjoint alignments at once.
- **No better merger has been published since.** Warnow's group has put out nothing new on merging after the 2023 TCBB paper on the MWT-AM problem. Her 2025 talk lists *"What are the best ways to merge disjoint alignments?"* as an open problem.
- **Outside that job, MAGUS is not the state of the art:**
  - *Fragmentary or length-heterogeneous data:* UPP2, WITCH/WITCH-NG and MAGUS+eHMMs do better. EMMA does better at adding sequences to an existing alignment.
  - *Very large or long nucleotide data:* TWILIGHT (2025) matches MAGUS's accuracy at 10k RNASim sequences and aligns 1M in about 30 minutes. MAGUS failed to finish 200k within 24 hours.
  - *Proteins:* learnMSA2 (which uses a protein language model), Muscle5, FAMSA2 (Nat. Biotechnol. 2026) and regressive T-Coffee. Language-model aligners (vcMSA, ARIES) lead at low sequence identity.
  - *Speed:* MAGUS is one of the slowest methods in recent benchmarks.
- **Best project niche: the merge step.** The three strongest ideas, all unpublished as far as I can find, are in Section 4:
  1. Exact pairwise merges by dynamic programming on GCM's edge weights.
  2. Better edge weights combined with backbones chosen to cover every column.
  3. A weight-aware trace search plus an exact ILP baseline that measures how far GCM is from optimal.

---

## 1. What MAGUS does

### 1.1 Pipeline

Paper defaults are shown alongside the current code defaults from `magus/main.py`.

1. **Guide tree.**
   - A "skeleton" of 300 roughly full-length sequences is selected. The paper says their lengths are within 25% of the 75th-percentile length.
   - The skeleton is aligned with MAFFT, and an HMM is built on that alignment with `hmmbuild`.
   - The remaining sequences are added with `hmmalign`, and FastTree2 builds the tree.
   - Alternatives (`-t`): FastTree without ML, Clustal Omega's mBed tree, MAFFT PartTree, or a user-supplied tree.
2. **Decomposition.**
   - The tree is split PASTA-style by repeatedly deleting centroid edges.
   - Splitting stops when either `--maxnumsubsets` (default 25) is reached or every subset is smaller than `--maxsubsetsize` (default 50).
   - The paper used 25 subsets, 50 for larger datasets, and 200 for 16S.B.ALL.
3. **Subset alignment.** Each subset is aligned with MAFFT L-INS-i, the same as PASTA.
   - In Recursive MAGUS, a subset larger than `--recursethreshold` (default 200) is aligned by MAGUS itself.
4. **Backbones (GCM step 0).**
   - `-r 10` backbones, each with at most `-m 200` sequences, are aligned with MAFFT.
   - Each backbone draws the same number of sequences from every subset: ⌊200/25⌋ = 8 per subset with the defaults.
   - The draw is **uniformly random** (`random.shuffle` in `graph_build/graph_builder.py`).
   - "Slow" mode extends the backbones to all sequences with HMMER (`--graphbuildhmmextend`).
   - Other backbone sources: `--graphbuildmethod mafftmerge|initial|subsethmm`, or user files with `-b`.
5. **Alignment graph (step 1).**
   - Each node is one column of a subset ("constraint") alignment.
   - The weight on edge (u,v) is the number of letter pairs, one letter from column u and one from column v, that share a column in some backbone, summed over all backbones.
   - The weights are raw counts and are not normalized.
6. **Clustering (step 2).** Markov Clustering (MCL) with inflation 4 (`-f`). The paper cites 1.4–6 as the recommended range.
7. **Trace: cleaning ("trimming") and ordering (step 3).** Default `--graphtracemethod minclusters`.
   - **(a) Cleaning.** Some clusters break the rules of an alignment column. The paper calls these vertical and horizontal violations. The code checks for two cases: two columns from the same subset in one cluster, and one column in more than one cluster. Violations are removed greedily, starting with the element that has the smallest total weight to the rest of its cluster. Removed elements become singletons (`clean_clusters.py`).
   - **(b) Ordering.** An A\* search orders the clusters. If no cluster can be placed next, because clusters cross each other, a cluster is split.
   - The code comment states the objective: *"A\* to search for the path of cluster breaks with the smallest number of clusters broken."* So the search minimizes the **number** of clusters broken, not the **weight** lost.
   - When the search stalls, it switches to weighted A\*, which gives up the optimality guarantee.
   - Alternatives come from the MWT-AM paper: `fm`, `mwtgreedy`, `mwtsearch`, `rg`, `rgfast`. An optional local-search optimizer (`--graphtraceoptimize`) can be run afterwards.
8. **Output.** By default (`-c true`) the output keeps every subset alignment intact, so the merge only adds homologies *between* subsets.
   - MAGUS runs a single iteration. Unlike PASTA, it does not re-estimate the tree and realign.
   - Building the backbones dominates the running time.

### 1.2 Headline results of the MAGUS paper

- **ROSE 1000-taxon data.** MAGUS beat PASTA on 9 of 10 model conditions and tied on 1.
  - The largest gap was on 1000L3: about 10–11% error for MAGUS versus 17–18% for PASTA (error = mean of SPFN and SPFP).
  - MAGUS was about twice as fast (about 40 → 15 minutes on 1000M2).
- **16S.** Tied on 16S.3 and was 1–2 points better on the other three datasets. On 16S.B.ALL it took 3 hours versus 6 for PASTA.
- **BAliBASE (8 protein sets).** A clear advantage on 3 sets and ties or near-ties on the rest.
- **PASTA+GCM** beat default PASTA on the hardest conditions. This shows the gain comes from the merger.

### 1.3 The MAGUS family (2021–2023)

- **Recursive MAGUS.** Smirnov, *PLOS Comput Biol* 17(10):e1008950, 2021, doi:10.1371/journal.pcbi.1008950.
  - Adds recursion, alternative guide trees, multi-node parallelism and compression of alignments over 100 GB.
  - HomFam (19 families of 10,099–93,681 sequences): average SP error 15.5% for MAGUS 1, 16.5% for recursive MAGUS with FastTree, 17.9% with the Clustal tree. MAFFT, UPP and PASTA were about 21–23%; Muscle (v3) was 46.6%.
  - RNASim 1M: only UPP(Fast) (about 77 h) and MAGUS (recursive, Clustal tree; about 128 h, about 8.3% error) finished within one week.
  - **Recursion buys scalability, not accuracy.**
- **MAGUS+eHMMs.** Shen, Zaharias & Warnow, *Bioinformatics* 38(4):918–924, 2022, doi:10.1093/bioinformatics/btab788.
  - MAGUS turns out to be "fairly robust" to fragmentary sequences.
  - Aligning the full-length "backbone" sequences with MAGUS and adding the rest with UPP's ensemble of HMMs beats UPP.
- **MWT-AM.** Zaharias, Smirnov & Warnow, "Large-Scale Multiple Sequence Alignment and the Maximum Weight Trace Alignment Merging Problem", *IEEE/ACM TCBB* 20(3):1700–1712, 2023, doi:10.1109/TCBB.2022.3191848. The conference version is AlCoB 2021.
  - Defines MWT-AM: find the merged alignment that maximizes the total edge weight placed in the same column. It generalizes Kececioglu's Maximum Weight Trace problem.
  - Proves MWT-AM is NP-hard.
  - Shows the MWT-AM score correlates with accuracy, but how strongly depends on difficulty:
    - With random decompositions, Spearman ρ was 0.851–0.997 on 1000M1 (hard) and 0.241–0.971 on 1000M4 (easier).
    - With MAGUS's default PASTA-style decomposition, ρ = 0.9 on 1000M1 but about 0 (−0.048) on 1000M4.
    - So the objective matters most on hard data.
  - Tests other clustering steps: MLR-MCL and Region Growing (RG).
  - Tests other trace steps: Fiduccia–Mattheyses (FM), MWTgreedy, MWTsearch and RG-fast, plus a greedy "optimizer" that moves nodes between columns.
  - **GCM(fm+opt) beat default GCM on both MWT-AM score and SP accuracy on every HomFam dataset.**
  - MAFFT-merge with L-INS-i was "generally as accurate as GCM" but more expensive, and it failed on the two largest HomFam sets (88,345 and 93,681 sequences).
  - T-Coffee did poorly at merging nucleotide alignments.
  - The paper also notes that MAGUS (non-recursive) was more accurate than recursive MAGUS.

### 1.4 Limitations

**Acknowledged by the authors (MAGUS, Recursive MAGUS, MWT-AM papers):**

- **Missing comparisons.** The MAGUS paper left Clustal Omega, regressive T-Coffee and Kalign for future work.
- **Assumes little length heterogeneity.** Biological test sets were filtered to within 20% of the median length. Fragmentary data was left to MAGUS+eHMMs.
- **"Slow" HMM-extended mode** adds running time with high variance.
- **The A\* ordering can stall** on "tangled" clusters. The README warns that more than 100 subalignments slow the ordering phase a lot, especially on heterogeneous data.
- **GCM optimizes no explicit objective.** MWT-AM is NP-hard, and there are no known approximation results.
- **Recursion does not improve accuracy.** Without recursion, MAGUS needs MAFFT-linsi to handle subsets of about 1,000 nucleotide or more than 4,000 amino-acid sequences.
- **Stated future work:**
  - Change how GCM defines its weights.
  - Use other base aligners: PROMALS, ProbCons, BAli-Phy.
  - Handle fragmentary and genome-scale data.
  - Find approximation algorithms for MWT-AM.

**Observed by others since:**

- **Speed.**
  - Slowest CPU method in learnMSA2's HomFam benchmark (0.19 h per family).
  - Timed out on the largest Pfam families (PF00005, PF07690) in the learnMSA (2022) tests.
  - Failed on 200k RNASim sequences within 24 hours in TWILIGHT's tests.
  - TWILIGHT's talk slides give 9.9 hours for MAGUS versus 44 minutes for MAFFT on 100k sequences.
- **Over-alignment (column inflation).**
  - learnMSA's authors note the "tendency of the divide-and-conquer methods (T-Coffee, MAGUS) to construct MSAs with much larger column counts".
  - On 100k RNASim sequences, TWILIGHT's alignment file was 1.6 GB versus 21.9 GB for MAGUS.

---

## 2. What came after MAGUS

### 2.1 Summary table

| Method | Citation (verified DOI) | Data | Claim relative to MAGUS / SOTA |
|---|---|---|---|
| WITCH | Shen, Park, Warnow, *J Comput Biol* 29(8):782–801, 2022, doi:10.1089/cmb.2021.0585 | NT/AA, fragmentary | Better than UPP. Merges the top-10 HMM extended alignments with a **weighted GCM** |
| WITCH-NG | Liu & Warnow, *Bioinform Adv* 3:vbad024, 2023, doi:10.1093/bioadv/vbad024 | same | Same accuracy as WITCH, ≥5× faster on 7/10 HomFam sets. Replaces MCL+A\* with an **exact Smith–Waterman MWT on two alignments** |
| UPP2 | Park, Ivanovic, Chu, Shen, Warnow, *Bioinformatics* 39(1):btad007, 2023, doi:10.1093/bioinformatics/btad007 | NT/AA | Best with fragments at high rates. **Ties MAGUS on large 16S; MAGUS slightly better on HomFam** |
| HMMerge | Park & Warnow, *Bioinform Adv* 3:vbad052, 2023, doi:10.1093/bioadv/vbad052 | very short seqs | Competitive with WITCH; better for adding very short sequences |
| EMMA | Shen, Liu, Williams, Warnow, *Algorithms Mol Biol* 18, 2023, doi:10.1186/s13015-023-00247-x (WABI 2023; bioRxiv 10.1101/2023.06.12.544642) | adding seqs to a backbone | Scales MAFFT-linsi--add by divide-and-conquer with **transitivity merging**. Beats WITCH-NG-add; up to about 186k sequences. Not tested de novo, and not compared with MAGUS |
| BAli-Phy as phylogeny-aware aligner | Gupta, Zaharias, Warnow, *Bioinformatics* 37:4677–4683, 2021, doi:10.1093/bioinformatics/btab555 | NT, ≤1000 seqs | More accurate than MAFFT and PRANK with a fixed tree |
| Muscle5 | Edgar, *Nat Commun* 13:6968, 2022, doi:10.1038/s41467-022-34630-w | AA, RNA | Highest accuracy on BAliBASE/Bralibase (developer's claim). Ensembles with column confidence; `-disperse` measures ensemble dispersion; Super5 scales to tens of thousands. No comparison with MAGUS |
| learnMSA | Becker & Stanke, *GigaScience* 11:giac104, 2022, doi:10.1093/gigascience/giac104 | AA, ≥1000 seqs | Matches SOTA on HomFam. MAGUS times out on the 2 largest Pfam families and wins on PF00096 |
| learnMSA2 | Becker & Stanke, *Bioinformatics* 40(S2):ii79–ii86, 2024, doi:10.1093/bioinformatics/btae381 | AA | ProtT5 embeddings: HomFam TC 66.00 / SP 87.30, "almost 6% points" more columns than the best amino-acid-only competitor (MAGUS was in the field). MAGUS's own numbers are only in Fig. 3a |
| vcMSA | McWhite, Armour-Garb, Singh, *Genome Res* 33(7):1145–1153, 2023, doi:10.1101/gr.277675.123 | AA, small sets | Clusters and orders language-model residue embeddings. Best on more QuanTest2 sets (147 sets × 20 seqs). Limited scalability |
| ARIES | Hoang, Armour-Garb, Singh, *Genome Res*, online 2026-10-01, doi:10.1101/gr.282254.126 (bioRxiv 10.64898/2026.01.02.697423) | AA | Builds a template embedding and aligns each sequence to it with dynamic time warping. Near-linear scaling; beats SOTA especially at low identity (abstract-level claim) |
| BetaAlign | Dotan et al., *Bioinformatics* 2025, doi:10.1093/bioinformatics/btaf009 (ICLR 2023 version) | small sets | Transformer treating alignment as translation, trained on simulations. Not a large-scale method |
| FAMSA2 | Gudyś, Zieleziński, Notredame, Deorowicz, *Nat Biotechnol* 2026-04-14, doi:10.1038/s41587-026-03095-3 (bioRxiv 10.1101/2025.07.15.664876) | AA, up to millions | extHomFam (390 families, up to 3M seqs): SP 79.6 / TC 61.8 in 15 h. Muscle5, T-Coffee and Clustal Omega needed 8–37 days and failed on the largest family. **MAGUS not tested** |
| Regressive T-Coffee | Garriga et al., *Nat Biotechnol* 37:1466–1470, 2019, doi:10.1038/s41587-019-0333-6 | AA | Up to 1.5M sequences; accurate at ≥10k. Predates MAGUS; could not be run in the Recursive MAGUS study |
| TWILIGHT | Tseng, Walia, Turakhia, *Bioinformatics* 41(S1):i332–i341, 2025 (ISMB), doi:10.1093/bioinformatics/btaf212 | NT, huge/long | See 2.3. **The only post-2023 paper to benchmark MAGUS on nucleotide reference data** |
| BAli-Phy 3 | Redelings, *Bioinformatics* 37:3032–3034, 2021, doi:10.1093/bioinformatics/btab129 | small | Statistical co-estimation of alignment and tree; small datasets only |

### 2.2 Warnow group after MAGUS: where the merging work stands

- **Change in direction.** After the MWT-AM paper, the group moved to **sequence-length heterogeneity**: WITCH, WITCH-NG, UPP2, HMMerge and EMMA.
- **Two results bear directly on merging:**
  - **WITCH uses GCM as a weighted consensus merger.** Edge weights come from HMM support rather than raw backbone counts.
  - **WITCH-NG shows that MWT restricted to two alignments is solvable exactly.** It uses a Smith–Waterman-style dynamic program with zero gap penalty and −∞ for unsupported pairs. The paper describes this as the same objective GCM approximates heuristically.
- **2024–2026: nothing new on merging.**
  - I checked the lab's papers page and a Europe PMC list of 30 Warnow-authored items from 2023–2026. The newest MSA paper is EMMA (December 2023).
  - The 2024–2026 output covers phylogenetic placement (SCAMPP, BSCAMPP), metagenomic profiling (TIPP3, TIPP-SD), phylogenetic networks (CAMUS) and seeding for viral genome alignment ("Min-frame transformation", bioRxiv 2026). None is about merging alignments.
  - **I found no method called "SPAM"**, and no GCM variants beyond those in the TCBB paper.
- **Warnow's 2025 talk** ("Advances in Large-scale MSA", tandy.cs.illinois.edu/warnow-MSA-2025.pdf) recommends:
  - Large numbers of sequences without length heterogeneity: *"MAGUS, TWILIGHT, FAMSA, and others."*
  - With length heterogeneity: *"UPP, WITCH(-ng), EMMA."*
  - Open challenges include *"What are the best ways to merge disjoint alignments?"* and ensemble use of multiple MSAs.

### 2.3 TWILIGHT (2025): the main newer competitor on nucleotides

- **Algorithm.**
  - Progressive alignment, run on CPU or GPU.
  - Columns that are more than 95% gaps are removed before profile alignment and restored afterwards.
  - Optional partitioning of the guide tree into subtrees, which are merged back by profile alignment or by **transitivity merging**.
  - An iterative mode re-estimates the tree with PartTree first, then FastTree/RAxML/IQ-TREE.
- **Head-to-head with recursive MAGUS:**
  - RNASim 10k: both about 8% error.
  - RNASim 1M: TWILIGHT reached 6.48% error in 32 minutes on CPU or 18 minutes on one GPU, under 16 GB of memory. MAGUS failed on 200k within 24 hours.
  - AliSim long-branched 10k: MAGUS had the lowest error (8.3%) versus about 9% for TWILIGHT, which was 4.6× faster.
- **What this means.** At moderate size and high rates, MAGUS is still at least as accurate. At very large scale, TWILIGHT wins easily.

### 2.4 Proteins

- **The most accurate protein aligners are protein-specific:**
  - Muscle5 for up to a few thousand sequences.
  - learnMSA2 for 1,000 or more sequences, with a GPU and a protein language model.
  - FAMSA2 for very large families; regressive T-Coffee for big families as well.
- **MAGUS is a strong amino-acid-only method on HomFam:**
  - Recursive MAGUS: 15.5% error versus about 21–23% for MAFFT, UPP and PASTA.
  - UPP2 paper: slight edge over UPP2.
- **But it has been beaten or skipped since:**
  - learnMSA2's language-model mode reportedly beats every amino-acid-only competitor, MAGUS included, by about 6 percentage points of TC.
  - FAMSA2 and Muscle5 never compared against MAGUS.
- **Caution:** most of these numbers come from each method's own developers.

### 2.5 MAFFT, reviews and independent benchmarks

- **MAFFT.** The latest version on the official page is 7.526 (April 2024). I found no new MAFFT journal paper from 2023–2026 **(not exhaustively verified)**.
- **Review of large-scale protein MSA.** Santus, Garriga, Deorowicz, Gudyś, Notredame, "Towards the accurate alignment of over a million protein sequences: current state of the art", *Curr Opin Struct Biol* 80:102577, 2023, doi:10.1016/j.sbi.2023.102577. They conclude a unified framework able to *consistently and efficiently* produce high-accuracy large alignments "is still lacking".
- **nf-core/multiplesequencealign.** Santus et al., *NAR Genom Bioinform* 7(3):lqaf104, 2025, doi:10.1093/nargab/lqaf104. A Nextflow benchmarking framework that bundles MAGUS, UPP2, learnMSA, FAMSA, Muscle5, Kalign3, regressive T-Coffee and others. Useful infrastructure, but it publishes no ranking.
- **MuSAlS.** Light et al., arXiv:2601.15458, January 2026. Reports MAGUS "failed to detect and align subalignments" on protein data. The evaluation does not use reference alignments, and the failure looks like a configuration problem. **Weak evidence.**
- **No independent comparison exists.** I found no 2024–2026 study that runs MAGUS, TWILIGHT, FAMSA2, Muscle5 and learnMSA2 on the same reference datasets. Filling that gap could be a modest add-on to your project.

---

## 3. Verdict by regime

| Regime | Strongest methods (evidence) | MAGUS's standing |
|---|---|---|
| Nucleotide, full-length, 1k–~50k seqs, high rate | MAGUS, UPP2 (tie on 16S), TWILIGHT (tie on RNASim 10k; MAGUS better on AliSim long-branch) | **Still SOTA-class for accuracy**, but much slower than TWILIGHT |
| Nucleotide, 100k–1M+, or very long sequences | TWILIGHT; UPP(Fast) and recursive MAGUS finish 1M only in days | Scales in principle but is impractical; **not SOTA** |
| Fragmentary or length-heterogeneous | UPP2, WITCH/WITCH-NG, MAGUS+eHMMs; EMMA for adding to a backbone; HMMerge for very short reads | MAGUS alone is not SOTA, but it is the default backbone aligner in these pipelines |
| Small (≤1000 seqs) | MAFFT L-INS-i/G-INS-i, BAli-Phy, T-Coffee/PROMALS (Warnow 2025 slides) | At this size MAGUS is close to MAFFT-linsi |
| Protein, 1k–100k | learnMSA2 + language model (GPU), Muscle5, MAGUS, UPP2 | Competitive among amino-acid-only methods; beaten by learnMSA2 per its authors |
| Protein, >100k to millions | FAMSA2, regressive T-Coffee, learnMSA2 | Times out on the largest Pfam families; **not SOTA** |
| Low-identity proteins | vcMSA, ARIES, learnMSA2 (language-model based) | Not competitive (MAFFT backbones are weak here) |
| Speed | FAMSA2, TWILIGHT, MAFFT-auto/PartTree, Kalign 3 | Among the slowest |
| **Merging k disjoint alignments** | GCM / GCM(fm+opt); MAFFT-merge (L-INS-i) similar accuracy but slower and fails at large k·n | **GCM is still the best published merger** |

**Bottom line for the student.** MAGUS is no longer the state of the art for MSA in general in 2026. It is still a top method for the hard case it targets: large, high-rate, full-length nucleotide data. Its merger is the main published tool for an open problem that Warnow herself highlights, which makes improving the merge step a well-motivated CS581 project.

---

## 4. Project ideas: improving the divide-and-conquer merge step

**Shared setup for every idea.**

- **Isolate the merger.** MAGUS accepts precomputed subset alignments (`-s`) and backbones (`-b`). Fix both and vary only the merge.
- **Swap components in code.** Clustering is in `merge/graph_cluster/`, the trace in `merge/graph_trace/`, the optimizer in `merge/optimizer.py`, and edge weights in `graph_build/graph_builder.py`.
- **Data.**
  - ROSE 1000M1–M4 and 1000L1–L4, 5 of the 20 replicates each. MAGUS runs in minutes per replicate.
  - RNASim 10k (3 replicates) and 16S.T.
  - For proteins, the 10 largest HomFam sets or the 8 BAliBASE sets from the PASTA/MAGUS papers.
- **Metrics.**
  - SPFN and SPFP via FastSP (doi:10.1093/bioinformatics/btr553).
  - Expansion ratio (estimated alignment length / true alignment length).
  - MWT-AM score, TC score for proteins, running time and memory.
  - Optionally, tree error from FastTree.
- **Statistics.** Paired Wilcoxon tests across replicates.

### P1. Exact pairwise merges on GCM's weights, applied progressively (top pick)

- **Idea.** Keep MAGUS's subsets, backbones and edge weights, but replace MCL plus A\* with a sequence of *exact* two-alignment MWT merges.
  - Merging alignments A and B is a dynamic program over their columns, scored by w(c,c′) with no gap penalty. This is the trick WITCH-NG used.
  - Merge in a chosen order: up the guide tree, or greedily by the largest total inter-alignment weight, UPGMA-style.
  - After each merge, a new column's edges are the sums of its member columns' edges.
  - Optionally finish with the MWT-AM optimizer.
- **Why it might help.** Every step is optimal and runs in polynomial time (O(L_A·L_B)). It avoids MCL's blindness to the alignment constraints and the A\* "tangled cluster" failure mode.
  - Risk: errors made early in a progressive merge are locked in. The scientific question is whether exact-but-greedy beats heuristic-but-global.
- **Experiment.** Compare default GCM, GCM(fm+opt) and pairwise-exact merging with 2–3 merge orders. Use identical subsets and backbones, and vary the number of subsets (10, 25, 50).
- **Already done?** Partially, in neighboring settings:
  - WITCH-NG solves the two-alignment MWT exactly, but for adding a query to a backbone, not for de novo merging of subsets.
  - PASTA merges pairs using sequence-based scores (Opal or Muscle), then applies transitivity. Opal is Wheeler & Kececioglu, doi:10.1093/bioinformatics/btm226.
  - TWILIGHT merges subtrees by profile alignment.
  - **I found no published method that combines backbone-derived MWT weights with exact pairwise merging in a de novo pipeline.**

### P2. Better edge weights: normalization, confidence and consistency

- **Idea.** Today an edge weight is a raw count of co-aligned letter pairs. Alternatives:
  - **(a) Normalize** each count by how often both columns had letters present in the same backbone. This stops gappy columns from being penalized just because they have fewer sampled letters.
  - **(b) Confidence-weight each backbone**, using agreement across backbones or Muscle5-style column confidence from backbone replicates.
  - **(c) Consistency transformation** before clustering: T-Coffee/ProbCons-style triplet extension through third subsets.
  - **(d) Posterior probabilities** instead of hard 0/1 alignments, for example HMMER posterior decoding in slow mode.
- **Why it might help.** MWT-AM tracks accuracy only as well as the weights reflect true homology. The MWT-AM paper names changing the weights as its first future direction.
- **Experiment.** Ablate each weighting with fixed subsets and backbones. Report SPFN and SPFP separately, and the correlation between MWT-AM score and accuracy under each scheme.
- **Already done?** Partially:
  - WITCH's weighted GCM uses HMM-derived weights, but for adding queries.
  - Muscle5 computes column confidence but never uses it as merge weights.
  - Consistency transformations are classic but have not been applied to MAGUS's column graph, as far as I found.

### P3. Coverage-aware backbone selection

- **Idea.** By default each backbone samples 8 random sequences per subset. A subset-alignment column that is gappy in all sampled sequences gets **no edges** and ends up a singleton. That plausibly drives the column inflation reported for MAGUS. Alternatives:
  - Choose backbone sequences with a greedy set cover, so every column of every subset is covered by at least t letters across all backbones.
  - Choose diverse sequences: farthest-point or medoid sampling within each subset.
- **Why it might help.** It targets missed homologies (SPFN) and over-alignment directly. It might also reach the same accuracy with fewer backbones, and backbones dominate the running time.
- **Experiment.**
  - Count singletons and measure the expansion ratio before and after the change.
  - Measure SPFN and SPFP while varying the number of backbones (5, 10, 15), comparing random and coverage-aware selection.
- **Already done?** The MAGUS paper tuned the *number and size* of backbones, always with random selection. Slow mode raises coverage with HMM extension, at high cost. **I found no published work on the selection strategy itself.** Low risk; pairs well with P2.

### P4. Weight-aware conflict resolution with an exact baseline

- **Idea.**
  - Replace minclusters' objective (fewest cluster breaks) with weight lost, using A\* or beam search.
  - Add an exact ILP for MWT-AM on small instances (k ≤ about 10 subsets, a few hundred columns each), using the trace/mixed-cycle formulation (Reinert et al., RECOMB 1997, doi:10.1145/267521.267845; Althaus et al., *Math Program* 2006, doi:10.1007/s10107-005-0659-3), or a simpler ordering-variable model solved with OR-Tools or Gurobi.
  - Add large-neighborhood search: free a window of columns and re-solve it exactly.
- **Why it might help.** Measuring GCM's optimality gap tells you whether the bottleneck is the *search* (P4) or the *weights* (P2). Nobody has measured this gap.
- **Experiment.** On ROSE data with 10–25 subsets, compare minclusters, mwtsearch, fm+opt, weight-aware search, the ILP (where it solves) and the neighborhood search on MWT-AM score, SP error and time.
- **Already done?** Partially. The MWT-AM paper tested FM, MWTgreedy, MWTsearch, RG, RG-fast and the optimizer. Exact branch-and-cut for classical MWT exists (Kececioglu, CPM 1993, doi:10.1007/BFb0029800). **I found no application of an exact solver to MAGUS's MWT-AM instances.** You must go beyond the TCBB paper.

### P5. Constraint-aware graph clustering in place of MCL

- **Idea.** Replace MCL with one of:
  - Leiden/CPM clustering.
  - Correlation clustering.
  - Agglomerative clustering with "cannot-link" constraints between columns of the same subset and order constraints.
  - Or choose MCL's inflation per dataset by maximizing the post-trace MWT-AM score.
- **Why it might help.** MCL ignores the rule of at most one column per subset and ignores column order, so the trace step has to repair the clusters afterwards.
- **Experiment.** Count the repairs needed, then compare MWT-AM score, SP error and time.
- **Already done?** Partially. The MWT-AM paper tested MLR-MCL, Region Growing (a constraint-aware greedy merge) and several inflation values. **I found no tests of Leiden or correlation clustering.**
  - Caution: in their results table, RG without the optimizer scored far lower MWT-AM (about 68M versus about 99M on one dataset). Constraint-awareness alone is not enough.

### P6. Adaptive or overlapping decomposition

- **Idea.**
  - Choose the number of subsets from an estimated rate of evolution, such as the median p-distance in the skeleton alignment. Easy data would get fewer, larger subsets, which MAFFT-linsi aligns well.
  - Or decompose by sequence similarity (k-mer or medoid clustering, as in FAMSA2's medoid trees).
  - Or use **overlapping** subsets that share anchor sequences, which directly supplies cross-subset homology edges.
- **Experiment.** Sweep k ∈ {5, 10, 25, 50} on a hard condition (1000M1) and an easier one (1000M4); the MWT-AM paper describes M1 as harder than M4. Show the optimum shifts with rate, then fit and test a simple rule for choosing k.
- **Already done?** Partially:
  - MAGUS used 25, 50 or 200 subsets.
  - Recursive MAGUS compared guide trees (FastTree, FastTree without ML, Clustal, PartTree, random).
  - The MWT-AM paper used random subsets of 50, 100 and 200.
  - **I found nothing published on adaptive k or overlapping subsets inside a GCM pipeline.** PASTA's pairwise-plus-transitivity merge is the closest relative.

### P7. Post-merge column compaction (fixing over-alignment)

- **Idea.** After the trace, merge adjacent columns that are *compatible* (no sequence has letters in both) when backbone evidence or residue identity supports it. Or realign small windows that contain many gappy columns.
- **Why it might help.** Inflated alignments have been reported independently by learnMSA's and TWILIGHT's authors. Homologies split across columns raise SPFN and can hurt tree estimation.
- **Experiment.** Measure the expansion ratio, SPFN, SPFP and ML tree error on ROSE data.
- **Already done?** The MWT-AM optimizer moves single nodes between columns; I found **no dedicated compaction step** for MAGUS. Small and safe: good as a second contribution.

### P8. Protein language-model similarity as extra graph edges (protein-only, higher risk)

- **Idea.** Add edges between subset columns weighted by the similarity of their mean ESM-2 or ProtT5 residue embeddings. Combine these with the backbone counts, or use embeddings to choose backbone sequences.
- **Why it might help.** Embeddings help most at low identity (learnMSA2, vcMSA, ARIES), which is where MAFFT backbones are weakest.
- **Experiment.** HomFam subsets (≤10k sequences) and BAliBASE. Compare MAGUS, MAGUS with embedding edges, and learnMSA2. Needs a GPU.
- **Already done?** Embeddings have been used for whole-alignment construction (vcMSA, learnMSA2, ARIES) and for pairwise alignment (PEbA/EBA), but **I found no use inside a divide-and-conquer merger.**

### P9. Theory (optional and risky)

- **Idea.** Approximation or inapproximability results for MWT-AM. The TCBB paper states none are known. Alternatively, an LP-relaxation-and-rounding heuristic evaluated empirically.
- **Feasibility.** Course-scale only if the scope is restricted, for example a fixed number of subsets or special graph structure.

**Recommendation.**

- Strongest course project: **P1**, or **P2 + P3** together.
- **P4** makes a good analytical companion to either.
- **P6 and P7** are safer fallbacks.
- Each fits in a semester on the 1000-taxon ROSE data plus one larger dataset, using the hooks MAGUS already has for supplying subsets and backbones.

---

## 5. Where the benchmark datasets live

**Caveat on access.** Illinois Data Bank pages returned 403 to automated fetching. I confirmed the dataset DOIs and descriptions through DataCite, but not the file sizes. Sizes marked "~" are estimates from sequence counts.

| Dataset | Content and size | Where |
|---|---|---|
| **ROSE 1000-taxon** (SATé, Liu et al., *Science* 2009, doi:10.1126/science.1171243) | 1000 sequences × 20 replicates per condition. 1000L1–L4, M1–M4 and S1–S4 were used in the MAGUS/MWT-AM papers. About 1 kb per sequence (**unverified**), so a few MB per replicate | MAGUS data: Illinois Data Bank IDB-2643961, doi:10.13012/B2IDB-2643961_V1 (`Datasets.zip`, which also has RNASim 1k/10k, 16S and BAliBASE). Original: https://sites.google.com/eng.ucsd.edu/datasets/alignment/sate-i (Google Drive) |
| ROSE fragmentary (HF/LF) | 7 ROSE conditions made fragmentary: 50% of sequences cut to about 25% length (HF), or 25% cut to about 50% (LF) | IDB-6128941, doi:10.13012/B2IDB-6128941_V1 |
| **RNASim** (Guo, Wang, Kim, arXiv:0912.2326; used from PASTA, doi:10.1089/cmb.2014.0156, onward) | 1,000,000 simulated RNA sequences, about 1,500 nt each (~1.5 GB unaligned), with a true tree and alignment. Common subsets: 1k (×20), 10k (×10) | Full set: `rnasim.tbz` from https://sites.google.com/eng.ucsd.edu/datasets/alignment/pastaupp. 10k ×10: IDB-4194451, doi:10.13012/B2IDB-4194451_V1. "RNASim2": Dryad doi:10.5061/dryad.95x69p8h8 (1.78 GB) |
| **CRW 16S** (Cannone et al., *BMC Bioinformatics* 2002, doi:10.1186/1471-2105-3-2) | 16S.3 (6,323 seqs), 16S.T (7,350), 16S.B.ALL (27,643; 24,246 after MAGUS's length filter), 16S.M (740 after filter) | Length-filtered versions in IDB-2643961. Cleaned versions plus high-fragmentation versions: IDB-2419626, doi:10.13012/B2IDB-2419626_V1. Unfiltered: IDB-1048258 and https://sites.google.com/eng.ucsd.edu/datasets/alignment/16s23s. 16S.B.ALL-100-HF: IDB-6604429 |
| **HomFam** (Sievers et al., *Mol Syst Biol* 2011, doi:10.1038/msb.2011.75) | Homstrad structural reference cores plus Pfam homologs. About 94 families (93–93,681 seqs); the 19–20 large ones have 10k–94k seqs | IDB-1048258, doi:10.13012/B2IDB-1048258_V1 (Recursive MAGUS: HomFam + 16S + RNASim). The 10 largest are also in the EMMA data, IDB-2567453 |
| extHomFam / extHomFam 2 | 2016 version about 415 MB (Harvard Dataverse doi:10.7910/DVN/BO2SVW). v2: Homstrad + Pfam 33.1, families from 200 to 3M sequences, about 3.9 GB | Zenodo doi:10.5281/zenodo.6524237 (a newer version exists) |
| **BAliBASE 3** (Thompson et al., *Proteins* 2005, doi:10.1002/prot.20527) | Reference protein alignments. MAGUS used 8 sets (195–732 seqs after filtering) | In IDB-2643961; "PASTA for Proteins" data IDB-4074787; http://www.lbgi.fr/balibase/ (**live status unverified**) |
| 10AA | 10 small biological protein datasets | EMMA data IDB-2567453; `small_10_aa.zip` on the PASTA Drive |
| EMMA extras | 5000M2/3/4-het (simulated with INDELible, length heterogeneity), serine recombinase Rec/Res (Res about 186k seqs) | IDB-2567453, doi:10.13012/B2IDB-2567453_V1 |
| QuanTest2 (Sievers & Higgins, *Bioinformatics* 2020, doi:10.1093/bioinformatics/btz552) | Accuracy scored by secondary-structure prediction (protein) | http://bioinf.ucd.ie/quantest2.tar (**live status unverified**) |

---

## 6. Key references (DOIs verified via Crossref or DataCite unless marked)

- **MAGUS lineage**
  1. Smirnov V, Warnow T. MAGUS. *Bioinformatics* 37(12):1666–1672, 2021. doi:10.1093/bioinformatics/btaa992
  2. Smirnov V. Recursive MAGUS. *PLOS Comput Biol* 17(10):e1008950, 2021. doi:10.1371/journal.pcbi.1008950
  3. Zaharias P, Smirnov V, Warnow T. Large-scale MSA and the MWT Alignment Merging problem. *IEEE/ACM TCBB* 20(3):1700–1712, 2023. doi:10.1109/TCBB.2022.3191848
  4. Shen C, Zaharias P, Warnow T. MAGUS+eHMMs. *Bioinformatics* 38(4):918–924, 2022. doi:10.1093/bioinformatics/btab788
- **Warnow-group follow-ups**
  5. Shen C, Park M, Warnow T. WITCH. *J Comput Biol* 29(8):782–801, 2022. doi:10.1089/cmb.2021.0585
  6. Liu B, Warnow T. WITCH-NG. *Bioinform Adv* 3:vbad024, 2023. doi:10.1093/bioadv/vbad024
  7. Park M et al. UPP2. *Bioinformatics* 39(1):btad007, 2023. doi:10.1093/bioinformatics/btad007
  8. Park M, Warnow T. HMMerge. *Bioinform Adv* 3:vbad052, 2023. doi:10.1093/bioadv/vbad052
  9. Shen C, Liu B, Williams KP, Warnow T. EMMA. *Algorithms Mol Biol* 18, 2023. doi:10.1186/s13015-023-00247-x
  10. Gupta M, Zaharias P, Warnow T. Phylogeny-aware alignment using BAli-Phy. *Bioinformatics* 37:4677–4683, 2021. doi:10.1093/bioinformatics/btab555
- **Other recent aligners**
  11. Edgar RC. Muscle5. *Nat Commun* 13:6968, 2022. doi:10.1038/s41467-022-34630-w
  12. Becker F, Stanke M. learnMSA. *GigaScience* 11:giac104, 2022. doi:10.1093/gigascience/giac104
  13. Becker F, Stanke M. learnMSA2. *Bioinformatics* 40(S2):ii79–ii86, 2024. doi:10.1093/bioinformatics/btae381
  14. McWhite CD, Armour-Garb I, Singh M. vcMSA. *Genome Res* 33(7):1145–1153, 2023. doi:10.1101/gr.277675.123
  15. Hoang M, Armour-Garb I, Singh M. ARIES. *Genome Res* 2026 (online 2026-10-01). doi:10.1101/gr.282254.126
  16. Gudyś A, Zieleziński A, Notredame C, Deorowicz S. FAMSA2. *Nat Biotechnol* 2026. doi:10.1038/s41587-026-03095-3
  17. Garriga E et al. Regressive T-Coffee. *Nat Biotechnol* 37:1466–1470, 2019. doi:10.1038/s41587-019-0333-6
  18. Tseng Y-H, Walia S, Turakhia Y. TWILIGHT. *Bioinformatics* 41(S1):i332–i341, 2025. doi:10.1093/bioinformatics/btaf212
  19. Dotan E et al. BetaAlign. *Bioinformatics* 2025. doi:10.1093/bioinformatics/btaf009
  20. Redelings BD. BAli-Phy v3. *Bioinformatics* 37:3032–3034, 2021. doi:10.1093/bioinformatics/btab129
- **Reviews and benchmarking frameworks**
  21. Santus L et al. Current state of the art for over a million protein sequences. *Curr Opin Struct Biol* 80:102577, 2023. doi:10.1016/j.sbi.2023.102577
  22. Santus L et al. nf-core MSA benchmarking framework. *NAR Genom Bioinform* 7(3):lqaf104, 2025. doi:10.1093/nargab/lqaf104
- **Classic merging and trace work**
  23. Kececioglu J. The maximum weight trace problem in MSA. CPM 1993, LNCS. doi:10.1007/BFb0029800
  24. Reinert K et al. Branch-and-cut for MSA. RECOMB 1997. doi:10.1145/267521.267845
  25. Althaus E, Caprara A, Lenhof H-P, Reinert K. Branch-and-cut for MSA. *Math Program* 105:387–425, 2006. doi:10.1007/s10107-005-0659-3
  26. Wheeler TJ, Kececioglu JD. Multiple alignment by aligning alignments (Opal). *Bioinformatics* 23(13):i559–i568, 2007. doi:10.1093/bioinformatics/btm226
  27. Modzelewski M, Dojer N. MSARC (residue clustering + FM). *Algorithms Mol Biol* 9:12, 2014. doi:10.1186/1748-7188-9-12
- **Foundational methods and evaluation**
  28. Mirarab S et al. PASTA. *J Comput Biol* 22(5), 2015. doi:10.1089/cmb.2014.0156
  29. Nguyen N, Mirarab S, Kumar K, Warnow T. UPP. *Genome Biol* 16:124, 2015. doi:10.1186/s13059-015-0688-z
  30. Mirarab S, Warnow T. FastSP. *Bioinformatics* 27(23):3250–3258, 2011. doi:10.1093/bioinformatics/btr553
- **Not peer-reviewed**
  31. Warnow T. "Advances in Large-scale Multiple Sequence Alignment" (talk slides, 2025, undated PDF). https://tandy.cs.illinois.edu/warnow-MSA-2025.pdf
  32. Light EG, Prior ME, Daniels NM, Ishaq N. MuSAlS. arXiv:2601.15458, 2026 (preprint; weak evidence about MAGUS)
