## Literature check: merging two multiple sequence alignments

I checked every DOI below against Crossref in Oct 2026. Where I describe how a method works, the source is the full text (PMC/Europe PMC or the author PDF) unless I say it comes from documentation or slides. Anything I could not confirm is marked **[unverified]**.

### 1. PASTA and SATé-II
- **Mirarab S, Nguyen N, Guo S, Wang L-S, Kim J, Warnow T.** PASTA: Ultra-large multiple sequence alignment for nucleotide and amino-acid sequences. *J Comput Biol* 22(5):377–386, 2015 (online Dec 2014). doi:10.1089/cmb.2014.0156 (verified).
  - **How it merges:**
    - A guide tree is split into subsets of at most 200 sequences, and each subset is aligned with MAFFT.
    - A spanning tree is built over the subsets. For each spanning-tree edge (v,w), the two subset alignments are aligned with OPAL. The merger "must not change the alignments on the type 1 subalignments", so the input columns stay fixed.
    - The overlapping pairwise ("type 2") alignments are then combined by a transitivity merge. This is an equivalence-closure over letters that share a column, and the result does not depend on the order of edges.
    - PASTA runs 3 iterations by default.
  - **Default merger:** The paper's methods name OPAL. The PASTA README says "the default is OPAL for dna and Muscle for protein" (switch with `--merger`).
  - **OPAL vs MUSCLE:** The PASTA paper does **not** compare them as mergers. It only says it ran SATé-II with MUSCLE as the merger "due to the high computational costs of running OPAL on large datasets."
  - **Accuracy:** No accuracy numbers are reported for the merge step alone.
- **Liu K, Warnow TJ, Holder MT, Nelesen SM, Yu J, Stamatakis AP, Linder CR.** SATé-II: Very fast and accurate simultaneous estimation of multiple sequence alignments and phylogenetic trees. *Syst Biol* 61(1):90–106, 2012. doi:10.1093/sysbio/syr095 (verified).
  - SATé uses tree-based decomposition, aligns subsets, then merges them with MUSCLE profile-profile (the SATé software also offers OPAL).
  - Warnow-group slides describe "using Opal instead of Muscle to merge" in the "next SATé" design. **[unverified:** I found no OPAL-vs-MUSCLE merger comparison in the SATé-II paper itself; the full text is not open access, so check it directly.**]**

### 2. OPAL and exact "aligning alignments"
- **Kececioglu J, Zhang W.** Aligning alignments. *CPM 1998*, LNCS 1448:189–208. doi:10.1007/BFb0030790 (verified).
  - Introduced optimistic and pessimistic gap counts for aligning alignments or profiles under affine gaps.
  - Gave a polynomial exact algorithm for aligning a single sequence against an alignment.
  - Conjectured that alignment-vs-alignment is NP-complete.
- **Ma B, Wang Z, Zhang K.** Alignment between two multiple alignments. *CPM 2003*, LNCS 2676:254–265. doi:10.1007/3-540-44888-8_19 (verified). Gives an independent NP-completeness proof.
- **Kececioglu J, Starrett D.** Aligning alignments exactly. *RECOMB 2004*, pp. 85–96. doi:10.1145/974614.974626 (verified).
  - **Problem:** Align the columns of A against the columns of B. Input columns are fixed; the alignment only substitutes, inserts or deletes whole columns. The objective is sum-of-pairs with "linear" gap costs, which in their terms means affine (γ + λ·length).
  - **Complexity:** NP-complete when gap counts are exact, and polynomial-time without them. So the hardness comes from counting gap openings exactly.
  - **Algorithm:** An exact DP over (i, j, shape), where a "shape" is an ordered partition of rows. It is solved as a shortest-path problem with bound pruning (using optimistic/pessimistic gap counts) and dominance pruning, and it has a linear-space version.
  - **Reported speed (talk slides):** real benchmarks (25 sequences, length 900) in about 1 s; simulated (200 sequences × 1000 columns) in about 3 min.
- **Wheeler TJ, Kececioglu JD.** Multiple alignment by aligning alignments. *Bioinformatics* 23(13):i559–i568, 2007 (ISMB). doi:10.1093/bioinformatics/btm226 (verified).
  - Opal builds the MSA by form-and-polish, merging subalignments with the exact aligning-alignments algorithm (SP objective, affine gaps counted exactly).
  - It uses no consistency and no external evidence.
  - The biggest gains came from distance estimation and 3-cut polishing, not from the exact merge.
  - Matches top tools on benchmarks "without employing alignment consistency."
  - The older AlignAlign software page says Opal supersedes it.

### 3. MUSCLE
- **Edgar RC.** MUSCLE: multiple sequence alignment with high accuracy and high throughput. *Nucleic Acids Res* 32(5):1792–1797, 2004. doi:10.1093/nar/gkh340 (verified).
  - Progressive profile-profile DP using a log-expectation profile score and position-specific affine gaps.
  - v3.x has `-profile -in1 A -in2 B`, which keeps the columns of each input MSA intact.
  - It is a heuristic score with no consistency and no external evidence.
- **Edgar RC.** Muscle5: High-accuracy alignment ensembles enable unbiased assessments of sequence homology and phylogeny. *Nat Commun* 13:6968, 2022. doi:10.1038/s41467-022-34630-w (verified).
  - The v5 manual's command list (align, super5, the ensemble tools, cloak) has **no profile-profile command**.
  - Internally, v5 "PPP" merges profiles with probabilistic consistency (HMM posteriors), but this is not exposed as a user-facing "merge two MSAs" mode. **[I did not check the v5 source; treat "removed" as likely but unconfirmed.]**
  - PASTA's MUSCLE merger refers to MUSCLE 3.8.

### 4. MAFFT
- **Katoh K, Frith MC.** Adding unaligned sequences into an existing alignment using MAFFT and LAST. *Bioinformatics* 28(23):3144–3146, 2012. doi:10.1093/bioinformatics/bts578 (verified).
  - Covers only `--add` and `--addfragments`, where the existing MSA is fixed and placement uses an inferred phylogeny.
  - It does **not** describe `--merge`, so this is the wrong citation for `--merge`.
- **Katoh K, Standley DM.** MAFFT multiple sequence alignment software version 7. *Mol Biol Evol* 30(4):772–780, 2013. doi:10.1093/molbev/mst010 (verified).
  - `mafft-profile` turns two MSAs into profiles and aligns them. It assumes the two groups are phylogenetically separate and is warned against for adding sequences.
  - `--addprofile aln1 aln2` requires aln1 to be monophyletic within aln2's tree.
  - The paper does not describe `--merge`.
- **`--merge` (documentation only, mafft.cbrc.jp/alignment/software/merge.html):**
  - Merges two or more sub-MSAs; "each sub-MSA is preserved."
  - Each sub-MSA is forced to be a monophyletic cluster of the guide tree (since v7.043), and progressive alignment then follows that tree. It works with L-INS-i (`--localpair`).
  - The docs say `--merge` "completely covers" `mafft-profile`, which "will be deleted."
  - No peer-reviewed algorithm description exists; this is the "MAFFT-merge" baseline used by Warnow's group.

### 5. T-Coffee and M-Coffee
- **Notredame C, Higgins DG, Heringa J.** T-Coffee: a novel method for fast and accurate MSA. *J Mol Biol* 302:205–217, 2000. doi:10.1006/jmbi.2000.4042 (verified). Progressive alignment over a library extended by consistency (triplet transitivity).
- **Wallace IM, O'Sullivan O, Higgins DG, Notredame C.** M-Coffee: combining multiple sequence alignment methods with T-Coffee. *Nucleic Acids Res* 34(6):1692–1699, 2006. doi:10.1093/nar/gkl091 (verified).
  - Turns several MSAs of the **same** sequences (from different aligners) into a weighted primary library, extends it by consistency, and re-aligns progressively.
  - This is consensus or meta-alignment, not merging disjoint alignments.
  - Input columns are not fixed; the output can differ from every input.
- **T-Coffee profile mode:** `-profile` / `-profile1 -profile2`, with `-profile_mode` defaulting to `cw_profile_profile` and `-profile_comparison` to `full50`, per the man page. **[The exact consistency treatment of profiles is unverified.]**
- Zaharias et al. (2023, item 6) used T-Coffee as one of only two available disjoint-alignment mergers. It performed poorly on nucleotide data.

### 6. MAGUS / GCM and the MWT-AM problem
- **Smirnov V, Warnow T.** MAGUS: Multiple sequence Alignment using Graph clUStering. *Bioinformatics* 37(12):1666–1672, 2021. doi:10.1093/bioinformatics/btaa992 (verified).
  - **Inputs:** Subset ("constraint") alignments are treated as absolute constraints, so input columns are fixed.
  - **External evidence:** 10 "backbone" alignments of 200 sequences each, built with MAFFT-linsi on random samples from every subset. Two columns are joined by an edge weighted by how many backbones align their letters.
  - **Steps:** MCL clustering (inflation 4) of the column graph, then resolution of conflicts (two columns from the same subset, crossing clusters) into a valid trace.
  - **Defaults:** 25 subsets, one iteration.
  - **Accuracy:** More accurate than PASTA on 9 of 10 1000-sequence model conditions (tied on 1), and faster.
  - The paper calls MSA merging a "straightforward generalization of the well-studied problem of aligning two alignments (Wheeler & Kececioglu 2007)."
- **Zaharias P, Smirnov V, Warnow T.** Large-scale multiple sequence alignment and the Maximum Weight Trace Alignment Merging problem. *IEEE/ACM TCBB* 20(3):1700–1712, 2023. doi:10.1109/TCBB.2022.3191848 (verified).
  - Defines **MWT-AM**: given disjoint MSAs A1..Ak and weights w(x,y) ≥ 0 on column pairs, find a trace (a merged MSA that induces each Ai) of maximum total weight.
  - **Theorem 1:** "MWT-AM is NP-hard but can be solved in O(2^B L^B B^2) for B alignments of total length L." For **B = 2** this is a polynomial (quadratic) DP. The paper does not develop the B = 2 case.
  - Weights come from MAGUS backbones. Variants: MCL, MRL-MCL, region growing (RG / RG-FAST), and Fiduccia-Mattheyses (FM), each followed by a greedy "optimizer". **GCM(fm+opt)** gives the best MWT-AM scores.
  - Reported findings:
    - MWT-AM scores correlate well with SP accuracy.
    - GCM(fm+opt) beats GCM-default in both MWT-AM score and SP score on every large protein dataset.
    - MAFFT-merge failed on the 88k and 93k-sequence sets.
    - GCM variants solve MWT-AM better than MAFFT-merge and T-Coffee.
  - Explicitly contrasts PASTA's pairwise OPAL/MUSCLE merging plus transitivity with GCM's all-at-once merge.
- Related, verified: **Smirnov V.** Recursive MAGUS. *PLoS Comput Biol* 17(10):e1008950, 2021. doi:10.1371/journal.pcbi.1008950. **Shen C, Zaharias P, Warnow T.** MAGUS+eHMMs. *Bioinformatics* 38(4):918–924, 2022. doi:10.1093/bioinformatics/btab788.

### 7. WITCH-NG and EMMA
- **Liu B, Warnow T.** WITCH-NG: efficient and accurate alignment of datasets with sequence length heterogeneity. *Bioinform Adv* 3(1):vbad024, 2023. doi:10.1093/bioadv/vbad024 (verified). This is the closest prior art.
  - Adding one query to a fixed backbone is framed explicitly as **the trivial two-alignment case of MWT-AM**: "solvable in O(mn) time and O(mn) space by a simple DP algorithm."
  - Weights come from the top-k (default 10) HMMs of the ensemble (eHMM), weighted by each HMM's support. Pairs no selected HMM aligns are set to −∞.
  - It runs Smith-Waterman with **zero gap penalty**, replacing WITCH's MCL + A* heuristic.
  - **Accuracy:** Identical to WITCH (within 0.001 SPFN/SPFP on biological data); for example, HomFam error is 0.225 for WITCH, 0.226 for WITCH-NG and 0.229 for UPP.
  - **Speed:** At least 5× faster on 7 of 10 HomFam datasets.
- **Shen C, Liu B, Williams KM, Warnow T.** EMMA: a new method for computing multiple sequence alignments given a constraint subset alignment. *Algorithms Mol Biol* 18:21, 2023. doi:10.1186/s13015-023-00247-x (verified).
  - Decomposes the constraint (backbone) alignment and its tree into small subsets, assigns queries to subsets, and runs MAFFT-linsi `--add` on each subproblem (at most about 500 sequences).
  - It then **merges the extended subalignments by transitivity through the shared backbone columns**. The constraint alignment is kept fixed, and "the transitivity merge never merges columns in the constraint alignment."
  - Insertion columns from different subproblems are never homologized; the paper notes a realignment stage as future work.
  - Scales to hundreds of thousands of sequences and is more accurate than WITCH-ng-add and MAFFT `--add` in many conditions.

### 8. Recent work, 2022–2026 (verified to exist)
- **Tseng Y-H, Walia S, Turakhia Y.** Ultrafast and ultralarge multiple sequence alignments using TWILIGHT. *Bioinformatics* 41(Suppl 1):i332–i341, 2025 (ISMB). doi:10.1093/bioinformatics/btaf212.
  - Progressive profile-profile Needleman-Wunsch with affine gaps (ClustalW-style profiles), TALCO tiling, adaptive X-drop banding and GPU support.
  - Gappy columns (above 95% gaps) are removed before profile alignment and added back afterwards.
  - Subtree MSAs are merged by progressive profile alignment or a transitivity merger. The authors report shorter alignments than PASTA's transitivity merger "with minimal impact on alignment quality."
  - Aligned 1M RNASim sequences in under 30 min using less than 16 GB.
- **Gudyś A, Zielezinski A, Notredame C, Deorowicz S.** Fast and accurate multiple-protein-sequence alignment at scale with FAMSA2. *Nat Biotechnol*, 2026. doi:10.1038/s41587-026-03095-3; preprint doi:10.1101/2025.07.15.664876.
  - Progressive alignment with a medoid guide tree.
  - The FAMSA tool (since v2.2.0) supports profile-profile alignment of two input FASTA MSAs while keeping each profile's columns. **[The paper's merge internals are unverified.]**
- **Sievers F et al.** Clustal Omega. *Mol Syst Biol* 7:539, 2011. doi:10.1038/msb.2011.75 (verified).
  - Progressive HMM-HMM alignment via HHalign.
  - `--profile1 --profile2` aligns two MSAs with each one's columns fixed.
- **Park M, Ivanovic S, Chu G, Shen C, Warnow T.** UPP2. *Bioinformatics* 39(1):btad007, 2023. doi:10.1093/bioinformatics/btad007 (verified). Adds sequences to a backbone; does not merge two MSAs.
- I found no 2024–2026 paper or preprint specifically on exact or MWT merging of **two** alignments as a pairwise merger inside PASTA.
  - HMMerge (Park & Warnow, *Bioinform Adv* 2023, vbad052) and the ARIES protein-language-model preprint (bioRxiv 2026) exist but address adding sequences or building MSAs from scratch.
  - I could not confirm a "learning-based profile-profile merger".

### Is the idea novel? ("Exact MWT DP for merging two alignments, with backbone-derived edge weights, used as PASTA's pairwise merger")
- **Not novel (the algorithm and the objective):**
  - Two-alignment MWT-AM is polynomial: Theorem 1 of Zaharias et al. 2023 with B = 2.
  - WITCH-NG already implements the O(mn) DP for the sequence-vs-alignment special case with ensemble-derived weights and zero gap cost.
  - Backbone-derived column-pair weights are exactly MAGUS/GCM Step 0.
  - Exact two-alignment merging inside PASTA already exists as OPAL, under an SP/affine objective rather than MWT. Extending WITCH-NG's DP from 1×m to m×n columns is a natural, small step.
- **Plausibly novel (the combination and the study):**
  - I found no paper that replaces PASTA's OPAL/MUSCLE pairwise merger with an MWT objective using external backbone evidence, or that compares that against OPAL/MUSCLE and GCM on the same decomposition.
  - Also untested: pairwise-exact MWT plus transitivity versus GCM's k-way heuristic for NP-hard MWT-AM, i.e., whether exactness on pairs beats a heuristic on k-way.
- **Caveats to state in the report:**
  - Being exact for each pair does not make the overall k-way merge optimal: transitivity can propagate pairwise errors, and the MWT-AM optimum is NP-hard for k ≥ 3.
  - Accuracy will depend mostly on the weights (number and size of backbones), as MAGUS showed.
  - Unsupported pairs must be handled deliberately: weight 0 (allowed) or −∞ (forbidden, as in WITCH-NG). With zero gap cost the DP is a weighted LCS / global trace with no gap penalties, which differs from OPAL's affine SP objective.
  - Strongest baselines: PASTA with OPAL, PASTA with MUSCLE, MAGUS GCM-default, GCM(fm+opt) and MAFFT `--merge`.

### Sources
- [PASTA (PMC4424971)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4424971/)
- [PASTA README](https://github.com/smirarab/pasta)
- [MAGUS (PMC8289385)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8289385/)
- [Zaharias et al. TCBB author version](https://par.nsf.gov/servlets/purl/10415277)
- [WITCH-NG](https://academic.oup.com/bioinformaticsadvances/article/3/1/vbad024/7068428)
- [EMMA (PMC10704716)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10704716/)
- [TWILIGHT (PMC12261412)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12261412/)
- [MAFFT merge documentation](https://mafft.cbrc.jp/alignment/software/merge.html)
- [MAFFT v7 (PMC3603318)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3603318/)
- [MUSCLE5 command list](https://drive5.com/muscle5/manual/commands.html)
- [MUSCLE v3 profile-profile](https://drive5.com/muscle/manual/profprof.html)
- [Kececioglu–Starrett RECOMB 2004 talk slides](https://www2.cs.arizona.edu/classes/cs550/fall24/handouts/recomb2004-talk.pdf)
- [AlignAlign page](https://alignalign.cs.arizona.edu/)
- [Ma–Wang–Zhang CPM 2003](https://link.springer.com/chapter/10.1007/3-540-44888-8_19)
- [M-Coffee (PMC1410914)](https://pmc.ncbi.nlm.nih.gov/articles/PMC1410914/)
- [FAMSA repository](https://github.com/refresh-bio/FAMSA)
- [T-Coffee man page](https://manpages.org/t-coffee)
- Crossref API, for all DOIs