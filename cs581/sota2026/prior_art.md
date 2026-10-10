# Prior art: do 2023–2026 benchmarks pit the fast aligners against MAGUS/PASTA/UPP2 on nucleotide data?

Checked 2026-10-10 (~35 min). This file extends `../literature/sota_review.md` (sections 2–3). DOIs marked [CR] were verified today through api.crossref.org. Those marked [EPMC] were checked through Europe PMC. Unmarked DOIs are reused from sota_review.md. Full texts read: TWILIGHT, Muscle5, learnMSA2, UPP2, Nute & Warnow 2016 and PASTA 2015 (Europe PMC/PMC), Gupta 2021 (OUP). **The FAMSA2 full text could not be fetched** (bioRxiv returned HTTP 429 and Nature is paywalled). Its entry relies on the abstract, the bioRxiv metadata and the GitHub README.

## Short answer
**No.** No 2023–2026 paper runs TWILIGHT, FAMSA2, Muscle5/Super5, Clustal Omega, the MAFFT modes and learnMSA2 against MAGUS/PASTA/WITCH/UPP2 on the MAGUS nucleotide suite (ROSE 1000, RNASim, CRW 16S).
- **TWILIGHT** is the only post-2023 paper with any of the "fast" tools against MAGUS/PASTA on nucleotides. It covers RNASim and AliSim only, with no ROSE, no CRW and no FAMSA.
- **UPP2** is the only 2023+ paper with Clustal Omega and MAFFT against MAGUS/PASTA on ROSE, CRW 16S and RNASim. It used MUSCLE **3.8**, not Muscle5.
- **FAMSA2 and learnMSA2** are protein-only.
- **Muscle5's** only nucleotide evidence is the small-RNA Bralibase set, and it never compared against MAGUS.

## Paper table

| # | Paper (DOI) | Nucleotide data | Competitors relevant here | Nucleotide result (headline) |
|---|---|---|---|---|
| 1 | Tseng, Walia, Turakhia. TWILIGHT. *Bioinformatics* 41(S1):i332–i341, 2025. doi:10.1093/bioinformatics/btaf212 [CR] | RNASim 10k→1M (~1,500 nt); AliSim 10k long- and short-branch (up to 1 Mb length); 8.1M SARS-CoV-2 (no reference) | Clustal Omega, MAFFT (mode not given in main text), Muscle5, regressive T-Coffee, PASTA, recursive MAGUS. **No FAMSA, UPP2, WITCH or learnMSA** | RNASim: TWILIGHT 6–8% error; MAGUS and Muscle5 ~8% at 10k; PASTA ~3 pts worse; Clustal/MAFFT/T-Coffee >20% |
| 2 | Gudyś, Zieleziński, Notredame, Deorowicz. FAMSA2. *Nat Biotechnol* 2026-04-14. doi:10.1038/s41587-026-03095-3 [CR]; bioRxiv 10.1101/2025.07.15.664876 [CR] | **None found.** Published title: "multiple-**protein**-sequence alignment" | Clustal Ω 1.2.4, MAFFT **DPPartTree** 7.526, Kalign 3.4.1, Muscle5 5.3 **-super5**, regressive T-Coffee, FAMSA 1.1 (from README). **MAGUS, TWILIGHT and learnMSA not listed** | n/a |
| 3 | Edgar. Muscle5. *Nat Commun* 13:6968, 2022. doi:10.1038/s41467-022-34630-w [CR] | Bralibase (small structural RNA sets) | Clustal Omega, MAFFT (mode not stated in main text), ProbCons | "higher than state-of-the-art": **numbers only in Fig. 2 and Supp. Table S2** |
| 4 | Becker & Stanke. learnMSA2. *Bioinformatics* 40(S2):ii79–ii86, 2024. doi:10.1093/bioinformatics/btae381 [CR] | **None.** "Currently, learnMSA supports only protein alignments" | Clustal Ω, MAFFT-sparsecore, Muscle5, regressive T-Coffee, FAMSA, MAGUS (HomFam, protein) | n/a. CPU runtime given (see below) |
| 5 | Park, Ivanovic, Chu, Shen, Warnow. UPP2. *Bioinformatics* 39(1):btad007, 2023. doi:10.1093/bioinformatics/btad007 [EPMC] | ROSE 1000 S/M/L (+HF), RNASim1000 (+HF), 10 CRW sets (16S.3, 16S.T, 16S.B.ALL + 7 small 16S/23S/5S) | MAGUS, PASTA, MAFFT (auto on large; linsi/ginsi/xinsi/qinsi on small), Clustal Ω 1.2.4, **MUSCLE 3.8.31**, regressive T-Coffee | UPP2 ≈ MAGUS best on 16S; Clustal Ω worst. **All numbers in figures** |
| 6 | Gupta, Zaharias, Warnow. BAli-Phy as phylogeny-aware aligner. *Bioinformatics* 37(24):4677–4683, 2021. doi:10.1093/bioinformatics/btab555 [CR] | ROSE (SATé) 25–1000 seqs, INDELible 25, RNASim1000 ×2 | MAFFT L-INS-i (`--maxiterate 1000 --localpair`), PRANK, PASTA | BAli-Phy (fixed tree) beats L-INS-i significantly at 25, 100 and 500 seqs |
| 7 | Nute & Warnow. Scaling statistical MSA to large datasets. *BMC Genomics* 17(S10):764, 2016. doi:10.1186/s12864-016-3101-8 [CR] | ROSE 1000 L1/M1/S1, INDELible M2, RNASim 1000 (×10 reps) | PASTA(default), **PASTA+BAli-Phy subsets**, PASTA+MAFFT-L (subset 100), MAFFT L-INS-i, MAFFT default | PASTA with a non-MAFFT subset aligner: best TC everywhere (Table 2) |
| 8 | Collins & Warnow. PASTA for proteins. *Bioinformatics* 34(22):3939–3941, 2018. doi:10.1093/bioinformatics/bty495 [CR] | **None** (BAliBASE 4, 224 sets with ≥50 seqs) | PASTA with MAFFT L-INS-i, G-INS-i, MAFFT-homologs, CONTRAlign, ProbCons as subset aligners | Protein-only. **Confirms PASTA with ProbCons subsets** |
| 9 | Mirarab et al. PASTA. *J Comput Biol* 22(5), 2015. doi:10.1089/cmb.2014.0156 | RNASim 10k–200k, 16S.3/16S.T/16S.B.ALL, ROSE 1000 (supp. E.2) | SATé-II, MAFFT, Clustal Ω, MUSCLE (v3) | Only MAFFT L-INS-i subsets and an Opal merger in the main text. **No alternative subset aligners tested** |
| 10 | Lassmann. Kalign 3. *Bioinformatics* 36(6):1928–1929, 2020. doi:10.1093/bioinformatics/btz795 [CR] | Claims "high alignment accuracy on both protein and nucleotide" (datasets not checked) | (unverified) | Abstract-level claim only |

## What each found on nucleotide data

**TWILIGHT (2025).** This is the key paper. Full text read via Europe PMC.
- **Datasets.** RNASim subsamples of 10k, 100k, 200k and 1M, following Recursive MAGUS; AliSim 10k sets; SARS-CoV-2. **No ROSE, CRW 16S or BAliBASE.**
- **Guide trees.** Clustal Ω, PASTA, regressive T-Coffee and MAGUS (recursive) were given the true tree and run for one iteration. MAFFT and Muscle5 built their own trees.
- **Missing details.** Tool versions and the MAFFT mode are only in the Supplementary PDF (not read; **unverified**). At ≥10k sequences MAFFT was presumably in auto/FFT-NS mode (**unverified**).
- **Errors** (FastSP = mean of SPFN and SPFP):
  - TWILIGHT reaches 6.48% error at 1M.
  - At 10k, MAGUS and Muscle5 are "around 8%", and PASTA is ~3 pts higher.
  - Fig. 2 legend: TWILIGHT 6–8%; Muscle5, MAGUS and PASTA 8–13%; Clustal Ω, MAFFT and T-Coffee >20%. **Per-size values are only in Fig. 2a.**
- **Speed.**
  - TWILIGHT on CPU is ~13× faster than MAFFT at 100k sequences and 47–176× faster than the other tools.
  - 1M sequences take 32 min on CPU or 18 min on 1 GPU.
- **Failures.** MAGUS did not finish 200k in 24 h. Muscle5 "did not scale" past 10k. All other tools failed at 1M.
- **Output size at 100k** (true alignment 1.4 GB): TWILIGHT 1.6 GB, MAFFT 1.9, T-Coffee 2.2, PASTA 18, MAGUS 21.9.
- **Iterative mode on AliSim long-branch 10k** (Fig. 4 only): MAGUS 8.3%, TWILIGHT ~9%, PASTA 12%. TWILIGHT was 4.6× faster than MAGUS.
- **Table 1 (TWILIGHT alignment error by guide tree).**
  - RNASim: true tree 7.83%, MAFFT tree 8.68%, PartTree 11.58%.
  - AliSim long-branch: 8.41%, 8.13%, 13.42%.
  - These rows change only TWILIGHT's guide tree. **They are not MAFFT-PartTree alignments.**

**FAMSA2 (2025/2026).**
- Protein-only by title and README. Benchmarks: extHomFam v37, afdb_clusters, plus simulated MSAs.
- Competitor list (from README): Clustal Ω, MAFFT `--dpparttree`, Kalign3, Muscle5 `-super5`, regressive T-Coffee, FAMSA1.
- **MAGUS, TWILIGHT and learnMSA are absent.** Whether the paper mentions nucleotides at all is **unverified** (full text not read).

**Muscle5 (2022).**
- The only nucleotide benchmark is Bralibase RNA: small sets, each well under 200 sequences.
- Muscle5 beats Clustal Ω and MAFFT there, but **the numbers are in Fig. 2 and Supp. Table S2 only**, and the MAFFT mode is not stated in the main text (**unverified**).
- The "10,000 sequences" result (59% vs 52% for Clustal Ω and 47% for MAFFT, columns correct) is **protein**.
- No ROSE, RNASim or CRW, and **no MAGUS or PASTA**. Super5 was never evaluated on nucleotides here.

**learnMSA2 (2024).**
- Protein-only: "Currently, learnMSA supports only protein alignments. We may extend it in the future to align DNA or RNA".
- CPU runtime on HomFam with protein language-model (pLM) embeddings is ~3 h per family. That is 15× slower than on GPU (0.2 h) and 16× slower than MAGUS (0.19 h), the slowest CPU method.
- Hardware: A100 80 GB GPUs; 32 cores and 250 GB per job.
- It needs a few hundred sequences or more; "very small numbers of sequences" give "very poor alignments".

**UPP2 (2023).** This is the closest thing to a broad nucleotide comparison, but it comes from the Warnow group.
- **ROSE 1000 (14 conditions, full-length and HF).** UPP2 ≈ UPP(MAGUS) > MAGUS > PASTA ≈ MAFFT > MUSCLE3 > Clustal Ω ≈ T-Coffee.
- **Large CRW 16S.**
  - UPP2 and MAGUS tie for most accurate. PASTA matches them on 16S.B.ALL and 16S.3 but not on 16S.T.
  - Clustal Ω has the highest error on all three. MAFFT-auto beats MUSCLE3 and T-Coffee on 16S.B.ALL and 16S.T.
  - MAFFT-xinsi, qinsi and ginsi failed on these sets.
  - On 16S.B.ALL, UPP2 took ~10 h versus ~110 h for UPP(MAGUS)+adj.
- **7 small CRW sets.** UPP2, UPP(MAGUS) and MAGUS are best, then PASTA, MAFFT-linsi and MAFFT-ginsi "fairly closely behind". Among the MAFFT variants, ginsi was slightly ahead of linsi. MUSCLE, T-Coffee and Clustal Ω are less accurate.
- **All numbers are in figures and supplementary figures (S28–S30).**

**Gupta et al. (2021).** Small-n nucleotide data with MAFFT L-INS-i as the baseline.
- BAli-Phy with a fixed estimated tree is significantly more accurate than L-INS-i at 25 and 100 sequences, and at 500 sequences (500M1 and 500M1-P).
- On 100M1 (ROSE), L-INS-i has ~25% error.
- Caveat from the authors: most MCMC effective sample sizes (ESS) are below 200, i.e. no evidence of convergence.

**Nute & Warnow (2016).** Pre-2023, but it is the main "PASTA with a non-MAFFT subset aligner" result on the MAGUS-style simulations (Table 2, 1000 sequences, precision/recall).
- **INDELible M2:** P+BAli-Phy 98.7/98.7 vs P(default) 95.1/94.6 vs MAFFT-L 80.2/75.0.
- **RNASim:** P+BAli-Phy 92.1/92.1 vs P(default) 90.3/90.4 vs MAFFT-L 91.8/91.5.
- **ROSE M1:** 79.8/79.6 vs 79.7/79.0 vs 74.9/63.3.
- **ROSE S1:** P(default) is best (85.3/85.1 vs 84.3 for P+BAli-Phy).
- **TC** is the largest gain everywhere, e.g. ROSE L1: 33.2% vs 15.9%.
- Cost: the BAli-Phy subsets needed 24 h on 32 cores each.

**MAGUS with a non-MAFFT base aligner.** **No paper found.** The MAGUS paper lists PROMALS, ProbCons and BAli-Phy as future work. Recursive MAGUS only replaces MAFFT with MAGUS on large subsets.
- The PASTA analogues exist: BAli-Phy subsets (Nute 2016, nucleotide) and ProbCons, CONTRAlign and G-INS-i subsets (Collins 2018, protein).

**Beating MAFFT L-INS-i on small (≤200) nucleotide sets.** Evidence:
- BAli-Phy with a fixed tree (Gupta 2021; ROSE/INDELible, 25–100 seqs).
- Muscle5 on Bralibase (versus "MAFFT", mode **unverified**).
- MAFFT-ginsi slightly beats linsi on the small CRW sets (UPP2 supplement).
- MAGUS and UPP2 are marginally better on small CRW (figures only). MAGUS on small inputs is essentially MAFFT-linsi (**unverified** whether it decomposes at ≤200).
- No 2023–2026 paper tests TWILIGHT, FAMSA2 or Super5 at this size on nucleotides.

## Gaps (what nobody has measured)
1. **TWILIGHT against the MAGUS nucleotide suite.** TWILIGHT has not been run on ROSE 1000 (S/M/L) or on CRW 16S.3, 16S.T, 16S.B.ALL or 16S.M. It has not been compared with UPP2 or WITCH on any nucleotide data.
2. **FAMSA2 on nucleotides.** No reference-based nucleotide benchmark of FAMSA2 (or FAMSA1) exists. FAMSA2 has never been compared to MAGUS, PASTA or TWILIGHT.
3. **Muscle5 and Super5 on large nucleotide sets.** Neither has been benchmarked on ROSE, RNASim ≥1k or CRW 16S against MAGUS. The only point is TWILIGHT's ~8% at RNASim 10k.
4. **MAFFT modes on the same data.** No paper reports PartTree/DPPartTree, auto and L-INS-i side by side on RNASim or 16S against MAGUS. TWILIGHT's PartTree rows are guide trees only.
5. **learnMSA2 has no nucleotide mode.** BAliBASE comparisons of MAGUS against Muscle5, FAMSA2 and learnMSA2 under the MAGUS filtering also do not exist.
6. **MAGUS with alternative base aligners.** MAGUS has never been run with Muscle5, FAMSA2 or BAli-Phy subsets, and its GCM merger has never been fed TWILIGHT subalignments. This is an open, cheap experiment: the `-s` and `-b` hooks already exist.
7. **No independent nucleotide benchmark.** Every nucleotide comparison above was run by a method's own developers. nf-core/multiplesequencealign (Santus 2025, doi:10.1093/nargab/lqaf104) bundles the tools but publishes no nucleotide ranking.
