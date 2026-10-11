# Prior art: sample complexity and missing data for MSC summary methods

Notation: f is the shortest internal branch length in coalescent units, n is the number of taxa, m is the number of genes/loci, and k is the sequence length per locus.
Verification: DOIs were checked against Crossref unless marked otherwise. Scaling statements come from arXiv/ar5iv full text where noted. "(unverified)" means I did not see the primary text.

## 1. Sample-complexity results under the MSC

| Paper | DOI / arXiv | Method | What is proven (m as a function of f, n) |
|---|---|---|---|
| Shekhar, Roch, Mirarab 2018, IEEE/ACM TCBB 15(5) | 10.1109/TCBB.2017.2757930; arXiv:1704.06831 | ASTRAL* (exact max-quartet-support), true gene trees | **Upper bound:** m > (9/2) log(4·C(n,4)/ε) / (1−e^{−f})², i.e. for small f, m > 20 log(n/(6ε)) / f², so **O(f⁻² log n)** (Thm 2.1). **Lower bound for ASTRAL\*:** Θ(f⁻²) is tight for n=4 (Thm 2.2). For general n there is a tree on which ASTRAL\* fails w.p. ≥ 1−ρ when m ≤ (a/5)(log n)/f² (Thm 2.3), so **ASTRAL\* needs Θ(f⁻² log n)** in the worst case. ASTRAL-I/II with constrained search: the paper suggests Θ(f^−(n−3)) via bipartition coverage (very conservative). Under a toy i.i.d. quartet-error model with e₀ < 2/3, the bounds still hold. |
| Mossel & Roch 2010, IEEE/ACM TCBB 7(1):166–171 | 10.1109/TCBB.2008.66; arXiv:0710.0262 | GLASS (min pairwise coalescence time across genes; STAR-like), uses gene-tree branch lengths | Consistency, with tolerance to moderate estimation error. As summarized by Dasarathy et al. 2015 and Shekhar et al. 2018, GLASS needs **m = O(f⁻¹)** genes when gene trees and branch lengths are exact. With sequences, k = Ω(f⁻²) is also needed, so total data is about f⁻³. The exact log n factor in the original theorem is (unverified). |
| STEAC / STAR (Liu et al. 2009) | — | STEAC (average coalescence times) | Per Dasarathy et al. 2015 and Shekhar et al. 2018, STEAC needs m ∝ **f⁻²** (true gene trees with branch lengths). STAR has no stated sample-complexity result (Shekhar et al.: "no separate result"). |
| Dasarathy, Nowak, Roch 2015, IEEE/ACM TCBB 12(2) (ISIT 2014 version: 10.1109/ISIT.2014.6875191) | 10.1109/TCBB.2014.2361685; arXiv:1404.7055 | METAL (NJ on concatenated log-det/JC distances) | **Upper bound:** for any k ≥ 1, m ≥ C·f⁻²·log(n/ε) suffices (Thm 4; the constant depends on mutation-rate bounds and depth Δ). So total data m·k = O(f⁻² log n) even with constant-length genes. **Lower bounds for any method:** m = Ω(f⁻¹) regardless of k (the branch leaves no trace otherwise), and m·k = Ω(f⁻²) (Steel–Székely), so m ≥ C·max{f⁻¹, f⁻²/k}. |
| Mossel & Roch 2017, Ann. Appl. Probab. 27(5) | 10.1214/16-AAP1273; arXiv:1504.05289 | Information-theoretic, distance-based / sequence data | To detect a branch of length f you need **m = Θ(1/(f²√k))** loci (sparse signal detection). This sharpens the m–k trade-off left open by Dasarathy et al. |
| Roch 2018, RECOMB-CG 2018, LNCS 11183 | 10.1007/978-3-030-00834-5_11; arXiv:1812.08357 | Internode distance (NJst/ASTRID) | **No sample-complexity bound is proven.** The paper proves a **lower bound on worst-case variance**: there are trees with Var[δ̂_int(ℓ,ℓ′)] ≥ C·d_S(ℓ,ℓ′)/m and max Var ≥ C·n/m. The author notes this suggests NJ over all distances needs m that is at least linear in n, compared with log n for ASTRAL. Fast-converging methods (distances within depth about log n) might match ASTRAL. Correlations between distance entries could lower the requirement. Tight bounds are stated as an open problem. |
| Allman, Degnan, Rhodes 2018, IEEE/ACM TCBB 15(1):337–342 | 10.1109/TCBB.2016.2604812; arXiv:1604.05364 | NJst / unrooted STAR | **Consistency only** (arbitrary n). No sample complexity. |
| Roch & Warnow 2015, Syst Biol 64(4):663 | 10.1093/sysbio/syv016 | Summary methods with estimated gene trees | Studies statistical guarantees under gene-tree estimation error. Only the abstract was seen: theorem details are (unverified). The bounded-sequence-length inconsistency result is usually credited to Roch, Nute & Warnow 2019 (next row). |
| Roch, Nute, Warnow 2019, Syst Biol 68(2):281 | 10.1093/sysbio/syy061; arXiv:1803.02800 | Fully partitioned ML; topology-based summary methods (e.g. ASTRAL, NJst) | **Inconsistency** when sequence length per locus is bounded, even with highly constrained heterogeneity, driven by long-branch attraction. This is a negative result, not a sample-complexity bound. |
| Hill, Legried, Roch 2022, Ann. Appl. Probab. 32(6) | 10.1214/22-AAP1799; arXiv:2007.06697 | ASTRAL-one under DLCoal (ILS + GDL) | m ≥ C′ f⁻² · e^{C\|μ−λ\|Δ} / (1−(λ/μ ∧ μ/λ))^C · log(n/ε), i.e. still **O(f⁻² log n)** with GDL-dependent constants. ASTRAL-multi: consistency only. No lower bound (listed as open). |
| Legried, Molloy, Warnow, Roch 2021, J Comput Biol 28(5) | 10.1089/cmb.2020.0424; bioRxiv 10.1101/821439 | ASTRAL-multi under GDL | Identifiability, plus **consistency** of ASTRAL-multi. No sample complexity. |
| Chan, Li, Scornavacca 2022, J Math Biol | 10.1007/s00285-022-01786-4 | Max-quartet-support (ASTRAL) | Error probability decays exponentially in m. Gives a closed form for n=4 and numerical bounds tighter in practice than Shekhar et al. Only the abstract was seen; the exact form is (unverified). |
| Roch & Snir 2013, J Comput Biol 20(2) | 10.1089/cmb.2012.0234 | Quartet-based, under lateral gene transfer (not ILS) | Sample-complexity analysis under an LGT model. Exact scaling is (unverified). |
| Hill & Roch 2025 preprint | bioRxiv 10.1101/2025.01.28.635331 | Any method, error-free gene trees, Gamma rate variation across loci | Information-theoretic lower bounds. The "impossibility zone" expands when rates vary across loci. Exact theorems are (unverified), and the journal version is (unverified). |

**Bottom line for Q1.**
- **ASTRAL\*:** the only summary method with matching upper and lower bounds, Θ(f⁻² log n).
- **Topology-only methods:** Shekhar et al. state f⁻² is "the best demonstrated" requirement. I found no published Ω(f⁻²) lower bound for *every* topology-based method. The heuristic is that quartet-frequency gaps are Θ(f), so Ω(f⁻²) should hold for n=4 by a standard testing argument, but no citation was found: treat as (unverified).
- **Methods using gene-tree branch lengths:** GLASS gets O(f⁻¹), which matches the any-method lower bound Ω(f⁻¹) from Dasarathy et al.
- **With sequence data:** METAL needs O(f⁻² log n) loci with constant k. The Mossel & Roch 2017 Θ(1/(f²√k)) result describes the m–k trade-off.
- **NJst/ASTRID:** no proven sample-complexity bound (Shekhar et al. conjecture O(f⁻²)). Roch 2018 gives only a variance lower bound, about n/m.

## 2. Missing data

| Paper | DOI / arXiv | Finding |
|---|---|---|
| Nute & Chou (+ Molloy, Warnow) 2017, RECOMB-CG, LNCS | 10.1007/978-3-319-67979-2_15 | Conference version claiming consistency of ASTRAL, ASTRID and NJst under i.i.d. taxon deletion (M_iid). |
| Nute, Chou, Molloy, Warnow 2018, BMC Genomics 19(Suppl 5):286 | 10.1186/s12864-018-4619-8 | Journal version, with the same claim plus a simulation study of missing-data models. |
| Correction, 2020, BMC Genomics 21 | 10.1186/s12864-020-6540-1 | Retracts the NJst/ASTRID part (Theorem 11) after Rhodes' counterexample: NJst and ASTRID are **not** consistent under MSC + M_iid and can be positively misleading. ASTRAL's consistency under M_iid stands. |
| Rhodes, Nute, Warnow 2020 | arXiv:2001.07844 | Gives the counterexample: expected internode distances under i.i.d. deletion converge to an additive matrix for a **wrong** tree. **Still only on arXiv** (v1 only, no journal-ref as of this search). The published record is the BMC Genomics correction above. |
| Morel, Williams, Stamatakis 2023, Bioinformatics 39(1):btac832 | 10.1093/bioinformatics/btac832 | **Asteroid**: the missing-data-corrected internode-distance method. Each gene's internode matrix D_k is compared with the **induced** species tree S\|L(G_k) through a per-gene length Σ 2^{−M_k(i,j)} D_k(i,j), summed over genes. The authors state they *prove* consistency under the MSC with any random taxon deletion independent across genes (proof in the supplement, which I did not read). Empirically it beats ASTRAL and ASTRID when more than 80% of data is missing. **Most directly relevant prior art for any "normalized internode distance" idea.** |
| Liu & Warnow 2023, Algorithms Mol Biol 18:6 | 10.1186/s13015-023-00230-6 | Weighted ASTRID (wASTRID): weights internode distances by gene-tree branch support. Does **not** fix missing data; the authors list robustness to missing data as future work, citing Asteroid's criterion. |
| Han & Molloy 2023, Genome Res 33(7):1042 (TREE-QMC) | 10.1101/gr.277629.122 | The n1/n2 normalizations reweight quartets involving *artificial taxa* during divide-and-conquer. They are **not** a missing-data correction; the TREE-QMC-v2 README reportedly warns against the normalization option when data are missing (unverified wording). Empirically robust to missing data. No consistency claim for missing data found. |
| Christensen, Molloy, Vachaspati, Warnow 2018, Algorithms Mol Biol 13:6 (OCTAL) | 10.1186/s13015-018-0124-5 | Completes incomplete gene trees optimally (minimum RF to a reference tree) in polynomial time. TRACTION (Christensen et al., WABI 2019; doi (unverified)) generalizes this to non-binary trees. |
| ASTRAL under missing data | (as in the Nute et al. rows) | ASTRAL remains consistent under M_iid (Nute et al. 2018, not retracted). I found no separate "taxon deletion" paper. |

I found no 2020–2026 paper other than Asteroid that proposes a normalization of internode distances for missing taxa, nor a proof that any ASTRID variant is consistent under missing data. Searches covered arXiv, bioRxiv and AMB, but the search was not exhaustive.

## 3. Empirical "genes needed" for ASTRAL vs ASTRID

- **Shekhar, Roch, Mirarab 2018** (10.1109/TCBB.2017.2757930) is the only direct measurement found.
  - **Setup:** n = 8 with caterpillar (anomaly zone), balanced and double-quartet trees. f has 9 values in [0.005, 0.1] and m goes up to 10⁵. There are 401 replicates per cell.
  - **Measurement:** the smallest m with ≤ 40/401 incorrect trees, i.e. error ε = 0.1 (P(correct) about 0.9, **not 0.95**), found by binary search and fit linearly against 1/f².
  - **ASTRAL-II results:** balanced 206 genes at f = 0.1 and 48,948 at f = 0.005. Caterpillar 255 and 59,528. Double-quartet 297 and 93,750.
  - **NJst (ASTRID) results:** also scales as about 1/f².
  - **ASTRAL vs ASTRID:** ASTRAL-II needed about 10–20% fewer genes on the caterpillar and balanced trees. NJst/ASTRID needed **fewer** genes on the double-quartet tree, which has long root branches, so ASTRAL is more sensitive to long branches.
  - **Error-free gene trees:** I believe the simulation used error-free gene trees (the abstract-level text suggests true MSC gene trees), but this is (unverified).
- **Vachaspati & Warnow 2015** (ASTRID, BMC Genomics 16(Suppl 10):S3, 10.1186/1471-2164-16-S10-S3) and **Mirarab & Warnow 2015** (ASTRAL-II, Bioinformatics 31(12):i44, doi 10.1093/bioinformatics/btv234, unverified by Crossref here) report accuracy vs. number of genes and ILS level. Neither reports a "genes needed for 95% recovery of a branch of length f" quantity (from memory: (unverified)).
- **Molloy & Warnow 2018** (Syst Biol 67(2):285, 10.1093/sysbio/syx077) is about gene filtering. It compares accuracy at fixed gene counts, not a per-branch sample-size threshold (unverified detail).
- **Chan et al. 2022** compares numerical sample-size bounds with simulations, for ASTRAL only (unverified detail).

**Gap:** no paper found measures m₀.₉₅(f) per internode for ASTRID vs ASTRAL, under missing data, or at n > 8. A per-branch empirical sample-complexity curve (ASTRID vs ASTRAL vs Asteroid, with and without M_iid deletion) appears to be open.
