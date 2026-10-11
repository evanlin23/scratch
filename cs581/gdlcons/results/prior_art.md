# Prior art: consistency of species-tree methods under GDL and DLCoal (search date 2026-10-10)

Legend: DOI "Y" means I resolved it on api.crossref.org (title and authors matched). [A] means I read it in the paper's text (the full text, or the abstract where I say so). [M] means memory or inference that I did not check against the text.

## Table

| paper | venue/year | DOI (verified) | method(s) | model | result | rooting/tagging assumptions |
|---|---|---|---|---|---|---|
| Zhang, Scornavacca, Molloy, Mirarab, "ASTRAL-Pro" | MBE 2020 (+ corrigendum to Def. 1) | 10.1093/molbev/msaa139 (Y) | ASTRAL-Pro (MLQST objective) | GDL (Arvestad et al. 2009 birth-death) | **Proved consistent** (Thm 2) [A] | **Correctly rooted and tagged** gene trees. "Partially correct rooting" is enough (Claim 1). Being consistent with imperfect rooting/tagging is only "suspected". DLCoal is explicitly left open [A] |
| Zhang & Mirarab, "ASTRAL-Pro 2" | Bioinformatics 2022 | 10.1093/bioinformatics/btac620 (Y) | ASTRAL-Pro 2 (placement-based) | — | No new theory (it is about algorithms and speed) [M] | Same objective as before [M] |
| ASTRAL-Pro3 (ASTER README only) | software, no paper found | — | A-Pro3 + CASTLES-Pro | — | No paper or theory found [A: search only] | — |
| Legried, Molloy, Warnow, Roch | JCB 2021 | 10.1089/cmb.2020.0424 (Y) | ASTRAL-one, ASTRAL-multi | GDL (Arvestad) | **Proved consistent** [A, via the DISCO and wQFM-DISCO texts. M for details] | **None**: unrooted, untagged MUL-trees (copies treated as alleles) |
| Markin & Eulenstein | Bioinformatics 2021 (arXiv 2004.04299) | 10.1093/bioinformatics/btab414 (Y) | Quartet methods, ASTRAL-one/ASTRAL-multi | **DLCoal** | **Proved consistent**. Key lemma: for a uniformly sampled gene per species, P(ab\|cd) > P(ac\|bd) = P(ad\|bc) [A, via wQFM-DISCO's restatement] | None (no rooting, no tagging) |
| Hill, Legried, Roch | Ann. Appl. Probab. 2022 | 10.1214/22-aap1799 (Y) | Quartet-based (ASTRAL-one/ASTRAL-multi [M]) | DLCoal-type (coalescent + branching) | **Sample-complexity bounds** in the subcritical and supercritical regimes [A: abstract] | None [M] |
| Legried, "Anomaly zones for uniformly sampled gene trees under GDL" | arXiv 2309.01663 (2023, v3 2024) | arXiv (Y via arXiv API) | Rooted uniformly sampled gene trees | GDL | No rooted anomaly zone for balanced 4-taxon trees. Rooted anomaly zones exist for 4-taxon caterpillars [A: abstract via search] | Not applicable |
| Willson, Roddur, Liu, Zaharias, Warnow, "DISCO" | Syst Biol 2022 (online 2021) | 10.1093/sysbio/syab070 (Y) | ASTRAL-DISCO (also DISCO+ASTRID, CA-DISCO) | GDL | **Proved consistent** (Thm 1) for ASTRAL-DISCO [A]. DISCO+ASTRID and CA-DISCO: empirical only [M] | "provided that ASTRAL-Pro correctly roots and tags each gene family tree". The authors call this requirement "nontrivial" and "far more restrictive" than the model conditions [A] |
| Hakim, Ratul, Bayzid, "wQFM-DISCO" | Bioinformatics Advances 2024 | 10.1093/bioadv/vbae189 (Y) | ASTRAL-DISCO, wQFM-DISCO | **DLCoal** | **Claims a proof**: ASTRAL-DISCO (Thm 3) and wQFM-DISCO (Thm 4) are consistent under DLCoal [A] | "provided that the input gene trees are correctly rooted and tagged", but "correct tagging" under ILS is **never formally defined** [A]. I think the proof is shaky, see the note below [M] |
| Parsons, Liu, Dua, Markin, Molloy, "On (the) correctness of gene tree tagging (and the consistency of ASTRAL-pro) under a unified model of GDL and coalescence" | bioRxiv preprint. v1 2026-01-21, v2 2026-04-12 | 10.64898/2026.01.20.700722 (Y. Crossref and the bioRxiv API show the v1 and v2 titles) | ASTRAL-Pro objective; TREE-QMC with duplication quartets excluded | **DLCoal** | v1 (per search snippets and a UMD talk abstract): consistent for an **"exclusion-only" version** of the A-Pro objective. **v2 drops "consistency of ASTRAL-pro" from its title**, and its abstract only says "study some statistical properties". Consistency of the full A-Pro objective is still open [A: abstracts. Theorems NOT read, because bioRxiv returned HTTP 429] | Correctly rooted and tagged, under a **new definition**: a vertex is a correct duplication if it is the MRCA of ≥1 pair of copies related by a duplication. They also evaluate the accuracy of the A-Pro tagging algorithm in simulation [A: abstract] |
| Molloy & Warnow, "FastMulRFS" | Bioinformatics 2020 (ISMB) | 10.1093/bioinformatics/btaa444 (Y) | FastMulRFS (MulRF objective) | Generic GDL | **Consistent only if "adversarial GDL" is prohibited** (Thm 6). Conj. 7 (low-probability adversarial GDL) is unproven. MulRF has no guarantee (heuristic search) [A] | None (unrooted). "Adversarial GDL" means the gene family tree has a bipartition incompatible with the species tree [A] |
| Wehe, Bansal, Burleigh, Eulenstein, "DupTree" | Bioinformatics 2008 | 10.1093/bioinformatics/btn230 (Y) | GTP, duplication cost | — | No consistency result. The FastMulRFS paper lists it as future work [A]. Willson et al. 2021 say DupTree, iGTP and DynaDup are not proven consistent under GDL [A: snippet] | Rooting is implicit in reconciliation [M] |
| Sapoval & Nakhleh, "On the consistency of duplication, loss, and deep coalescence GTP costs under the MSC" | RECOMB-CG 2026, LNCS (+ bioRxiv 10.64898/2026.02.20.707019) | 10.1007/978-3-032-26891-4_9 (Y) | GTP with any linear combination of dup, loss and DC costs | **MSC only** (ILS, single-copy) | **Proved inconsistent** for all such linear combinations [A: abstract] | Rooted gene trees [M] |
| Than & Rosenberg | JCB 2011 | 10.1089/cmb.2010.0102 (Y) | MDC (minimize deep coalescence) | MSC | **Inconsistent** (anomaly zones for asymmetric 4-taxon trees and all trees with ≥5 taxa) [A: snippet] | Rooted |
| Alanzi & Degnan | PLOS ONE 2021 | 10.1371/journal.pone.0251107 (Y) | Unrooted MDC | MSC | **Inconsistent** [A] | Unrooted |
| Parsons & Bansal, "DupLoss-2" | Syst Biol 2025 | 10.1093/sysbio/syaf073 (Y) | GTP (dup+loss) | GDL simulations | **Empirical only**. Best on most benchmarks [A: abstract] | Rooted gene trees via reconciliation [M] |
| Morel et al., "SpeciesRax" | MBE 2022 | 10.1093/molbev/msab365 (Y) | ML under DTL (UndatedDTL), MiniNJ start | DTL | **Empirical only**. I found no consistency claim [A: search. M] | Rooting comes from the model |
| Emms & Kelly, "STAG" | bioRxiv 2018 | 10.1101/267914 (Y) | STAG (all-species subtrees + distance) | — | **Empirical only** [M] | None |
| Vachaspati & Warnow, "ASTRID" | BMC Genomics 2015 | 10.1186/1471-2164-16-S10-S3 (Y) | ASTRID / NJst | MSC | Consistent under the MSC with complete data. **Not consistent** under i.i.d. taxon deletion (Rhodes, Nute, Warnow, arXiv 2001.07844). No GDL proof for ASTRID-multi. The weighted ASTRID paper notes there are no GDL proofs for distance methods [A: snippets] | None |
| Willson, Roddur, Warnow, "Comparing methods for species tree estimation with GDL" | bioRxiv 2021 | 10.1101/2021.02.05.429947 (Y) | Many methods | DLCoal simulations | Empirical. States that GTP methods are not proven consistent under GDL [A: snippet] | — |
| Yan, Smith, Du, Hahn, Nakhleh | Syst Biol 2022 | 10.1093/sysbio/syab056 (Y) | ILS methods applied to paralogs | DLCoal simulations | Empirical only | None |
| Mishra & Hahn, "Distribution of gene tree topologies with duplication, loss, and coalescence" | bioRxiv 2026 | 10.64898/2026.01.19.700405 (bioRxiv API) | Reconciliation / topology distribution | MSC-DL (unrestricted) | Theory: counts and probabilities of topologies. Methods that ignore coalescence infer spurious duplications and losses [A: abstract] | Relevant to tag errors under ILS |
| Hernandez-Rosales, Hellmuth, Wieseke, Huber, Moulton, Stadler, "From event-labeled gene trees to species trees" | BMC Bioinf 2012 | 10.1186/1471-2105-13-S19-S6 (Y) | Triple-based species tree from event labels | Combinatorial, no stochastic model | Characterization: S is consistent with an event-labeled gene tree iff S displays all species triples rooted at speciation vertices. Polynomial-time [A: summary] | Correct event labels assumed |
| Hellmuth et al., "Orthology relations, symbolic ultrametrics, and cographs" | J Math Biol 2013 | 10.1007/s00285-012-0525-x (Y) | Orthology to event-labeled tree | Combinatorial | Cograph characterization [M] | — |
| Hellmuth et al., "Phylogenomics with paralogs" (ParaPhylo) | PNAS 2015 | 10.1073/pnas.1412770112 (Y) | Cograph editing + max consistent triples + least resolved tree | Combinatorial, with errors | Handles errors by optimization (NP-hard, ILP). No statistical guarantee [A: snippet] | Repairs erroneous orthology |
| Lafond, Dondi, El-Mabrouk | AMB 2016 | 10.1186/s13015-016-0067-7 (Y) | Correcting orthology relations | Combinatorial | NP-hardness and inapproximability of corrections [A: snippet] | Erroneous relations |
| Lafond & Hellmuth, "Reconstruction of time-consistent species trees" | AMB 2020 | 10.1186/s13015-020-00175-0 (Y) | Event-labeled trees with HGT | Combinatorial | Feasibility / algorithms [A: snippet] | Correct labels |

## Answers to the four questions

1. **ASTRAL-Pro definitions and theorem.**
   - [A] Def. 1 (corrected in the corrigendum): node u with children u1, u2 "can be tagged as speciation only if the sets αGu1 and αGu2 are mutually exclusive" (the species sets). All other nodes are tagged duplication.
   - [A] Def. 2: Q is an SQ iff its four leaves come from four distinct species and "the LCA of any three out of four leaves of Q is a speciation node". The anchors are the two degree-3 nodes of G restricted to Q, and they may be duplications.
   - [A] Def. 3: the anchor LCA ψ_G(Q) is the LCA of the two anchors.
   - [A] Def. 4: two SQs on the same four species are equivalent iff they have the same anchor LCA.
   - [A] Def. 5: the per-locus quartet score is the number of equivalence classes whose topology matches S, so each class counts once. MLQST maximizes the sum over loci.
   - [A] Prop. 3: "Under the GDL model, every SQ in every correctly tagged rooted gene tree is isomorphic in topology to the species tree."
   - [A] Thm 2: "Under the GDL model, the solution to the MLQST problem is a statistically consistent estimator of the species tree for correctly rooted and tagged gene trees." The proof is just Prop. 3: every class matches S, so the result is deterministic per locus.
   - [A] Claim 1: "If all nodes on the path between the root r and a node u are tagged as speciations, changing the root to any branch on the path does not alter the PL quartet score," so partially correct rooting suffices.
   - [A] The authors "suspect" consistency holds "even when gene trees are imperfectly rooted and tagged" and say "we leave it to the future to study whether ASTRAL-Pro is statistically consistent under the DLCoal model." So the result is GDL only.
   - [M] ASTRAL-Pro 2 adds no theory. I found no ASTRAL-Pro3 paper, only an ASTER README mention [A].

2. **Rooting/tagging error, LCA tags, DLCoal.**
   - [A] I found **no paper that proves consistency of ASTRAL-Pro or DISCO under rooting or tagging error**, or with the species-overlap (LCA) tags treated as a random estimate. DISCO (2022) says only that a low probability of rooting error "might" still allow a proof.
   - [A] Under DLCoal there are two results:
     - wQFM-DISCO (2024) claims ASTRAL-DISCO and wQFM-DISCO are consistent, assuming correct rooting and tagging. It never defines correct tagging when ILS is present.
     - Parsons et al. (2026 preprint) give the first general definition of correct tagging under DLCoal. Their v1 claims consistency for an "exclusion-only" A-Pro objective only. In v2 the consistency wording is gone from the title and abstract, which suggests the claim was weakened or reframed [inference]. They also measure the accuracy of A-Pro's tagging heuristic in simulation.
   - [A] Markin & Eulenstein prove DLCoal consistency for ASTRAL-one/ASTRAL-multi (no tags), **not** for ASTRAL-Pro. ASTRAL-Pro's text cites them as having shown "that method" (ASTRAL-multi) consistent.
   - [A] **Consistency of the full ASTRAL-Pro objective under DLCoal has not been proved**. A UMD talk abstract (2026) calls it "a major open question".
   - [M/inference, useful for the pilot] Given a correct rooted gene tree under pure GDL, species-overlap tagging never mislabels a true speciation, because speciation children have disjoint species sets. It can only miss a true duplication, after duplication followed by complementary losses. So under GDL the tag errors are one-sided, while under ILS both directions are possible (cf. Mishra & Hahn 2026, spurious duplications).
   - [M] Weak point in the wQFM-DISCO proof: it treats the pruned DISCO subtrees as single-locus MSC trees. Under DLCoal the species-overlap duplication nodes need not correspond to locus births, and that is exactly the gap Parsons et al. address.

3. **Other methods.**
   - ASTRAL-one/ASTRAL-multi are proved consistent under GDL (Legried et al. 2021) and under DLCoal (Markin & Eulenstein 2021), with sample complexity from Hill, Legried, Roch [A].
   - FastMulRFS is consistent only when adversarial GDL is excluded. MulRF has no guarantee [A].
   - ASTRAL-DISCO is consistent under GDL with correct rooting and tagging [A]. DISCO+ASTRID has no proof [M].
   - ASTRID/NJst are consistent under the MSC only, and inconsistent under random missing data [A: snippet]. ASTRID-multi has no GDL proof [A: snippet].
   - STAG, SpeciesRax and DupLoss-2 are empirical only [M / A: abstract].
   - wQFM: DLCoal consistency for wQFM-DISCO is claimed [A]. TREE-QMC: no GDL theorem found, except as the vehicle in Parsons et al. (TREE-QMC-Pro, `--gdl --tagged` [A: snippet]).
   - The Hellmuth/Lafond line is combinatorial: it gives characterizations and NP-hard error-correction problems, with no statistical-consistency results [A: summary].
   - I found no "Allman/Baños/Rhodes identifiability under GDL" paper. The identifiability results under GDL I know of are Legried et al. and Legried's anomaly-zone paper [M].

4. **Inconsistency of gene tree parsimony.**
   - [A: abstract] Under the **MSC**, GTP with any linear combination of duplication, loss and deep-coalescence costs is inconsistent (Sapoval & Nakhleh, RECOMB-CG 2026). Earlier, MDC was shown inconsistent (Than & Rosenberg 2011; unrooted: Alanzi & Degnan 2021).
   - [A: snippet] Under **GDL**, I found **no proof of inconsistency (or consistency)** for duplication-cost GTP (DupTree, iGTP, DupLoss-2). Willson et al. 2021 and Molloy & Warnow 2020 both say it is an open problem. This looks like a genuine gap.

## Could not verify
- **The theorem statements and definitions in Parsons et al. 2026 (both versions).** bioRxiv returned HTTP 429 every time I tried the full text, the PDF, the JATS XML and the supplement. What I report comes from the Crossref and bioRxiv-API abstracts, search-engine snippets and a UMD talk abstract I could not open (404).
- The Legried et al. 2021 and Hill et al. 2022 theorem texts, which I took from other papers' restatements and the abstract.
- STAG and SpeciesRax full texts.
- The exact citation for the "prior fine-grained analysis" of duplication versus DC cost that Sapoval & Nakhleh mention.
