# Forest+DTM pilot: do DMR forest components + GTM beat NJ / FastME at short sequence lengths?

CS581 (Fall 2026) project pilot, candidate idea B from `cs581/literature/wide_scan.md` §2.1.
Everything here is reproducible from `cs581/forest/code/` (see "Reproducing" at the end).

## 1. The question and where it comes from

* Divide-and-conquer deck (`581-Divide-and-Conquer-trees-2023`), summary slide: *"The Guide Tree
  Merger (GTM) is the current leading DTM technique, based on empirical performance … However, GTM
  does NOT allow blending … Open problem: Develop a better DTM approach."*
* Lectures 1–2 (statistical consistency, sequence-length requirements, absolute fast converging
  methods): NJ can need exponentially long sequences; AFC methods (short quartets, DCM-NJ) only
  trust short distances.
* Kim, Lokhov, Vuffray, Romero-Severson & Goldberg, *BMC Bioinformatics* 27:186 (2026),
  doi:10.1186/s12859-026-06488-y, implemented the Daskalakis–Mossel–Roch forest algorithm
  (*SIAM J. Discrete Math.* 2011, doi:10.1137/09075576X) in `lanl/distphylo`. Forest returns only the
  reliably reconstructable part of the tree as leaf-disjoint components; it beat NJ only in limited
  cases. Their Discussion proposes, as future work, feeding forest components to a DTM.

**Pilot question.** If we take the forest components as constraint trees and merge them with GTM
under an NJ or FastME guide tree (a modern DCM-NJ), do we get a *full* tree that is more accurate
than NJ / FastME at short sequence lengths, in particular in hard (long-branch, deep) regimes?

## 2. Prior art (≈20 min search; details and DOIs in `lit/prior_art.md`)

* **Not done.** Forest+DTM appears only as a suggestion in Kim et al. 2026 (0 citing works found on
  Semantic Scholar / OpenAlex). None of ~50 works citing DMR 2011 merges forest components with a
  DTM or supertree method.
* Closest ancestors: **DCM-NJ / DCM1** (Huson, Nettles & Warnow, *JCB* 1999,
  doi:10.1089/106652799318337): threshold graph → overlapping small-diameter subsets → NJ per
  subset → strict-consensus merge; **DCM2** (Huson, Vawter & Warnow, ISMB 1999; no DOI found);
  **DCM-NJ+MP** (Nakhleh et al., *Bioinformatics* 2001, doi:10.1093/bioinformatics/17.suppl_1.S190);
  **Short Quartet Method** (Erdős, Steel, Székely & Warnow 1999,
  doi:10.1002/(SICI)1098-2418(199903)14:2<153::AID-RSA3>3.0.CO;2-R and 10.1016/S0304-3975(99)00028-6).
  These use *overlapping* subsets and consensus/dyadic-closure merges, not a disjoint forest + DTM.
* Other "reliable partial tree" methods, none of which merges components: Daskalakis et al. RECOMB
  2006 (doi:10.1007/11732990_24), Mossel TCBB 2007 (doi:10.1109/TCBB.2007.1010), Gronau, Moran &
  Snir RSA 2012 (doi:10.1002/rsa.20372), Mihaescu, Hill & Rao Algorithmica 2013
  (doi:10.1007/s00453-012-9644-4), Brown & Truszkowski AMB 2012 (doi:10.1186/1748-7188-7-32).
* DTMs: NJMerge (doi:10.1186/s13015-019-0151-x), TreeMerge (doi:10.1093/bioinformatics/btz344),
  Constrained-INC (doi:10.1186/s13015-019-0136-9), GTM (Smirnov & Warnow, *BMC Genomics* 2020,
  doi:10.1186/s12864-020-6605-1), Park, Zaharias & Warnow, *Algorithms* 2021 (doi:10.3390/a14050148).

So novelty is not the risk; usefulness is.

REPRO_PLACEHOLDER
