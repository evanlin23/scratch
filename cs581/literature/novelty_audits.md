# Novelty audits of the top candidates (2026-10-10)

Deeper prior-art checks (2015-2026, cited-by lists of the key papers) for the ideas ranked highest in
`ranking_v1.md`, because the slides' open questions may be stale. Verdicts: NOVEL / PARTIALLY KNOWN /
ALREADY ESTABLISHED.

## Fragment-aware two-phase ML (branch claude/cs581-ml) - PARTIALLY KNOWN

Method (frag.py): backbone = sequences >= 50% of median length; FastTree on the backbone; RAxML-NG
over all sequences with the backbone tree as a non-comprehensive `--tree-constraint`; optional
unconstrained polish.

Closest prior work:
- Smirnov & Warnow 2021, Syst Biol 70:268, doi:10.1093/sysbio/syaa058: full-length tree + UPP/SEPP +
  pplacer/APPLES placement on ROSE 1000M1-M4 and RNASim; no constrained ML search or post-placement ML.
- Berger, Krompass & Stamatakis 2011, Syst Biol 60:291, doi:10.1093/sysbio/syr010: names constrained
  search with the reference tree as backbone as an existing approach (then proposes EPA).
- Balaban et al. 2024 (uDance), Nat Biotechnol 42:768, doi:10.1038/s41587-023-01868-8: backbone ->
  place -> local re-estimation; filters fragmentary sequences.
- Springer et al. 2001 (molecular scaffold), PNAS 98:6241; Catalano et al. 2025, Syst Biol 74:672,
  doi:10.1093/sysbio/syaf025: backbone-constrained placement of fragmentary taxa is routine.
- Sayyari, Whitfield & Mirarab 2017, MBE 34:3279, doi:10.1093/molbev/msx261: fragments hurt gene trees;
  filtering proposed, no constrained re-insertion.
- Park, Zaharias & Warnow 2021, Algorithms 14:148, doi:10.3390/a14050148: 1000M1-HF baselines.

Still novel: a controlled benchmark of backbone-constrained joint ML insertion (with budgeted polish) vs
unconstrained RAxML-NG / IQ-TREE / FastTree on 1000M1-HF and the Smirnov-Warnow conditions (FN + CPU).
Caveats: pilot used the TRUE alignment, n = 5, one condition, +0.5 FN gap is 2/0/3 (n.s.). Do not call
the method new.

Safe phrasing: "We evaluate a known but unbenchmarked strategy for alignments with many fragments ...
(cf. Berger et al. 2011; uDance) ... measuring the accuracy/CPU trade-off against unconstrained ML."

Search gaps: Park 2021 full text (403); Smirnov-Warnow supplement; theses; PUmPER (Izquierdo-Carrasco
et al. 2014); IQ-TREE 3 tree-updating features.

## GTM blending: ML-scored, constraint-preserving SPR from GTM (branch claude/cs581-gtm) - PARTIALLY KNOWN (novel core)

The open problem is NOT stale: Warnow's 2023, ICERM 2024 and IMSI 2025 talks still list "Develop a better
DTM approach that allows blending" (tandy.cs.illinois.edu/Warnow-DTMs-IMSI-2025.pdf).

Closest prior work:
- Smirnov & Warnow 2020 (GTM), BMC Genomics 21(S2):235, doi:10.1186/s12864-020-6605-1: blending is NP-hard,
  unblended is polynomial; calls for blended DTMs.
- Molloy & Warnow 2019 (NJMerge, AMB 14:14, doi:10.1186/s13015-019-0151-x; TreeMerge, Bioinformatics
  35:i417, doi:10.1093/bioinformatics/btz344): blending mergers, distance-based / greedy.
- Zhang, Rao & Warnow 2019 (Constrained-INC, AMB 14:2, doi:10.1186/s13015-019-0136-9); Le et al. 2021
  (INC-ML, TCBB, doi:10.1109/TCBB.2020.2990867): full blending by constrained insertion, distance/quartet voting.
- Park, Zaharias & Warnow 2021, Algorithms 14:148: benchmarks GTM/TreeMerge/C-INC; no constrained search after merging.
- Rec-I-DCM3 (Roshan et al. 2004): merge then unconstrained refinement (overlapping subsets).
- Chernomor, Minh & von Haeseler 2015, J Comput Biol 22:1129, doi:10.1089/cmb.2015.0146; Gentrius (MBE 2024,
  msae219): conditions for SPR/NNI changing induced subtrees (our feasibility test is essentially known).
- MrBayes / TNT enforce several partial constraints; IQ-TREE / RAxML-NG take one constraint tree.

Still novel: likelihood on the full alignment to decide blending; SPR search keeping all k subset trees induced
as a GTM post-processor; evidence it repairs GTM when subsets are not clades; headroom analysis on published data.

Safe phrasing: "To our knowledge no published DTM uses sequence likelihood to decide blending. We test a simple
constrained ML local search started from GTM that accepts only SPR moves keeping every subset tree induced
(feasibility test adapted from Chernomor et al. 2015). In 200-taxon simulations where subsets are not clades it
reduces GTM's error; on published benchmarks the available gain is small and we observed none significant."
Avoid: "first blending DTM", "first search under multiple partial constraints", "solves the open problem",
unqualified "6-14 points" (two of three conditions used exact subset trees).

Search gaps: Google Scholar cited-by (used Semantic Scholar: 10 papers cite GTM, OpenAlex says 12); Park 2021
and Smirnov 2021 thesis full texts (403); RECOMB-CG/WABI/ISMB 2025-26 not searched systematically.

## Soft-constraint MAGUS (main track) - PARTIALLY KNOWN

MAGUS's own code (magus-msa 0.2.0; commit 2021-01-08 "Added -c flag") has `-c false`: every sequence becomes
its own group and the subset alignments are appended as backbones, i.e. our m = subset size. README: "drastically"
slower, "strongly not recommended above 200 sequences". No paper evaluates it. `--graphbuildmethod initial` is a
precedent for self-derived evidence (and has the skeleton bug we found).

Closest prior work:
- Smirnov & Warnow 2021, MAGUS, Bioinformatics 37:1666, doi:10.1093/bioinformatics/btaa992: subsets are absolute constraints.
- Zaharias, Smirnov & Warnow 2022, MWT-AM, TCBB, doi:10.1109/TCBB.2022.3191848: GCM variants (FM, MWTgreedy/search,
  optimizer), all with hard constraints (our fm+opt row).
- Smirnov 2021, Recursive MAGUS, PLoS Comput Biol 17:e1008950, doi:10.1371/journal.pcbi.1008950: hard constraints.
- Opal (Wheeler & Kececioglu 2007, doi:10.1093/bioinformatics/btm226) and MUSCLE (Edgar 2004, doi:10.1186/1471-2105-5-113):
  split-and-realign refinement (conceptual precedent).
- M-Coffee (Wallace et al. 2006, doi:10.1093/nar/gkl091): alignments as soft consistency evidence.
- Ian Chen, CS581 Sp2025 "Exploring the Effect of Iteration on MAGUS" (tandy.cs.illinois.edu/ian-may6.pdf): PASTA-style
  re-decomposition; no constraint relaxation, no self-evidence.
- No overlap (fixed constraints/backbones): EMMA, WITCH, WITCH-NG, HMMerge, MAGUS+eHMMs, UPP2, TWILIGHT 2025, FAMSA2, MuSAlS.

Still novel: partial relaxation via similarity groups (tractable, ~9% extra runtime); self-derived evidence (MAGUS
output -> 8 seqs/subset -> HMM-extended -> soft re-merge); first controlled held-out evaluation of relaxing GCM
constraints (random-split control and oracle decomposition support the mechanism).

Safe phrasing: "MAGUS's software includes an unconstrained mode (-c false) described as impractical beyond ~200
sequences; to our knowledge no published study evaluates it or any partial relaxation. We propose partial relaxation
by splitting each subset alignment into similarity groups, combined with evidence derived from MAGUS's own output,
connecting GCM to classical split-and-realign refinement (MUSCLE, Opal)." Avoid "first to relax GCM constraints".
TODO if chosen: run `-c false` as a baseline on the 1000-sequence conditions.

Search gaps: Google Scholar cited-by (Semantic Scholar: 75 cite MAGUS, 13 Recursive MAGUS, 4 MWT-AM); MAGUS_CP repo (503);
WABI/RECOMB-CG 2026 abstracts; other universities' course projects.

## ASTRAL-Pro inconsistency under GDL via its own rooting/tagging (branch claude/cs581-gdlcons) - PARTIALLY KNOWN (core negative result novel)

Known: ASTRAL-Pro is consistent under GDL given true roots and tags; its tagger can mislabel under adversarial
duplication-then-loss; consistency with imperfect rooting/tagging is an open conjecture that the authors expect to hold.
Not found anywhere: an inconsistency result caused by its own rooting/tagging, an exact limiting-score formula, or a
random-tag-error threshold.

Closest prior work:
- Zhang, Scornavacca, Molloy & Mirarab 2020, ASTRAL-Pro, MBE 37:3292, doi:10.1093/molbev/msaa139: Theorem 2 (consistent
  for correctly rooted/tagged trees); tagger "not guaranteed to find the correct tags or the root"; authors "suspect"
  consistency holds with imperfect rooting/tagging.
- Parsons, Liu, Dua, Markin & Molloy 2026, bioRxiv doi:10.64898/2026.01.20.700722: v1 Theorem 1 (exclusion-only objective
  consistent under DLCoal given correct tags) withdrawn in v2 (Apr 2026) as Conjecture 1, "an open question"; mention that
  adversarial GDL misleads A-pro's tagging; no inconsistency result.
- Willson et al. 2022, DISCO, Syst Biol 71:610, doi:10.1093/sysbio/syab070: consistency "provided that ASTRAL-Pro correctly
  roots and tags", "not likely to hold on many conditions"; speculate random low-probability tag error may still allow a
  proof (= the slide question, left open).
- Molloy & Warnow 2020, FastMulRFS, Bioinformatics 36:i57, doi:10.1093/bioinformatics/btaa444: adversarial GDL; GDL models
  with fixed (constant) rates across edges.
- Morel et al. 2022, SpeciesRax, MBE, doi:10.1093/molbev/msab365: ASTRAL-Pro "may mislabel paralogs as orthologs under
  high loss rates" (simulation).
- Background: Smith & Hahn 2022 (doi:10.1093/sysbio/syab097); Xiong et al. 2022 (doi:10.1093/sysbio/syac040); Legried et al.
  2021 (doi:10.1089/cmb.2020.0424); Markin & Eulenstein 2021 (doi:10.1093/bioinformatics/btab414); Hill, Legried & Roch 2022
  (doi:10.1214/22-aap1799).

Novel: explicit counterexamples (stock ASTRAL-Pro inconsistent on true gene trees via its own species-overlap tagging or
min-duplication rooting); exact 4-taxon condition O + H_AB > max(H_AC, H_BC); threshold q* for random hidden-paralog
mislabelling.
Caveats: all counterexamples use branch-specific, extreme rates (standard GDL models use constant rates) -> MUST check
homogeneous rates; predicted q* = 0.079 outside observed bracket 0.10-0.15 (approximate); numerical, not proofs.

Safe phrasing: "Zhang et al. (2020) proved ASTRAL-Pro consistent under GDL given correctly rooted and tagged gene trees and
conjectured this extends to imperfect rooting and tagging. We give evidence that, in general, it does not: under a GDL model
with branch-specific rates we derive exact limiting 4-taxon ASTRAL-Pro scores and exhibit configurations where its own
tagging or min-duplication rooting makes it converge to a wrong species tree from true gene-family trees. For random
hidden-paralog tag errors we characterize a threshold below which consistency is retained. Whether inconsistency occurs under
the constant-rate model remains open."

### Is ASTRAL-Pro state of the art? Default tool and standard baseline, not uniformly best.
- DISCO 2022: ASTRAL-Pro, ASTRID-DISCO, SpeciesRax top group; ASTRID-DISCO/CA-DISCO best with few genes/short seqs/missing data.
- SpeciesRax 2022: best in own simulations; no method dominates empirically.
- wQFM-DISCO 2024 (Bioinf Adv, doi:10.1093/bioadv/vbae189): matches or beats ASTRAL-Pro.
- AleRax 2024 (Bioinformatics, doi:10.1093/bioinformatics/btae162): more robust than SpeciesRax / ASTRAL-Pro 2; slower.
- Weiner et al. 2025 (Peer Community Journal, doi:10.24072/pcjournal.579): microbial + HGT: AleRax best, ASTRAL-Pro 2 worst.
- DupLoss-2 (Syst Biol 2025/26, doi:10.1093/sysbio/syaf073): ~10% less error than best existing method on most benchmarks.
- wQFM-GDL (RECOMB-CG 2026, doi:10.1007/978-3-032-26891-4_8): beats ASTRAL-Pro3, SpeciesRax, FastMulRFS, DupLoss-2 in 134/156
  conditions; ~25% less error than ASTRAL-Pro3 at 200-500 taxa.
- STAG absent from 2022-2026 benchmarks; FastMulRFS consistently outperformed. (Most benchmarks are by the methods' own authors.)

Search gaps: Legried 2021 / Hill 2022 theorem texts; RECOMB-CG 2026 TOC (dblp blocked); WABI/ISMB 2026; theses.

## ASTRID/NJst sample complexity vs ASTRAL (branch claude/cs581-samplecx) - theory OPEN; empirical PARTIALLY KNOWN

No 2019-2026 paper proves an upper or lower sample-complexity bound in n for ASTRID/NJst. The qualitative n-effect was
reported in 2015; a controlled k95-vs-n measurement (caterpillar vs balanced) appears unpublished.

Corrections to our claim:
- Not a formal conjecture: Roch (arXiv:1812.08357, sec. 3) says his bound "does not in fact lead to a bound on the sample
  complexity" and that correlations "could drastically lower" it; he only says m must be >= linear in n "to make all variances
  negligible". Do not call it "Roch's conjecture".
- The pilot's trees are not Roch's worst case (caterpillar with short f branches alternating with long 4 log n branches).
- n = 8 -> 64 cannot yet separate log n (x2) from n (x8): ASTRAL's constant grew ~1.9x, ASTRID's ~3.1x.
- "Nobody measured beyond n = 8" is too strong (Mirarab & Warnow 2015).

Closest prior work:
- Roch 2018, RECOMB-CG, LNCS 11183:196, doi:10.1007/978-3-030-00834-5_11: Var >= C d/m, max Var >= C n/m (caterpillar
  construction); suggests NJ on all distances needs m >= linear in n vs log n for ASTRAL; no simulations.
- Shekhar, Roch & Mirarab 2018, TCBB 15:1738, doi:10.1109/TCBB.2017.2757930: Theta(f^-2 log n) for ASTRAL*; simulations at
  n = 8 only; NJst also ~f^-2; never varies n.
- Mirarab & Warnow 2015, ASTRAL-II, Bioinformatics 31:i44, doi:10.1093/bioinformatics/btv234: NJst vs ASTRAL-II for n = 10-1000
  at 50/200/1000 genes; gap depends on n (P = 0.0004), not on genes.
- Vachaspati & Warnow 2015, ASTRID, doi:10.1186/1471-2164-16-S10-S3; Liu & Warnow 2023, wASTRID, doi:10.1186/s13015-023-00230-6:
  "dataset dependent", no n sweep.
- Allman, Degnan & Rhodes 2018, TCBB 15:337, doi:10.1109/TCBB.2016.2604812: consistency only.
- Background: Lacey & Chang 2006 (doi:10.1016/j.mbs.2005.11.003), NJ needs sequence length growing with tree diameter;
  Dasarathy et al. 2022 (doi:10.1007/s00285-022-01731-5, 3 taxa); Hill & Roch 2025 (doi:10.1007/s11538-025-01533-y).

Novel: k95(n, f) for ASTRID/NJst vs ASTRAL under controlled MSC, caterpillar vs balanced, to n = 128; testing Roch's own
construction; log vs polylog vs linear growth.
Safe phrasing: "Roch (2018) proved a variance lower bound suggesting internode-distance methods may need a number of genes
growing at least linearly with n, versus logarithmically for ASTRAL, and left the sample complexity open. Shekhar et al. (2018)
measured genes needed only at n = 8, and Mirarab & Warnow (2015) observed an NJst-ASTRAL gap growing with n at fixed gene
counts. To our knowledge the growth of the genes-needed threshold with n has not been measured; we measure it on caterpillar,
balanced and Roch's worst-case trees."
Add: Roch's alternating short/long caterpillar; a fast-converging (local-distance) ASTRID variant, as Roch proposes.

Search gaps: Google Scholar cited-by (Semantic Scholar/OpenCitations list one citer of Roch 2018); NJMerge, TREE-QMC and FASTRAL
experiments not read in full; theses by title only.

## ASTRID-Pro: GDL-corrected internode distance (branch claude/cs581-gdl) - NOVEL (narrowly)

No 2018-Oct 2026 paper proves consistency or inconsistency of ASTRID, NJst, STAR, STAG, MiniNJ or any averaged
internode-distance method under GDL or DLCoal. DISCO (2022): "to date, no distance-matrix method has been proven
statistically consistent under a model of GDL"; weighted ASTRID (2023) notes "the lack of such proofs". Parsons et al.
2026 is about ASTRAL-pro only.

Closest prior work:
- Willson et al. 2022, DISCO, Syst Biol 71:610, doi:10.1093/sysbio/syab070: ASTRAL-DISCO consistent given correct
  rooting/tagging; consistency of ASTRID-DISCO, ASTRID-multi and MiniNJ open.
- Zhang et al. 2020, ASTRAL-Pro, MBE 37:3292, doi:10.1093/molbev/msaa139: consistent via tagging + speciation-driven quartets
  (ASTRID-Pro is essentially its distance analogue).
- Rhodes, Nute & Warnow 2020, arXiv:2001.07844: NJst/ASTRID inconsistent under MSC + i.i.d. taxon deletion (additive matrix,
  wrong topology); our lemma is a missing-data internode expectation; GDL twist: supercritical branch lets s_yx > s_xA + s_xB.
- Morel, Williams & Stamatakis 2023, Asteroid, Bioinformatics 39:btac832, doi:10.1093/bioinformatics/btac832: proven fix for
  missingness bias under any gene-independent deletion model (multi-copy via MiniNJ); competes with our survival-weighted
  correction; not claimed for GDL.
- Legried et al. 2021, JCB 28:452, doi:10.1089/cmb.2020.0424 (ASTRAL-one/multi consistent under GDL; ASTRID-multi empirical
  only); Legried 2023, arXiv:2309.01663 (GDL anomaly zones for rooted caterpillar quartets at large lambda; none for balanced,
  mirroring our pattern).
- Hernandez-Rosales et al. 2012, BMC Bioinf 13(S19):S6, doi:10.1186/1471-2105-13-S19-S6: with correct event labels, no ILS,
  the species tree is identifiable (contribution is about the estimator, not identifiability).

Novel: additivity lemma and four-point analysis for a speciation-only, ortholog-only internode distance under GDL; sufficient
condition lambda_e <= mu_e on every branch; 4-taxon supercritical counterexample; survival-reweighted correction; held-out
pre-registered gain over ASTRID-multi. ASTRID-multi itself stays unresolved (do not claim to answer that slide question).

Safe phrasing: "To our knowledge, we give the first consistency analysis of an ASTRID/NJst-style averaged internode distance
under GDL: assuming no ILS and correct rooting and tagging, an orthology-restricted, speciation-node distance converges to a
tree metric whose topology is the species tree whenever no branch is supercritical, and can be positively misleading otherwise.
We compare a survival-reweighted correction with Asteroid-style induced-length correction. Consistent quartet-based pipelines
(ASTRAL-Pro, ASTRAL-DISCO) already exist in this setting."

Search gaps: Google Scholar; citers of Rhodes et al. arXiv; theses (Willson, Legried, Hill, B. Liu); RECOMB/WABI/ISMB 2026;
STAG full text; (Parsons v2 full text read 2026-10-10 16:25 UTC: DLCoal + correct tagging, consistency left as Conjecture 1; no inconsistency result; no overlap); whether Asteroid's theorem covers GDL-induced (dependent) deletion.

## Summary

| candidate | verdict | novel core | must avoid / must add |
|---|---|---|---|
| ASTRID-Pro GDL distance theorem | NOVEL (narrowly) | first consistency analysis of an internode-distance method under GDL; supercritical counterexample; correction | do not claim ASTRID-multi answered; compare with Asteroid |
| ASTRAL-Pro inconsistency via own tagging | PARTIALLY KNOWN, core novel | counterexamples to Zhang et al.'s conjecture; exact 4-taxon formula; threshold | must test constant-rate GDL; numerical, not proofs |
| GTM blending (ML-scored constrained SPR) | PARTIALLY KNOWN, core novel | likelihood-decided blending; open problem still listed 2025 | not "first blending DTM"; gains only in simulation |
| Soft-constraint MAGUS | PARTIALLY KNOWN | partial relaxation + self-derived evidence; first held-out test | must run MAGUS `-c false` baseline; not "first to relax constraints" |
| ASTRID/NJst sample complexity | theory OPEN; empirical PARTIALLY KNOWN | k95 vs n, caterpillar vs balanced, Roch's construction | not "Roch's conjecture"; use Roch's worst-case tree |
| Fragment-aware ML | PARTIALLY KNOWN (method is known practice) | a benchmark only | used true alignment, n = 5; do not call the method new |

## Consistency-filtered GCM evidence (branches claude/cs581-bbevidence, claude/cs581-protcons) - PARTIALLY KNOWN (novel application)

Audit 2026-10-10 (literature only).
- **(a) Self-consistency mask on MAGUS's backbones:** the score is old (GUIDANCE residue-pair score,
  doi:10.1093/molbev/msq066; trimAl consistency, doi:10.1093/bioinformatics/btp348; TCS,
  doi:10.1093/molbev/msu117, which keeps residues scoring > 0.6 for 80.4% specificity / 73.9% sensitivity on
  BAliBASE/PREFAB). Applying it to prune GCM's alignment-graph evidence: not found.
- **(b) Two-aligner intersection as GCM evidence:** the multi-aligner idea is M-Coffee's (doi:10.1093/nar/gkl091;
  also Heads-or-Tails, MergeAlign doi:10.1186/1471-2105-13-117, MSARC doi:10.1186/1748-7188-9-12). Inside GCM:
  not found. Zaharias, Smirnov & Warnow (TCBB 2022, doi:10.1109/TCBB.2022.3191848) list as future work "explore
  modifying how GCM defines the weights on the pairs of columns" and other backbone aligners (PROMALS, ProbCons).
  MAGUS issue #34 (Apr 2025, C. Shen, "we would like to extend MAGUS") may indicate concurrent lab work.
- **(c) "GCM is precision-limited on proteins":** not stated anywhere and contradicted as a data-type rule (MAGUS on
  HomFam SPFN 0.222 vs SPFP 0.088, Recursive MAGUS supplement; 16S.3 SP 92.1 vs Modeler 86.5, TCBB Table 3). Only
  safe as "on BAliBASE RV100, GCM's errors were dominated by SPFP".
- **MAGUS itself:** edge weight = number of supporting residue pairs; backbones "not constrained to be consistent
  ... with each other"; MCL is described as extending "the consistency principle ... to longer paths", so a
  reviewer will ask why an explicit filter is not redundant with MCL. `-b` already accepts external backbones.
- **Safe:** "To our knowledge, agreement-based filtering has not been applied to GCM's alignment-graph evidence";
  "we bring GUIDANCE/TCS-style residue-pair reliability into MAGUS's graph-building step".
- **Avoid:** "first consistency-based ...", "new consistency score", a general protein-vs-nucleotide mechanism,
  implying filtering helps trees (Tan et al. 2015, doi:10.1093/sysbio/syv033).
- **Baselines a reviewer will demand:** tuned MAGUS (more/larger backbones, MCL inflation, fm+opt trace); a
  global edge-support threshold (>= k of 10 backbones); TCS-filtered/weighted backbones; GUIDANCE2-style perturbed
  backbones; mixed-aligner union (M-Coffee style); M-Coffee/T-Coffee/MergeAlign/MAFFT L-INS-i on the full data
  (feasible at BAliBASE sizes); SPFN and SPFP separately, with HomFam and CRW; the gate vs MUMSA overlap score,
  Muscle5 confidence and TCS as reference-free difficulty predictors.
- **Not accessed:** MAGUS supplement, Smirnov thesis, Illinois Data Bank pages (403), full texts of UPP2/EMMA/WITCH-NG.

## Quick checks, 2026-10-10 19:40 UTC (check-in 8)

**ASTRID-Pro (orthology-restricted speciation-node distance + FastME).**
- No method of that name was found.
- Closest prior work, to cite:
  - SpeciesRax's MiniNJ (Morel et al., MBE 2022): a distance-based starting tree for multi-copy gene families. It is
    fast and has no consistency analysis.
  - CASTLES-Pro (GBE 2025, doi:10.1093/gbe/evaf200): estimates species-tree *branch lengths* under GDL+ILS from
    orthologous quartets (leaf pairs meeting at speciation nodes). Related use of tags; not a topology method and no
    consistency theorem for a distance.
  - ASTRID-multi, ASTRID-DISCO, Asteroid.
- Search found no consistency proof for any GDL distance method, so the theorem's novelty holds.
- Safe claim: "first consistency theorem (with an exact condition) for a GDL distance method, given correct tags."

**GCM edge-support threshold (delete cross-subset edges supported by < k of the 10 backbones).**
- Not described in MAGUS (Smirnov & Warnow 2021) or GCM-improvements (Zaharias, Smirnov & Warnow, TCBB 2022). The
  latter lists changing GCM's edge weights as future work.
- Still to check: MAGUS's code (any `minWeight`-style option), and T-Coffee/M-Coffee library pruning, which is the
  likely closest analogue (low-weight library pairs).
