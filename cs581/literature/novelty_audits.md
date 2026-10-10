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
