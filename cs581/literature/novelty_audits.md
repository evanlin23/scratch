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
