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
