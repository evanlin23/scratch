---
title: "CS581 Project Proposal: Fast, Provably Consistent Species Trees Under Gene Duplication and Loss (and Where ASTRAL-Pro Fails)"
author: "Evan Lin · [NetID] · October 2026 (draft, pilot results provisional)"
---

**Problem.** Gene families evolve by gene duplication and loss (GDL) as well as speciation, so gene trees have
many copies per species and paralogs. ASTRAL-Pro (Zhang et al., MBE 2020) is the standard quartet method for such
data. It is proved statistically consistent under GDL *if the gene trees are correctly rooted and every node is
correctly tagged as a duplication or a speciation*, but it roots and tags the trees itself with a heuristic.
Distance methods (ASTRID, NJst) are much faster, but no consistency result exists for them under GDL. The
phylogenomics lecture (part 2) lists both questions as open: *"Is ASTRAL-Pro statistically consistent for GDL
under some random model of error for rooting and tagging?"* and *"Can we find a distance correction for GDL so
that ASTRID and NJst are statistically consistent under GDL models?"* (also: "Unknown if distance-based species
tree estimation (e.g., ASTRID-multi) is statistically consistent under GDL models"). The most recent work (Parsons,
Liu, Dua, Markin & Molloy, bioRxiv 2026, v2) studies ASTRAL-Pro with *correct* tags under a model that adds deep
coalescence, leaves its consistency as a conjecture, and defers specific rate settings to future work.

**Pilot findings (done; proofs written, numbers provisional).**
*(1) ASTRAL-Pro with its own tagging is not consistent.* On the 4-taxon tree (((A,B),C),D) under GDL with
branch-specific rates and true gene trees, we derived the exact limiting ASTRAL-Pro score of each topology with the
correct root and ASTRAL-Pro's species-overlap tagging. Proved:
- *Sign law.* At a duplication with unequal numbers of copies on its two sides, the wrong pairing gains iff C's
  copies survive better than B's. So ASTRAL-Pro is consistent whenever C survives no better than A and B, whatever
  the duplication process.
- *Counterexample.* A rigorous instance where every wrong topology outscores the true one (2.4025e-4 vs 3.2663e-4
  per copy; closed form checked with 50-digit interval arithmetic). In the same tree, with hidden paralogs mislabelled at random, it
  fails above a threshold of q* = 0.188.
- *Re-tagging cannot fix it.* Re-tagging by reconciliation against any tree T makes ASTRAL-Pro return T, so tagging
  against a first-pass tree cannot correct it.

In simulation, ASTRAL-Pro's own rooting makes it 2.6× worse (z = −117). Measured: ASTRAL-Pro3, DISCO+ASTRAL and
wQFM-GDL return the wrong tree in 10/10 blocks of 1,000 families, also with FastTree-estimated gene trees (5/5),
while duplication-loss parsimony is correct (0/10 wrong). The failure needs high rates with ~3× rate differences
between branches: none at published simulation rates.

*(2) A distance method with an exact consistency condition.* "ASTRID-Pro" averages, over orthologous gene copies,
the number of speciation nodes on the path between them. With true trees and tags, we proved:
- *Limit.* It converges to a metric that is additive on the species tree, with explicit edge lengths
  β(c) = ½[s(c₁) + s(c₂) + s(sib c) − s(c)], where s is a copy's survival probability.
- *Exact condition.* It is consistent iff every β(c) > 0, which holds whenever no branch is supercritical
  (duplication rate > loss rate). A counterexample with β = −0.267 matches simulation to 0.003.
- *Correction.* Survival-reweighting is consistent under any rates, given a correct first-pass tree.

Every published benchmark condition has β ≥ 0.46. *Empirically, it is the fast method:* a streaming C++
implementation (MinDup rooting, species-overlap tags, FastME) on the FastMulRFS (= ASTRAL-Pro S100) and DISCO
simulations (100 taxa, 12+ conditions, estimated gene trees), paired by replicate:
- *vs ASTRAL-Pro3:* lower FN rate (−0.011, 46/14/11, p = 2e-6). Our ASTRAL-Pro3 reproduces the published
  ASTRAL-Pro trees on the same inputs.
- *Speed:* 0.3 s median vs 12.9 s; on 1000-gene DISCO inputs 1–40 s vs 6–29 min, where ASTRAL-Pro3 timed out on the
  three heaviest conditions; 14 s at 1,000 taxa.
- *vs other methods:* it beats ASTRID-multi (−0.004, p = 0.019), DISCO+ASTRAL and FastMulRFS (−0.02). It ties
  wQFM-GDL, ASTRID-DISCO and Asteroid, which have no GDL consistency proof.
- *Empirical data:* on 1KP it matches ASTRAL-Pro3 (5 FN vs 10 for ASTRID-multi).

**Research questions.**
(1) Can the stock-rooting result (the simulation certificate) be proved, and does the 4-taxon failure lift to
n taxa, constant-rate models, or the DLCoal model with deep coalescence?
(2) For ASTRID-Pro, what happens with *estimated* tags and roots (hidden paralogy, MinDup rooting), and does the
reweighting correction stay consistent when the first-pass tree is estimated?
(3) Which practical methods inherit each failure (ASTRAL-Pro3, DISCO+ASTRAL, wQFM-GDL, SpeciesRax, FastMulRFS,
ASTRID-multi), and does any published or empirical dataset sit near the failure region?
(4) Do the two answers share a cause? Both failures are driven by unequal survival across lineages; we will state
one survival condition under which both methods are consistent, if it exists.

**Approach and evaluation.** Theory first: complete the proofs (4-taxon closed forms, sign law, the additive-limit
theorem), then try the extensions in RQ1/RQ2 with the same machinery (copy-number generating functions for
linear birth-death). Every theorem gets a numerical check: exact formulas vs simulation with SimPhy-style GDL
simulators (true gene trees, then FastTree gene trees). Methods: ASTRAL-Pro3, DISCO+ASTRAL, wQFM-GDL,
ASTRID-multi, ASTRID-DISCO, FastMulRFS, our ASTRID-Pro. Data: the DISCO and FastMulRFS simulation conditions
(published rates) and our failure-region conditions. Metrics: species-tree FN rate, fraction of replicates
returning the wrong tree as the number of gene families grows, and runtime.

**Timeline.** Week 1: write the proved results cleanly; attempt the stock-rooting proof. Week 2: n-taxon and
constant-rate extensions; ASTRID-Pro with estimated tags. Week 3: method comparison on published and
failure-region conditions; scan published rates for near-failure branches. Week 4: write-up.

**Risks.** (a) *Practical reach.* Both failures need extreme or heterogeneous rates, so the theory changes no
answer on standard benchmarks. The method contribution is speed at ASTRAL-Pro3-level accuracy; it does not beat
the fast heuristics ASTRID-DISCO or Asteroid. Most of the gain over ASTRAL-Pro3 comes from the S100 data (on DISCO
only 14 pairs, n.s.), so RQ3 must test it more widely. (b) The stock-rooting proof and the DLCoal extension may not close in 4 weeks; the fallback is the proved
correct-root theorems plus rigorous numerics. (c) Novelty is checked against Zhang et al. 2020, Willson et al.
2022 (DISCO), Molloy & Warnow 2020, Legried et al. 2021, Markin & Eulenstein 2021 and Parsons et al. 2026: none
gives an inconsistency result for ASTRAL-Pro's own tagging or a consistency theorem for a GDL distance method; the
reconciliation observation may be folklore and will be presented as such.

**Course connection.** Both questions are printed verbatim in the phylogenomics lecture (part 2); ASTRAL,
ASTRID/NJst, statistical consistency and GDL are core course material, and the project uses only the course's
analysis tools (quartet probabilities, additive metrics and the four-point condition).
