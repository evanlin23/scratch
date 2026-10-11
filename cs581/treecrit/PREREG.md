AI-assisted (Claude), exploration code for CS581 project

# Pre-registration: posterior-guided column splitting for trees (treecrit Part B)

Written and pushed **before any tree was computed on a split alignment**. Part A (REPORT.md) found:

- Within a draw, one point of SPFP (false homology pairs) costs FastTree about 0.15 nRF points; one point of
  SPFN (missed pairs) costs about 0.03. A false pair is about 5x as expensive as a missed one.
- The merge-step filters (es4, recipe, hard-bb) remove about 13 SPFN points but only about 1 SPFP point,
  which the prices turn into about −0.3 to −0.5 nRF, as observed.
- In the hard-bb merge, the false pairs sit mostly on edges the vote model already trusts: 35–55% come from
  impure subset columns joined by correct edges (not fixable in the merge), 35–50% from wrong edges, and the
  wrong edges that survive hard-bb have posteriors of 0.5–0.99. Merged pairs held together only by edges with
  posterior ≤ 0.99 have precision of about 0.35–0.75, below the break-even precision of about 0.84 implied
  by the price ratio.

## The change (one method)

`code/splitcols.py VOTEDIR TAU`: take the hard-bb merge (vote.py `hard-bb`, the variant gcmvote pre-registered).
Inside every final column, link two nodes (subset columns) when their cross-subset GCM edge has beta-binomial
posterior > TAU. Each connected component becomes its own column; components are written next to each other,
so every sequence keeps its residue order. Nothing is deleted. Columns whose nodes are all strongly linked are
unchanged. This converts likely false merges into over-splits, which Part A says are cheap for trees.

## Tuning (training: SIMHIGH R1–R10, done before this file)

- Training draws: the gcmvote bank draws SIMHIGH R3–R10 (`_gv`), the ones with a hard-bb merge. (gcmtrees draws
  R1–R8 have no hard-bb merge; they are not used for tuning.)
- Grid TAU ∈ {0.9, 0.99, 0.999, 0.9999}. Criterion fixed in advance of any split tree: the Part A price model,
  predicted ΔnRF = 0.028·ΔSPFN + 0.151·ΔSPFP (vs the unsplit hard-bb), averaged over training draws
  (`code/tune_split.py`, `data/tune_split.jsonl`).
- Result: 0.9 → −0.12, **0.99 → −0.31**, 0.999 → −0.13, 0.9999 → +1.18. **Selected TAU = 0.99.**
  At 0.99: ΔSPFN +4.1, ΔSPFP −2.8, at most +10% columns.
- Training check (reported, not used to re-select): FastTree on hard-bb+split0.99 for the 8 training draws.

## Test (held out)

- Draws: SIMHIGH R11–R20 from the gcmvote bank that have a hard-bb FastTree row: R11, R12, R14, R15, R16, R17,
  R18, R19, R20 (R13 has no bank), plus any SIMHIGH R21–R50 replicate that has a bank tarball and a hard-bb tree
  row on a `claude/cs581-gcmvote-*` branch when the test trees start. No other draws are added afterwards.
- Trees: FastTree 2.1 `-lg -gamma -nosupport` (re-running MAGUS on two draws reproduced the logged trees
  exactly), one thread, 4 at a time. nRF vs the true tree.
- **Primary:** ΔnRF = (hard-bb+split0.99) − hard-bb, paired per draw. Mean, median, W/T/L with a 0.1-point
  tie band, two-sided Wilcoxon signed-rank. Success = mean ΔnRF < 0 with p < 0.05.
- Power is low and stated in advance: Part A's noise analysis puts the SD of paired ΔnRF at about 1 point,
  so with n ≈ 9 only effects near −1 point are detectable; the predicted effect is about −0.3. A null here
  would not rule out the predicted effect.
- Secondary (descriptive): split − MAGUS; observed vs price-predicted ΔnRF per draw (does the price model
  hold causally?); ΔSPFN, ΔSPFP, column counts; pooled training + test.
- If time allows, secondary arm: the same split applied to MAGUS's own merge (vote.py `magus`, TAU 0.99)
  vs MAGUS on the test draws.
