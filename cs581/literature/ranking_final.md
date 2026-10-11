# Ranked comparison of explored project ideas (final, 2026-10-10 06:50 UTC)

Supersedes `ranking_v1.md`. Sources: every exploration session's `cs581/<dir>/REPORT.md` (branch
`claude/cs581-<dir>`), the novelty audits (`novelty_audits.md`), the MAGUS-variant fan-out (54 replicates,
`../experiments/variants/SUMMARY.md`) and the measured end-to-end benchmark (all 14 datasets; PASTA failed on RNASim,
`claude/cs581-e2e-*`). **Nothing is chosen here**; the overnight runs (bottom) may move rows 3 and 4.

"Course" = how easily the choice is justified from CS581 material: **high** = an open question printed
verbatim in the slides or on the instructor's project list, about methods taught in lecture; **medium** =
methods taught, question ours; **low** = needs material not covered.

## Overnight update (2026-10-10 08:05 UTC): two new leads, one weakened

**New top candidate: EPA-ng large-tree bug** (`claude/cs581-epang`, audit `epang_upstream_audit.md`).
Answers the open problem in Wedell, Shen & Warnow (BSCAMPP, TCBB 2025): why EPA-ng's placement error more
than doubles above 2,000 leaves. Root cause: with more than 2,000 tips EPA-ng switches on per-rate scalers,
and its premasking code (`shift_partition_focus`) shifts the scaler buffers by `offset` instead of
`offset × rate_cats`, so fragments (queries not starting at column 0) get wrong likelihoods and falsely
confident placements. Unreported since v0.3.3 (2018); present in bioconda 0.3.8 (bundled by BSCAMPP),
reaches TIPP3-fast and PICRUSt2 (~26.9K-tip default tree). Evidence (RNASim 10K, 303 paired fragment
queries, true query alignment): stock mean delta error 1.09/0.97/0.82 at 500/1k/2k leaves, then
1.74/1.65/1.71 at 3k/5k/9k; forcing scalers on at 500-2k reproduces the jump; `--no-pre-mask` removes it;
the 5-line patch gives 0.78/0.79/0.73 at 3k/5k/9k (keeps improving); full-length queries unaffected (as
predicted). End to end: BSCAMPP with the patched EPA-ng at subtree size 9,000: mean delta 0.661 vs 0.798
for stock BSCAMPP at its default 2,000 (−17%, 1,000 fragment queries), at 1.7× wall-clock and 12.9 vs
2.7 GB. The second hunk (`|=`) only restores SIMD kernels (speed). Downstream (BSCAMPP at more sizes,
TIPP3/PICRUSt2, speed) running on `claude/cs581-epangdown`. Course: phylogenetic placement
(pplacer/EPA-ng/SEPP, Warnow lab's SCAMPP/BSCAMPP/TIPP). Novelty: unreported bug; the diagnosis is done,
so a project would be the characterisation, re-tuning and downstream impact, plus an upstream fix.

**fragml2 FINAL (01:51 UTC, stopped early; `claude/cs581-fragml2`).**
- *1000M1-HF, true alignment:* primary pipeline 23.7% FN vs RAxML-NG 24.4%. ΔFN +0.11 (3/1/4, p = 0.84, n = 8) at
  0.34× CPU. Constrained ML +0.40. Graft-only +9.1 (Holm p = 0.03). uDance 27.2% vs 22.1% for the primary on the
  same leaves.
- *EPA-ng fix:* identical at 1K; on RNASim10K-HF (n = 2), graft-only 47.7% stock vs 35.6% patched.
- Session verdict: "unclear, leaning not promising at 1K; promising only if reframed to 10K+".
- **Student's decision (02:00 UTC): the project is the MAGUS support filter.** Next steps:
  - threshold × backbone-count sensitivity;
  - a reference-free vote-mixture model (two-component binomial with per-edge exposure) to pick the cutoff.
- The student will reproduce the results and write the proposal themselves. The course policy (CS581-Fall2026.txt):
  AI-written code is allowed with disclosure; AI-written text is not.

**aproroom FINAL (01:40 UTC, `claude/cs581-aproroom`).**
- *Error decomposition:* confirms the interim. Gene-tree estimation error is 80–95% of all error, including at 1000
  species. It does not shrink with more genes (0.026 at both 1k and 10k estimated genes vs 0.000 on true trees).
- *Where methods separate:*
  - high duplication: ASTRID-Pro beats Asteroid (−0.011) and ASTRID-multi (−0.013); on gdl_1e-9_0 it is 0.058 vs
    ASTRID-DISCO 0.077;
  - high ILS: ASTRAL-Pro3 is best;
  - 100 genes: wQFM-GDL is best.
- *Scaling:* ASTRID-Pro takes 14 s at 1000 species with constant memory. ASTRAL-Pro3 (4 threads) takes more than 25 min
  at 500 species; Asteroid runs out of memory at 500; wQFM-GDL takes about 29 min at 1000 genes.
- *Hybrid:* an ASTRID-Pro guide makes ASTRAL-Pro3 5–10× faster at equal accuracy.
- *Headroom:* a better objective can gain at most about 0.01 FN. The lever is robustness to gene-tree error (up to
  0.03–0.09 FN), but wapro's weighting, the standard approach to that, did not help.

**Check-in 11 (01:40 UTC): the MAGUS fix does NOT improve trees; GCM clustering swap is not useful on its own.**
- *gcmtrees, pooled from the main session + h1–h3:* 12 complete simulated protein sets, FastTree nRF.
  - Pre-registered primary, recipe minus MAGUS: −0.15 (6/2/4 at a 0.1 band, Wilcoxon p = 0.64). SIMHIGH (n = 5):
    −0.40.
  - `linsi#es3` −0.03 (p = 0.05, but the mean is ≈ 0); hard filter −0.33 (p = 0.13).
  - Mean room (MAGUS − true alignment) is 1.6 RF overall, only on SIMHIGH. The ~5-point SP gain does not reach
    FastTree trees.
  - IQ-TREE on SIMHIGH and the last 3 datasets are still pending; the conclusion is unlikely to change.
- *gcmclust (REPORT draft, held-out n = 12):*
  - On the raw graph no clustering beats MCL.
  - With the support filter: es4 + Leiden-CPM −1.38 vs MAGUS (proteins −4.24, DNA +0.05) and es4 + MCL −1.30. The
    swap itself (CPM vs MCL on the filtered graph) is −0.08 (n.s.); over all 10 protein sets −0.52 (p = 0.01).
  - The es4 filter makes the trace 2–5× faster (27 s → 7 s).
  - Session verdict: "unclear, leaning not promising" stand-alone. It confirms the damage is in the graph.
  - FINAL (02:02 UTC, `claude/cs581-gcmclust`) confirms the draft. MCL clusters break the one-column-per-subset
    rule in 48–73% of clusters (CPM 2–20%). CPM has a cliff on DNA (γ = 0.03 → +1.6, γ = 0.05 → +39). FastTree nRF
    on 4 AliSim sets: MAGUS 0.090, es4 + MCL 0.082, es4 + CPM 0.084 (n = 4, not tested).
- Combined MAGUS picture: the support filter is the active ingredient (alignments better and the merge faster); the
  clustering and subset-aligner swaps add nothing; there is no tree gain.

**basemeth INTERIM 2 (00:58 UTC): other subset aligners lose on simulated proteins.**
- *Δ SP error vs L-INS-i subsets:*
  - SIMMOD_R1: ProbCons +0.95, MUSCLE5 +1.03, FAMSA +3.6, Clustal +7.6.
  - SIMHIGH_R1: G-INS-i −0.8, ProbCons +3.1, MUSCLE5 +3.5, FAMSA +5.2, Clustal +9.5.
- The BAliBASE gain (−0.3 to −0.5) reverses on simulated data and HomFam.
- Primary over 10 sets (Holm): no aligner beats L-INS-i. The instructor's "MAGUS with other base methods" item is not
  promising with current aligners.

**gcmtrees INTERIM (00:43 UTC; SIMHIGH R1–R4, FastTree nRF %).**
- Recipe minus MAGUS: +0.40, −1.81, −0.80, +0.90. Mean −0.33, 2 better / 2 worse, despite about −5 SP points of
  alignment gain (R1: SPFN 31.2 → 21.5).
- Room (MAGUS minus true alignment) is 3–6 RF points.
- *Correction to check-in 10's reading of the oracle diagnostic:*
  - split(MAGUS) 12.24 ≈ MAGUS 11.94, so MAGUS's own false positives are not what costs trees.
  - split(recipe) is 8.32 vs 12.34 on R1 but 8.93 vs 7.52 on R2.
  - So the diagnostic is noisy and inconclusive, not evidence that "wrong pairings drive tree error".
- Waiting on SIMHIGH R5–R8 and SIMMOD R1–R8 (helpers h1–h3).

**fragml2 INTERIM (00:34 UTC, `claude/cs581-fragml2`; 1000M1-HF, true alignment, paired vs RAxML-NG, n = 5).**
- RAxML-NG: 23.5% FN, 34.6 CPU-min.
- *Primary pipeline (FastTree backbone → patched EPA-ng → graft → RAxML-NG fast polish):* 24.8%; ΔFN +0.44
  (1/1/3, p = 0.31); about 1/3 of the CPU.
- *Constrained RAxML-NG:* ΔFN +0.02 (3/0/2); 0.36× CPU.
- *Graft only, no polish (Smirnov–Warnow style):* +9.0 (0/0/5). uDance ≈ graft-only quality on the leaves it keeps.
- **EPA-ng fix inside the pipeline:** placements are identical at ~500 tips, as expected. On RNASim10K-HF (5,040
  tips), graft-only FN is 48.8% with stock EPA-ng vs 35.9% with the patched one (−12.9 points); after a FastTree
  polish, 32.2% vs 30.6%.
- So the bug silently costs a placement-based tree pipeline ~13 points at a realistic size, and polishing recovers
  most but not all of it.
- Possible combined project: fast ML trees for fragmentary data (ties RAxML-NG at 1/3 CPU, beats uDance), with the
  EPA-ng bug as the new finding.

**basemeth INTERIM (00:25 UTC, `claude/cs581-basemeth`).** Other subset aligners inside MAGUS (merge-only, same
subsets and backbones), Δ SP error vs L-INS-i subsets:
- *BAliBASE (n = 5):* MUSCLE5 −0.47 (4/0/1), ProbCons −0.37, Clustal −0.31, FAMSA2 +0.26, G-INS-i +0.05. Gains are
  mostly SPFP.
- *HomFam (n = 3):* MUSCLE5 +2.05, ProbCons +2.27, Clustal +1.85, FAMSA2 +0.83 (all worse).
- *CPU cost relative to L-INS-i:* MUSCLE5 3.9×, ProbCons 7.7×.
- Early read: small, inconsistent gains; not a strong lead.

**Check-in 10 (00:10 UTC), interim.**
- *gcmtrees (FastTree nRF %, only 2 SIMHIGH datasets complete; the helper sessions h1–h3 are still running):*
  - SIMHIGH_R1: true 6.12, MAGUS 11.94, recipe 12.34, es3 11.63, hard 12.34.
  - SIMHIGH_R2: true 5.22, MAGUS 9.33, recipe 7.52, es3 8.73, hard 7.82.
  - So far mixed: recipe +0.40 and −1.81.
  - Oracle diagnostic: removing only the false-positive pairs (using the true alignment) gives the recipe 8.32 vs
    12.34 on R1. Wrong pairings, not missed ones, drive tree error. Not usable directly.
- *fragml2 (1000M1-HF, true alignment, 4–6 reps), mean FN and CPU:*
  - FastTree 48.1% (107 s); IQ-TREE `--fast` 37.5% (424 s); RAxML-NG 23.5% (2,081 s, n = 4);
  - constrained ML from a FastTree backbone 25.2% (743 s);
  - EPA-ng placement + graft, no polish: fixed 34.6% / stock 36.7% (about 150 s);
  - placement + RAxML-NG fast polish 24.6% (752 s);
  - uDance (on the leaves it keeps) 27.2% vs constrained ML 21.9% on the same leaves.
  - Two-step ≈ RAxML-NG at about a third of the CPU, better than uDance, IQ-TREE and FastTree, but not better than
    RAxML-NG.
  - The stock vs fixed EPA-ng gap (2.1 points) at under 2,000 tips is unexpected; it awaits the session's
    explanation.

**aproroom INTERIM (23:45 UTC, `claude/cs581-aproroom`): almost no headroom between GDL summary methods.**
- *Data:* FastMulRFS, 100 taxa, 100 genes, 60 reps.
- *FN rate by input:* true trees + true tags / true trees / RAxML 100 bp / 25 bp:
  - ASTRID-Pro 0.009 / 0.008 / 0.042 / 0.106;
  - ASTRID-DISCO 0.009 / 0.009 / 0.041 / 0.104;
  - Asteroid – / 0.010 / 0.040 / 0.105;
  - ASTRAL-Pro3 – / 0.008 / 0.048 / 0.130;
  - wQFM-GDL – / 0.007 / 0.040 / 0.116.
- *Where the error comes from:* tagging error ≈ 0. Gene-tree estimation error is 80–95% of the total.
- *Headroom:* picking the best of 5 methods per input after the fact (an oracle) gains only 0.001 / 0.008 / 0.012
  over the best single method. On DISCO with 1000 true gene trees every method has FN 0.
- *ASTRAL-Pro3:* the true species tree scores below ASTRAL-Pro3's own tree in 60/60 inputs, so its extra error is
  in its objective, not its search.
- **Conclusion:** a new summary method cannot gain much over the best existing ones. The remaining error is in the
  gene trees. Drop ASTRID-Pro as a "beat the state of the art" project.

**wapro INTERIM (23:04 UTC, `claude/cs581-wapro`): weighting does not help ASTRID-Pro.** FastTree-SH gene trees.
- *Train, FastMulRFS reps 01–03 (n = 36), FN rate:* wQFM-GDL 0.0703 < ASTRID-DISCO 0.0753 ≈ wASTRID-Pro 0.0756 ≈
  ASTRID-Pro 0.0759 < DISCO+wASTRID 0.0765 < Asteroid 0.0782 < ASTRAL-Pro3 0.0899 ≈ DISCO+wASTRAL 0.0901.
- Every weighting, contraction or normalisation lever is within ±0.001 of unweighted ASTRID-Pro, or worse.
- *Held-out interim (9 sets), wASTRID-Pro minus method:* vs wQFM-GDL +0.008 (n.s.), vs ASTRID-DISCO +0.003, vs
  ASTRAL-Pro3 −0.0115.
- The most accurate method here is wQFM-GDL, a quartet method that is slower (~70 s median), not the ASTRID family.
- Early read: "likely not promising" for beating the best methods.

**magusgen final (22:16 UTC, `claude/cs581-magusgen`): self-soft + consensus as "one general MAGUS" is killed.**
- 15 sets: 6 DNA/RNA and 9 protein.
- *Best combination `ss:wsoft0.03:linsi&fftns2`:* DNA −0.88 (6/0/0, p = 0.031); protein −0.27 (5/1/3, n.s.); about
  +20–30% wall time.
  - Held-out proteins: BBA0154 −1.34, BBA0190 −1.61, BBA0117 +1.02, SIMMOD_R1 +2.83; SIMHIGH_R1 timed out.
- *Self-soft alone:* DNA −0.74 (6/0/0), but protein +1.28. BBA0134 is +11.07 after a 2,292 s trace: the minclusters
  trace blows up on hard proteins.
- *Mechanism:* self-soft is the DNA ingredient and consensus/support is the protein ingredient, but they do not
  combine safely.
- *Session verdict:* "kill as one variant; unclear, leaning not promising".
- *For the MAGUS proposal:* use gcmgen's `wsoft0.03:linsi&fftns2#es4` (no self-soft). It is protein-better and
  DNA-neutral.

**Check-in 9 (21:40 UTC).**
- *magusgen, 6 sets, still running.* Self-soft MAGUS alone vs MAGUS:
  - 1000L1 −0.70, 1000M2 −0.73, 1000S1 −2.30;
  - BBA0039 −0.01, BBA0067 −0.31, BBA0101 −1.22.
- Self-soft + wsoft0.03 consensus: −0.62, −0.67, −2.28 on DNA and +0.01, −1.04, −1.83 on BAliBASE. This would be the
  general MAGUS: the DNA gain comes from self-soft and the protein gain from the consensus weights.
- But alncrit already found that self-soft's DNA alignment gain gives no tree gain (tree headroom on full-length
  DNA ≈ 0.2–0.6 FN).
- *Headroom map (for choosing):*
  - full-length DNA, alignment → tree: ≤ 1 FN point, even with the true alignment (alncrit, ml);
  - hard simulated proteins: MAGUS vs true-alignment trees 2–5 RF points (protcons);
  - fragmentary data: FastTree 49% / IQ-TREE 30% / GTM 28% / RAxML-NG 25% FN on 1000M1-HF (Park et al. 2021,
    reproduced). This is the largest gap.
- *New sessions:*
  - `claude/cs581-gcmtrees`: does the recipe improve FastTree/IQ-TREE trees? 16 simulated protein sets,
    pre-registered.
  - `claude/cs581-fragml2`: backbone → (constrained ML | fixed/stock EPA-ng placement + polish) vs RAxML-NG on
    fragmentary data, pre-registered.
- *bbtool-5:* phase 2 still running (1000M2); not needed for any decision.

**astridpro2 final (20:27 UTC, `claude/cs581-astridpro2`; verdict "promising", pre-registered go criterion met).**
- *Pooled simulated data (FastMulRFS 120 + DISCO 80 runs, 100 taxa, estimated gene trees):*
  - vs ASTRAL-Pro3: −0.0135 FN (84/30/26), Holm p = 6e-8, n = 140. On DISCO alone −0.004 (8/9/3, p = 0.054, n = 20).
  - vs ASTRID-multi: −0.0034 (p = 0.023, Holm 0.068); on high-duplication DISCO −0.014 (20/4/3, p = 0.003).
  - Ties ASTRID-DISCO and Asteroid; beats DISCO+ASTRAL and FastMulRFS (−0.024); ties wQFM-GDL.
- *Speed:* 100 taxa × 1,000 genes in 1–40 s vs 2.5–28 min for ASTRAL-Pro3, which timed out on 4 of 24 heavy runs.
  14 s at 1,000 taxa.
- *Theory:* proved (see check-in 7).
- *Caveats:* 500-gene runs not reached; ASTRAL-Pro3 on DISCO covers reps 01–02 only; ASTRAL-Pro's own Dryad gene
  trees were blocked (403).

**protcons final (20:17 UTC, `claude/cs581-protcons`; pre-registered; 28 held-out protein sets: 8 simulated, 10
HomFam, 2 10AA, 8 fresh BAliBASE draws).**
- *Pre-registered primary `linsi|cons0.7` + gate:* FAILS. −0.59 (15/5/8, p = 0.13).
  - BAliBASE replicates the pilot (−2.03, 7/1/0), but excluding BAliBASE it is −0.01.
  - The gate's protein accuracy is 54% (chance).
- *MAFFT-only `L ∩ FFT-NS-2 --op 3`:* −1.28 (23/1/4, p = 0.0003); excluding BAliBASE −0.97 (p = 0.007); simulated
  8/8 (−2.28); HomFam −0.15; nucleotides +4.8 (0/0/5).
- *Post-hoc edge support ≥ 5 of 10 backbones* (k chosen from 3 values on the same data):
  - proteins −1.90 (17/4/7, p = 0.004); simulated −5.81 (8/0/0, via SPFN −11.6); HomFam −0.11; BAliBASE −0.71;
  - nucleotides −0.03.
  - This independently agrees with gcmgen's edge-support result.
- *Trees (FastTree, 8 simulated sets):* Δ nRF +0.36 (primary), −0.11 (L∩F), all p > 0.4. **No tree gain.**
- *Runtime:* all recipes cost < 5% of MAGUS on proteins.
- *Session verdict:* consistency filtering as pre-registered is not promising. The pivot, support-aware GCM graphs
  (edge-support threshold/weighting, plus L∩F on proteins) pre-registered on fresh data, is "unclear to promising".
  The suspected mechanism is MCL fragmentation from low-support edges.
- *Consolidated picture with gcmgen:*
  - Edge support is the general lever: it helps proteins (mostly simulated, via recall) and is neutral on DNA.
  - The FFT-NS-2 confirmation adds a BAliBASE (precision) gain.
  - Neither helps DNA/RNA, and neither yet improves trees.

**bbtool-6 final (20:05 UTC; 16S.M and 1000L1, 3 draws each; Δ vs MAGUS's own merge).**
- *16S.M:* Clustal +5.4 to +7.1. `mafft --auto` backbones −0.39 / −0.84 / −0.32. Union with Clustal +0.1 to +0.3.
  L-INS-i without `--ep` and G-INS-i ≈ 0.
- *1000L1:* Clustal +23 to +25. `--auto` +40 to +44. Union −0.01 / −0.23 / −0.11. L-INS-i without `--ep` +0.4 to
  +1.7. G-INS-i −0.01 / −0.32 / −0.55.
- *Conclusion unchanged:* Clustal backbones are protein/BAliBASE-only.
- *Caveat:* the 16S.M state rows carry `"protein": true` (the old type-check bug). It is a recorded label only:
  `bbtool_bench` does not use it to choose any tool or option.

**gtmscale final (20:04 UTC, `claude/cs581-gtmscale`): GTM blending is NOT PROMISING as a better DTM.**
- Blending vs GTM over 46 cases: −0.45 FN points (CI −0.66 to −0.27), 32/7/7, p = 3e-6.
- IQ-TREE `-fast` is more accurate up to 5K taxa at similar time. At RNASim 10K, blending only ties FastTree
  (10.42), though IQ-TREE ran out of memory there.
- One re-decomposition gains 2–6 points; blending gains ~0.5. Hard constraints lock in subset-tree errors
  (an unconstrained FastTree polish beats blending with poor guides).
- Only "unclear" if reframed (soft constraints plus iteration, or an analytical negative result).
- Row 5 is downgraded to not promising.

**gcmgen final (20:02 UTC, `claude/cs581-gcmgen` REPORT.md; 24 replicates, one MAGUS draw each; recipe fixed on 10
training replicates).**
- *Recipe `wsoft0.03:linsi&fftns2#es4`:*
  - Proteins −2.76 (9/1/0, p = 0.002): BAliBASE −1.08, AliSim −5.28.
  - DNA/RNA +0.01 (6/4/4), worst +0.47 (RNASim_R1).
  - Held-out only: proteins −4.17 (4/0/0), DNA/RNA +0.06 (4/3/3).
  - Extra cost: one FFT-NS-2 run per backbone (~1% of L-INS-i).
- *Edge support alone (`linsi#es3`):* proteins −2.03 (6/3/1), DNA/RNA −0.03 (worst +0.05), BAliBASE ≈ 0.
- *Second opinions:*
  - G-INS-i is reliable on ROSE (hard filter +0.87 vs +32.5 with FFT-NS-2) but weaker on BAliBASE.
  - L-INS-i with a random guide tree is ~600× the cost for no gain.
- *Reference-free switch:* it only detects ROSE and loses on RNASim; the no-switch recipe is better.
- *Session verdict:* "yes, a general recipe: clearly helps proteins, neutral (not helpful) on DNA/RNA within ±0.5".
  Caveats: single draws; 4 simulated and 6 BAliBASE replicates; RNASim loses ~0.4.
- *Still missing:* end-to-end timing and tree accuracy for the recipe (protcons trees for the hard filter are within
  noise).

**Check-in 8 (19:30 UTC).**
- *ASTRID-Pro, empirical (`claude/cs581-astridpro2`; my paired recomputation from `results/{disco,fmrfs}_runs.jsonl`;
  FN rate, ASTRID-Pro minus the other method, negative favours ASTRID-Pro):*

  | vs | n | Δ FN | W/T/L | p |
  |---|---|---|---|---|
  | ASTRAL-Pro3 | 71 | −0.0113 | 46/14/11 | 2e-6 |
  | ASTRID-multi | 113 | −0.0042 | 44/41/28 | 0.019 |
  | DISCO+ASTRAL | 58 | −0.022 | — | 5e-7 |
  | FastMulRFS | 57 | −0.021 | — | 2e-6 |
  | DupLoss-2 | 54 | −0.108 | — | — |
  | wQFM-GDL | 54 | −0.004 | — | 0.09 (tie) |
  | ASTRID-DISCO | 113 | +0.001 | — | 0.36 (tie) |
  | Asteroid | 106 | −0.0006 | — | 0.26 (tie) |

  - The ASTRAL-Pro3 gain is mostly on the FastMulRFS (= ASTRAL-Pro S100) data: −0.0132, n = 57. On DISCO it is
    −0.0036, n = 14, n.s.; ASTRAL-Pro3 timed out on the 3 heaviest DISCO conditions.
  - Our ASTRAL-Pro3 reproduces the published ASTRAL-Pro trees (0.0855 vs 0.0863).
  - Runtime: median 0.3 s vs 12.9 s for ASTRAL-Pro3. On 1000-gene DISCO inputs it takes 1–40 s vs 6–29 min. At
    1000 taxa it takes 14 s, while ASTRAL-Pro3 and FastMulRFS take > 20 min at 500 taxa.
  - Empirical 1KP: 5 FN, equal to ASTRAL-Pro3, vs 10 for ASTRID-multi.
  - Session verdict: "promising; Go by the pre-registered criterion" (matches ASTRAL-Pro3 and is ≥ 10× faster).
  - Honest framing: a provably consistent (true tags), very fast distance method that matches or beats ASTRAL-Pro3,
    but ties the fast heuristics ASTRID-DISCO and Asteroid, which have no GDL proof.
  - **This upgrades the GDL line from "theory only" to "theory plus a faster, at-least-as-accurate method".**
- *gcmgen held-out, now n = 10–12:*
  - The selected recipe `wsoft0.03:linsi&fftns2#es4`: −0.82 overall (7/1/2). Proteins −2.77 (3/3); DNA/RNA +0.02
    (4/1/2), worst +0.34.
  - Edge support alone `linsi#es3`: −0.44 (4/8/0, p = 0.009), worst +0.04. Its protein gain is mostly SIMMOD_R2
    (−4.68); on the BAliBASE held-out sets it gives only −0.03 and −0.22.
  - Held-out proteins are only n = 3 (SIMHIGH_R2 pending).
- *protcons:*
  - Fresh BAliBASE draws so far: L∩FFT `--op 3` −0.10 (BBA0039) and −0.38 (BBA0067); `cons0.7` −1.87 on BBA0067.
  - 1000M2_R1: `cons0.7` +10.25 (DNA harm again).
  - FastTree trees on the 8 simulated protein sets: L∩FFT vs MAGUS RF −0.60, −0.90, +0.60, +0.50, +0.10, 0, −0.60, 0
    (mean −0.11). **No reliable tree gain**, within FastTree noise.
- *magusgen (5 sets):* self-soft alone is ≤ 0 on all five (1000L1 −0.70, 1000M2 −0.73, BBA0039 −0.01, BBA0067 −0.31,
  BBA0101 −1.22). Self-soft + wsoft0.03 on BBA0101 −1.83, vs −1.22 for self-soft alone.
- *gtmscale RNASim 10K (1 rep):* with a good FastTree guide, GTM, Blend and full FastTree are all within 0.1 FN. With
  a k-mer guide, Blend repairs GTM (−0.75), back to full-FastTree level. Same picture: blending only matters when
  the guide tree is bad.
- *bbtool phase 2:* MAFFT-only backbone modes (`--auto`, L-INS-i without `--ep`, G-INS-i) give noise-level changes on
  BBA0154/0190 (−0.7 to +1.3) and −0.5 to −1.9 on BBA0081/0101. Clustal backbones: −0.7 to −7.9 on all 12 BAliBASE
  draws, but catastrophic on DNA. Conclusion unchanged.

**19:05 UTC: a general recipe for GCM evidence (`claude/cs581-gcmgen`, report draft; train/held-out split fixed
in `SPLIT.md` before any results; no data-type switch).**
- *Mechanism, on every data type:* cross-subset GCM edges supported by fewer than 4 of the 10 backbones are almost
  all wrong. On BAliBASE, the precision of their evidence is 0.02–0.04 at support 1 and 0.05–0.10 at support 2–3, vs
  0.76–0.88 at support ≥ 4.
- *Recipe selected on 10 training sets:* `wsoft0.03:linsi&fftns2#es4`. L-INS-i pairs that FFT-NS-2 does not confirm get
  weight 0.03, then edges with support < 4 are deleted.
  - Train: −1.14 (7/2/1, p = 0.014); proteins −1.82; DNA/RNA −0.11; worst +0.41.
  - Held-out so far (n = 5): −1.60 (4/0/1). Proteins: BBA0154 −1.22, BBA0190 −1.98, SIMMOD_R2 −5.11. DNA/RNA:
    RNASim +0.34, 1000M4 −0.05.
- *The edge-support threshold alone (`linsi#es3`, no second aligner, MAGUS's own evidence only):*
  - held-out −0.87 (4/2/0, p = 0.031); no held-out set got worse (worst −0.03);
  - train proteins −1.27 with `#es3`, −1.56 with `#es4`.
- *This corrects check-in 7's "edge support ≈ 0":* thresholds 2–3 are ≈ 0 on DNA, but thresholds 3–5 help proteins.
- Still running: more held-out DNA/RNA replicates (1000S3, 1000M2_R1, 1000L1_R1, RNASim_R1) and SIMHIGH_R2. If it
  holds, this is the first MAGUS-line change that is better on proteins and neutral on DNA/RNA with one setting.

**Check-in 7 (17:30 UTC).** Δ = error points vs MAGUS's own merge on the same subsets and backbones (negative is better).
- *Consensus GCM evidence, held-out so far (`claude/cs581-protcons`):* MAFFT-only `L-INS-i ∩ FFT-NS-2 --op 3`:
  - simulated proteins: −2.27, 8/0/0, p = 0.008;
  - HomFam (10 families): −0.15, 6/1/3, p = 0.43 (Acetyltransf −1.60, aat −1.00, PDZ −0.65 … blmb +2.09, p450 +0.43);
  - 10AA: −0.14, +0.33;
  - all 20 held-out sets: −0.97, 15/1/4, p = 0.007.
  So it is real on simulated proteins, flat on HomFam (recall-limited there), and harmful on DNA. Still to come: fresh
  BAliBASE draws, trees, nucleotide controls.
- *Generalizing it to DNA (`claude/cs581-gcmgen`, one draw per set):* instead of deleting backbone pairs that FFT-NS-2
  does not also align, keep them at weight 0.03 (`wsoft0.03`).
  - DNA/RNA (8 sets): mean −0.08, 3/3/2, range −0.64 (16S.M) to +0.23. This makes it safe on DNA (the hard filter is
    +20 to +31 on ROSE).
  - BAliBASE (6 sets): mean −1.05, 6/0/0, including BBA0154 −1.64 and BBA0190 −1.84, which are in its held-out split.
  - But SIMMOD_R1 is +1.12, where the hard filter gave −2.12. The soft weight gives up the simulated-protein gain.
  - Edge-support thresholds (`#es2/3`) and duplicated backbones (`dup2`) ≈ 0 everywhere.
  - The session's pre-registered adaptive switch (train/held-out split in `SPLIT.md`) is still running.
- *General MAGUS (`claude/cs581-magusgen`, 3 sets so far):* self-soft + soft consensus 1000L1 −0.65, BBA0039 −0.05,
  BBA0101 −1.02; self-soft alone −0.70 / −0.01 / −1.22. So far the combination adds nothing beyond self-soft.
- *GTM blending at scale (`claude/cs581-gtmscale`), FN points:*
  - n = 2,000 with a poor FastTree guide (56% FN): Blend-FT vs GTM −0.75 (5/0/1, p = 0.06); it beats full FastTree by
    6.4 but loses to full IQ-TREE by 3.0.
  - Published conditions: 1000M1-HF −0.62 (FT guide) and 0.00 (IQ guide); Cox1-HET −0.22 (7/2/1, p = 0.023) and −0.05.
  - So the large round-1 gains (−6 to −14) need a bad guide tree; with the guides used in practice, blending adds
    ≤ 0.75. Downgrade to "modest". RNASim 10K still running.
- *ASTRID-Pro theory (`claude/cs581-astridpro2`, interim; empirical part running):*
  - Exact limit with true trees and tags: additive on the species tree with closed-form edge lengths
    β(c) = ½[s(c₁) + s(c₂) + s(sib) − s(c)], where s(·) is a survival probability.
  - Consistent iff every interior β > 0. β > 0 whenever the branch above c is not supercritical. A counterexample
    with β = −0.267 matches simulation to 0.003.
  - Survival reweighting is consistent under any rates, given a correct first-pass tree.
  - Every DISCO/FastMulRFS benchmark condition has β ≥ 0.46, so the failure does not occur in published benchmarks.
    That is the same "theory, low practical reach" profile as ASTRAL-Pro.

**ASTRAL-Pro under GDL, round 2 final (16:15 UTC, `claude/cs581-astralpro2`, proofs in `cs581/astralpro2/theory.md`).**
Updates row 1: the "numerical, not proofs" caveat is now largely removed.
- *Proved* on the 4-taxon caterpillar (((A,B),C),D):
  - closed-form limiting ASTRAL-Pro scores with the correct root and ASTRAL-Pro's own species-overlap tags;
  - a sign law: at a duplication with unequal copy numbers on its two sides, the wrong pairing wins iff the
    outgroup-side taxon C survives better than B. Corollary: if c ≤ min(a, b), it is consistent for any process on
    that branch;
  - a rigorous inconsistency instance (pure-birth λT = ln 20, survival a = b = 0.01, c = 0.2): correct 2.4025e-4 vs
    3.2663e-4 per wrong topology (50-digit interval arithmetic);
  - re-tagging by reconciliation against any tree T makes ASTRAL-Pro return T, so "re-tag against a first-pass
    tree" cannot fix it.
- *Not proved:* ASTRAL-Pro3's own min-score rooting makes it 2.6× worse (simulation certificate, z = −117); it turns
  7 of 30 consistent cells into failures.
- *Methods (10 blocks of 1,000 families):*
  - ASTRAL-Pro3, DISCO+ASTRAL and wQFM-GDL: wrong 10/10 (9/9 on a third pool);
  - duplication+loss parsimony and DupLoss-2: 0/10 wrong;
  - with FastTree-estimated gene trees, ASTRAL-Pro3 is wrong 5/5.
- *Practical reach is low:* no failures at published rates (S25, λH ≈ 0.93). Failures need λH ≥ 4 with ~3×
  per-branch rate heterogeneity, and then affect 0.2–0.4% of quartets.
- *Verdict:* "promising (theory project)". It gives a negative answer to slide 35's "Is ASTRAL-Pro consistent for GDL
  …?" Course: high. Risk: low practical relevance. Overlap check (16:25 UTC, Parsons, Liu, Dua, Markin & Molloy, bioRxiv 10.64898/2026.01.20.700722, v2 of Apr 12 2026, full text read): they work under DLCoal with *correct* tagging (their new definition), conjecture consistency (Conjecture 1, "an open question"), describe an adversarial scenario that breaks the exchangeability argument, and leave specific rates λ, μ to future work. No inconsistency result, nothing with ASTRAL-Pro's own tagging. So our result does not overlap; it complements theirs (with ASTRAL-Pro's own tags the answer is negative, even without ILS). Prop 5 (reconciliation against T returns T) is probably folklore; cite it as an observation.

**Interim, 15:40 UTC (round-2 sessions still running).**
- *Consensus GCM evidence on HomFam (`claude/cs581-protcons`):* held-out, so it counts against the proposal.
  - MAFFT-only `L ∩ FFT-NS-2 --op 3` is mixed: Acetyltransf −1.60, PDZ −0.65, aat −1.00, adh 0.00, blmb +2.09,
    p450 +0.43; the four other families are within −0.44 to +0.13.
  - Primary `linsi|cons0.7` hurts on four families (PDZ +3.54, blmb +3.42, p450 +1.10, aat +0.79). The losses come
    through SPFN (blmb +8.9), consistent with MAGUS on HomFam being recall-limited, and the gate said "filter" there.
  - 10AA: every variant within ±0.33.
  - So far, the effect holds on simulated proteins and BAliBASE but not on HomFam. Fresh BAliBASE draws and the
    nucleotide controls are still running.
- *magusgen (general MAGUS):* hard filters are catastrophic on 1000L1 DNA (consistency mask +9.2, SPFN 29); self-soft
  −0.70 there, as before; soft-weighted consensus variants queued next.
- *bbtool-5, nucleotide controls (17:12 UTC; 3 MAGUS draws each; Δ vs MAGUS's own merge, same subsets):*
  - Clustal backbones: RNASim +0.23 to +0.88; 1000M2 +17.9 to +19.9 (e2e +13.8 to +15.5). They are unusable on DNA.
  - Union of 10 L-INS-i + 10 Clustal backbones: RNASim −0.41 / −0.45 / −0.50 (3/3); 1000M2 +0.01 / +0.38 / +0.22.
    The RNASim gain may just come from having 20 backbones instead of 10. "20 L-INS-i backbones" is the control
    that separates the two.
- *gcmgen:* a soft down-weight (w = 0.03 on pairs without consensus) keeps the BAliBASE gain; it is being run on all datasets.

**Check-in 6 (15:30 UTC): DecoDiPhy identifiability, round 2 final (`claude/cs581-decodiphy2`, proofs in
`cs581/decodiphy2/theory.md`).** Source: Arasti, Şapcı, Rachtman, El-Kebir & Mirarab, RECOMB 2026, which states
"we suspect (with no proof) that conditions that break identifiability require extreme cases of symmetry … We
leave a full characterization of identifiability to future work."
- *Proved:* k = 1; k = 2 (non-adjacent: identifiable against every alternative with k' ≤ 2; adjacent: the edges,
  node measure and ȳ are identifiable, but (p, x) form a 1-dimensional continuum for every parameter value);
  k = 3 non-adjacent: non-identifiable ⇔ continuum ⇔ "closed claw"; more than n/2 non-adjacent placements are
  never identifiable.
- *General k:* a Hall-type criterion (some set U of internal nodes with |N[U] ∖ V(S)| < |U|). The direction
  "Hall deficiency ⇒ continuum" is proved; the converse held in 145,543 checks. The full equivalence is a
  conjecture with 0 counterexamples among ~140k MILP-certified cases (n ≤ 12, k = 4–5).
- *Novelty:* the supplement (SB.3) proves only the adjacency continuum, the x ∈ {0, 1} ambiguity, and k = 2 edge
  identifiability for p = (½, ½). Our k = 2 generalization, the k = 3 theorem, the Hall criterion and the refutation
  of the symmetry conjecture are not in it.
- *Practical consequence, on the authors' stored runs (4,171 runs, 12 trees):* 11.4% of noise-free fits (24.8% at
  k = 10) have a non-trivial equivalence class. On fits with exactly the true edges, abundance L1 error is 0.0003
  when the class is trivial vs 0.119 when it is not. Flagged fits carry 97% of all noise-free abundance error
  (30–38% with noise). Reporting the class as intervals covers the truth in 82–90% of flagged runs (0% for the point
  estimate); a central point estimate is no better (p ≥ 0.2). Post-processing cost is negligible.
- *Verdict:* "promising" as theory plus uncertainty-aware output, not a point-accuracy win. Course: medium
  (tree metrics and distance methods are taught; mixture deconvolution is not). Risks: the general-k proof may not
  close in 4 weeks (fallback: k ≤ 3 theorem plus conjecture); the strongest numbers are noise-free.

**Check-in 5 (13:50 UTC).** Held-out simulated proteins complete (`claude/cs581-protcons`, 8 AliSim replicates,
pre-registered): MAFFT-only `L-INS-i ∩ FFT-NS-2 --op 3` −0.89 / −4.17 / −3.25 / −1.60 (SIMHIGH R1-R4) and
−0.97 / −3.79 / −2.51 / −1.01 (SIMMOD R1-R4): 8/8 better, mean −2.27, p = 0.008; primary `linsi|cons0.7` mean −0.84
(5/8); `(L + Clustal)|cons0.7` worse (+0.3 to +12.1 on 6/8). With the in-sample BAliBASE pilot (−1.83, 7/1/0) this
makes "consensus evidence for GCM" (MAFFT only) the strongest MAGUS-line candidate; ~2 points is comparable to
MAGUS's own gain over PASTA (2.67). Pending: HomFam, fresh BAliBASE draws, tree accuracy, nucleotide controls
(filtering is expected to hurt there). Draft: `cs581/proposal/proposal_gcmcons.md`.

**Check-in 4 (12:20 UTC).**

*Clustal-backbone confirmation, BAliBASE final (`claude/cs581-bbtool-1..4`, 8 RV100 sets × 3 fresh MAGUS draws = 24 paired
runs)*: merge-only −1.50 points (median −0.95, 20/0/4, p = 0.0025); end-to-end MAGUS with Clustal backbones −1.75
(20/0/4, p = 0.0006) at 3.26× mean speed-up (1.85-4.47×). Per set: BBA0081 −7.95, BBA0154 −1.94, BBA0067 −1.31,
BBA0190 −1.13, BBA0101 −0.99, BBA0039 −0.21, BBA0117 +0.47, BBA0134 +1.04. MAFFT-only backbone swaps do not reproduce
it (`--auto` −0.29 n.s.; L-INS-i without `--ep` −0.22, G-INS-i −0.22, n = 12, n.s.); union L + Clustal −0.98
(20/1/3, p = 2e-5). But outside BAliBASE (protbench, 15 datasets) Clustal backbones hurt: +1.85 (1/3/11, p = 0.005),
worst on simulated proteins. So: better and 3× faster on MAGUS's own protein benchmark, not general.

*Consistency-filtered GCM evidence, held-out test (interim, `claude/cs581-protcons`; pre-registered; Δ error points vs
MAGUS's own merge, same subsets/backbones)*: SIMHIGH R1/R2, SIMMOD R1/R2, 10AA coli_epi.
- Primary `linsi|cons0.7`: +0.27 / −1.99 / +0.14 / −1.63 / +0.03 (mixed).
- MAFFT-only `L-INS-i ∩ FFT-NS-2 --op 3`: −0.89 / −4.17 / −0.97 / −3.79 / +0.33 (better on all 4 simulated replicates).
- `(L + Clustal)|cons0.7`: +11.71 / +3.57 / +1.51 / −2.75 / +0.11 (Clustal hurts simulated proteins again).
- The pre-registered gate (support of L-only pairs < 0.615 → filter) said "do not filter" everywhere (support 0.66-0.87),
  including where filtering helped by 2-4 points: the gate does not transfer. Edge-support threshold (≥ k of 10
  backbones) ≈ 0 on 10AA. Still to come: SIM R3/R4, trees, HomFam, fresh BAliBASE draws, nucleotides.

*EPA-ng final (`claude/cs581-epang`)*: 4 RNASim replicates × 1,000 fragments: stock BSCAMPP b2000 0.860, stock b5000
1.641, patched b5000 0.757 (−12%, p = 0.03), patched b9000 0.721 (−17%, p = 0.03, 13.8 GB; 1 of 4 runs OOM); 94% of
queries tie. Full 9K tree: stock 1.539 (50 s) vs patched 0.661 (36 s). Overconfidence: LWR ≥ 0.9999 on 59-62% of
fragments (stock, > 2k tips) vs 18% patched. Verdict "promising as diagnosis + 5-line fix + corrected BSCAMPP design".

*PICRUSt2 with metagenome ground truth (`claude/cs581-picrust`, the PICRUSt2 paper's mammal and ocean validation
sets)*: placements change a lot (mammal: 94 of 323 ASVs change edge, 73 change closest reference; ocean: 85 of 1,148,
9 change domain), but predicted KO/pathway accuracy vs the paired metagenomes is unchanged (KO Spearman 0.7816 vs
0.7819 mammal, 0.8088 vs 0.8073 ocean; precision/recall within 0.003). Hidden-state prediction smooths placement
errors: no downstream harm to PICRUSt2 accuracy.

*Other finals*: prost3di "promising (narrowly)" for small low-identity MSAs (RV11 +0.048 SP, p = 0.002; RV12 −0.027;
~300× L-INS-i CPU), "not promising as MAGUS evidence". iqstop "NOT PROMISING". PASTA reruns: RNASim 10.10% (published
10.08%), 16S.M 14.06% (published 12.99%).

**Check-in 3 (10:45 UTC).**

*EPA-ng downstream, final* (`claude/cs581-epangdown`): on nt78 (77K leaves) stock BSCAMPP 1.72 / 4.77 / 5.09
at 2k / 5k / 10k vs fixed 1.72 / 1.70 / 1.69; RNASim 50K 0.53 / 1.43 vs 0.53 / 0.51 (fix vs stock at 5k:
p = 1e-80, 1e-55). Bug 2 alone = all of the accuracy (matches the full fix on 98.3% of queries). But fixed
BSCAMPP at 5k-10k is **not better than stock at its default 2k** on these larger sets (−0.02, n.s.); only the
RNASim 10K R1 run showed a small gain (−0.125, p = 0.005). Standalone speed: the fix is 3.3× / 2.9× faster on
fragments at 5k / 10k tips (1.6-1.7× full-length). PICRUSt2 2.6.3 (12K-tip subset of its 26,868-tip tree, 749
V4 ASVs): 2.7% of ASVs change edge; 5 ASVs at NSTI ~32 (stock) vs ~1.3 (fixed); predicted-metagenome change per
sample median 0.18%, max 16%. Session verdict: "unclear, leaning not promising as a downstream-accuracy
project". Net: the EPA-ng project is a diagnosis + correctness + 3× speed story, not a BSCAMPP accuracy story.

*EPA-ng whole-tree result (11:30 UTC, epangdown follow-up)*: on 10K (nt78) and 8K (RNASim) sub-backbones with the
same 1,000 fragments, stock whole-tree EPA-ng 2.55 / 1.72 vs fixed 0.78 / 0.78 (p = 2e-82, 6e-52) vs BSCAMPP(e)
b2000 0.79 / 0.81; fixed whole-tree is faster (24 s vs 43-69 s), memory-limited (12-13 GB at 8-10K leaves). The
BSCAMPP paper's Exp. 5 conclusion (whole-tree EPA-ng "much more" error-prone) is the bug. BSCAMPP(p) vs fixed
BSCAMPP(e) on nt78: −0.016 (n.s.) at b2000, −0.064 (p = 0.006) at b5000, pplacer 2.6-3.1x slower. Revised
session verdict: "unclear, but better than first stated".

*Consistency-filtered GCM evidence* (`claude/cs581-bbevidence`, new lead, MAGUS line, MAFFT-only): on proteins
GCM is precision-limited (Δ evidence precision vs ΔSPFP ρ = −0.66, 182 pairs). One draw per BAliBASE set:
L-INS-i backbones with columns of cross-backbone consistency < 0.7 masked: −1.97 points (7/0/0, p = 0.016; no
new alignment, ~10-20 s); L-INS-i ∩ FFT-NS-2 `--op 3`: −1.83 (7/1/0); L ∩ Clustal −2.02; (L + Clustal) masked
−2.53 (7/0/1). Filtering hurts nucleotides (1000M2 +9.5 to +41); a reference-free gate was fit in-sample.
Pre-registered held-out test (simulated proteins with trees, HomFam, 10AA, fresh BAliBASE draws, nucleotide
controls) running on `claude/cs581-protcons`.

*Clustal backbones, protbench final*: outside BAliBASE Clustal backbones hurt: mean +1.85 points, 1/3/11,
p = 0.005 (simulated proteins +1.6 to +4.8; HomFam mixed, PDZ +7.0); end to end +1.44 (3/1/11), though 2.35×
faster. Dead as a general method.

*Others*: prost3di (predicted 3Di): RV11 (low identity) +0.048 SP over L-INS-i (48/3/25, p = 0.002) but RV12
−0.027 (24/0/64); costly (ProstT5 ~60 ms/residue); no verdict written; MAGUS-evidence effect within noise.
iqstop: KH stopping rules lie on the same speed/lnL frontier as plain `-nstop N` (no better); on empirical 16S
every faster rule loses lnL. WITCH-lite: "unclear, leaning not promising". decodiphy: promising (theory).

**EPA-ng follow-up (09:15 UTC, `claude/cs581-epang`, `claude/cs581-epangdown`).** Replicated on RNASim R1:
stock BSCAMPP 0.846 at 2,000 vs 1.682 at 5,000 (the paper's 2× jump), fixed 0.721 at 5,000 (p = 0.005).
On two more datasets (1,000 fragment queries each, BSCAMPP subtree sizes 2k/5k/10k): nt78, stock 1.72 /
4.77 / 5.09 vs fixed 1.72 / 1.70 / 1.69 (fix vs stock at 5k: p = 1e-80); 16S, stock 25.3 / 28.0 vs fixed
25.4 / 26.9 (5k: p = 0.009). Bug 2 alone carries the whole accuracy effect; bug 1 alone is speed only
(fixed EPA-ng 18-22% faster on 5k-10k subtrees). So the robust results are: the anomaly is explained, the
fixed placer is flat in subtree size, and it is faster; "larger subtrees are more accurate" holds on
RNASim 10K (−15 to −17%) but not on nt78 or 16S (≈ 0). The fix matters most where EPA-ng runs on whole
large trees with fragmentary queries (PICRUSt2's amplicons on a ~26.9K-tip tree; BSCAMPP paper Exp. 5).

**Distance-mixture deconvolution theory (`claude/cs581-decodiphy`, final).** Exhaustive LP checks (all
tree shapes to n = 12 for k ≤ 2, n = 9 for k = 3) plus ~54k random instances, zero exceptions: k = 1 and
non-adjacent k = 2 identifiable; adjacent k = 2 always a continuum; k = 3 non-identifiable exactly when a
"closed claw" exists. Refutes the paper's "needs extreme symmetry" conjecture (28% of generic k = 3 claw
cases; explicit n = 5 counterexample). Theorem: d determines the node measure modulo weighted-Laplacian
moves. k-selection: learned rule 0.43 vs 0.33 exact-k (p = 5e-18) but placement barely changes. Verdict:
promising as a theory project; source is a Mirarab-lab RECOMB 2026 paper; course link via tree metrics.

**Clustal Omega backbones: real on BAliBASE, not general.** Measured end to end with the threading fix,
MAGUS with Clustal backbones is 3-4× faster on BAliBASE (e.g. BBA0067 254 s vs 985 s; BBA0154 268 vs 940 s)
and more accurate on BBA0039 (−0.14 to −0.36, 3 draws), BBA0067 (−0.8 to −1.2), BBA0154 (−1.9); cached
inputs also BBA0101 (−2.8), BBA0190 (−1.2); worse on BBA0117 (+1.3) and BBA0134 (+1.5). But it does not
generalize: simulated proteins (AliSim LG+G4 with indels) +4.3 to +5.1 points worse, 10AA 1GADBL ≈ 0,
HomFam aat ≈ 0 paired / +0.8 end to end, and 16S.M (DNA) +7.0. `mafft --auto` backbones are mixed
(−1.5 on BBA0067, +0.1 to +0.25 elsewhere, −0.4 on 16S.M). Row 3 below is downgraded to "dataset-dependent".

**WITCH-lite** (`claude/cs581-witchlite`): BLAST-guided HMM selection matches WITCH on 16S (+0.01-0.02 SPFN
at 22-26% of the time) but fails on divergent ROSE (+10 SPFN); best divergent-data variant (beam-2
descent) +1.1 SPFN at 46% of the time. Leaning not promising.

Still running: prost3di, iqstop, decodiphy (theory: explicit k = 3 non-identifiable instance found),
bbevidence (mechanism), the bbtool confirmation draws, PASTA reruns on RNASim/16S.M.

| rank | idea (branch) | strongest result | beats strongest baseline? | runtime | novelty (audit) | course | verdict |
|---|---|---|---|---|---|---|---|
| 1 | **ASTRAL-Pro is inconsistent under GDL because of its own rooting/tagging** (gdlcons) | Exact 4-taxon formula; with true gene trees stock ASTRAL-Pro3 returns the wrong species tree in 40/40, 20/20, 4/4 datasets (500 / 2k / 10k families); correct with true tags. Holds **even with constant rates** (λ = μ on every branch), but only at extreme turnover (λ = μ = 8: wrong 20/20 and 4/4; never at λ ≤ 4). Tag-error threshold q* = 0.080 ± 0.005 observed vs 0.079 predicted | n/a (theory + simulation) | cheap (simulation) | core result novel: counterexamples to Zhang et al.'s consistency conjecture; partly known context | **high**: slide 35 asks verbatim "Is ASTRAL-Pro consistent for GDL under ... error for rooting and tagging?"; GDL lecture "Phylogenomics, part 2" | promising (theory-led). Caveats: numerical, not proofs; constant-rate failures need extreme rates; wQFM-GDL probably inherits the issue (unverified) |
| 2 | **ASTRID-Pro: a GDL-corrected internode distance** (gdl) | Theorem: ortholog-only internode distance is a tree metric with the species-tree topology iff no supercritical branch; counterexample matches to 3 decimals; held-out DISCO test on high-GDL data vs ASTRID-multi −0.0096 FN (16/12/5, Holm p = 0.021); ties ASTRID-DISCO / ASTRAL-Pro | only on high-GDL data; ties elsewhere | cheap | NOVEL (narrowly); do not claim ASTRID-multi is answered | **high**: slide 35 verbatim ("distance correction for GDL so ASTRID/NJst are consistent?"); ASTRID taught in ILS lectures | promising (theory-led) |
| 3 | **MAGUS with Clustal Omega backbones (proteins)** (sota2026 → bbtool-1..6, bbevidence, protbench) | Pilot: aligning only GCM's 10 backbones with Clustal Omega instead of MAFFT L-INS-i: −1.6 to −1.9 pts on BAliBASE (5/5 sets, up to −3.0), via SPFP; backbones much cheaper. **First confirmation runs are mixed**: BBA0117 draw 0 +1.32 paired / +0.29 end to end (worse), BBA0039 draw 0 −0.13 paired / −0.14 end to end (pilot on that set: −0.45). No effect on RNASim | open (confirmation running: 8 sets × 3 draws) | expected faster on proteins (MAGUS is slow there: BBA0190 41.6 min vs PASTA 12.6); measured end-to-end timings pending | NOVEL (no published non-MAFFT GCM backbones); mechanism claim partly known (M-Coffee: consensus needs decorrelated errors) | **high**: instructor's project list says verbatim "study MAGUS with different base methods" | lead, unconfirmed. Conflicts with "keep MAFFT a black box" unless a MAFFT-only recipe matches it (being tested) |
| 4 | **Soft-constraint MAGUS** (main track) | Measured end to end on 10 ROSE datasets: self-soft −0.77 pts vs MAGUS (9/1/0, p = 0.004) at +9% wall-clock; fan-out (54 reps): self-soft −0.55 (47/4/2, p = 7e-9), slow-soft-m3 −0.67 (47/1/5); **worse on BAliBASE** (self-soft 22.1 vs 20.6); no gain in ML trees (p = 0.90) | yes on DNA (≈ 29% of MAGUS's own gain over PASTA, −2.67); no on proteins | +9% vs MAGUS | partly known; must add MAGUS `-c false` baseline; not "first to relax constraints" | **medium-high**: MAGUS taught; "improve MAGUS / merging" on project list | modest |
| 5 | **GTM blending** (gtm) | GTM-Blend-ML vs GTM −14.3 / −6.1 / −10.1 FN pts (20/0/0, p = 9e-5) in simulation | in simulation only; loses to full IQ-TREE (+11); n.s. on published data | moderate | partly known; core (likelihood-decided blending) novel; open problem still listed in 2025 | **high**: D&C deck verbatim "develop a better DTM that allows blending" | promising if framed "blending fixes GTM when the guide tree is poor" |
| 6 | **Sample complexity of ASTRID/NJst vs ASTRAL** (samplecx) | caterpillar n = 64, f = 0.1: ASTRID needs 1189 genes vs ASTRAL 492; constant grows 4.9→15.3 vs 4.0→7.6 (n = 8→64); CIs overlap | n/a | cheap | theory open; empirical partly known; not "Roch's conjecture" | **high**: slide 35 verbatim | promising as evaluation, suggestive only |
| 7 | **Which alignment errors matter for trees** (alncrit) | at equal SPFN, real aligner errors cost 2–9× more tree accuracy than random errors (p = 0.020) | n/a | cheap | open | **high**: MSA deck verbatim "Which alignment criteria are predictive of tree accuracy?" | unclear, promising angle |
| 8 | **Long queries / second domain in UPP, WITCH, EMMA** (lenhet) | all backbone-anchored adders lose an extra shared region (long-query SPFN 0.49–0.50); trim → add → detect-and-align flanks fix: 6/6 wins over base method and over `mafft --add` (p = 0.031) | yes, but the scenario is partly true by construction | cheap | gap: no UPP-family paper benchmarks queries longer than the family | medium: UPP/WITCH taught; question ours | unclear, leaning promising as a narrow benchmark |
| 9 | **CAMUS: adaptive quartet filter / base tree** (camus) | base tree is the dominant error source (oracle −0.071, p = 1e-9) but every practical fix fails; sample-size-aware filter −0.014 held out (21/33/6, p = 0.013) | small | cheap | filter novel for CAMUS, anticipated by NANUQ/TINNiK | medium | unclear, leaning promising for a narrow project |
| 10 | **Faster MAGUS at equal accuracy** (fastmagus) | skeleton-100 FastTree guide tree + 4 backbones + one self-soft round: 1.82× faster, −0.16 pts vs our MAGUS (2/1/1, n = 4) | ties MAGUS at 1.8×; not ≥ 2× | 1.8× faster | not audited | medium | leaning not (but a useful speed knob) |
| 11 | Methods under new simulation models (models) | AliSim matched to ROSE statistics does not reproduce ROSE's difficulty; clock violation raises the tree cost of alignment error | n/a; run-to-run variance of PASTA/MAGUS swamps factors at 3 reps | ~60 machine-h for a full design | not audited | medium | unclear, leaning not (narrower clock-violation version promising) |
| 12 | Fragment-aware ML heuristic (ml) | matches RAxML-NG tree error on fragmentary data at ~1/3 the CPU (n = 5); alignment-confidence masking negative | ties at lower cost | cheaper | known practice (benchmark only) | medium | small effect; secondary experiment at best |
| 13 | Rogue taxa (rogue) | MAGUS/PASTA change < 0.1 pt; with MAFFT FFT-NS-2, injected rogues cost 0.8 FN pts (p = 0.03, n = 7), removed by an HMM filter | no for MAGUS/PASTA | cheap | not audited | medium | leaning not |
| 14 | DISCO-R (disco) | vs ASTRID-DISCO −0.002 RF, p = 0.32 | no | cheap | not found | high (slide 35) | not promising |
| 15 | Consensus alignment (main track) | −0.78 vs chosen input, ties MAGUS(Slow), worse than best input | no | slow | partly known | medium | modest |
| – | DL subset trees + GTM (dldtm) | pretrained DL subset trees much worse than cheap classical ones; GTM cannot repair | no | – | – | medium | not promising |
| – | Pairwise merging in PASTA (pairmerge) | exact DP vs OPAL n.s. | no | – | algorithm published | medium | not promising |
| – | Quartet amalgamation (quartets) | ASTRAL-IV already reaches the best quartet score in 377/380 | no | – | (a), (b) published | high | not promising |
| – | Supertrees at scale (supertree) | ties ASTRAL-III at ~1/30 memory; SCS dominates | no | – | low | medium | not promising |
| – | Forest+DTM (forest) | gains come from the decomposition; FastTree beats all | no | – | new but not useful | medium | not promising |
| – | Linguistic models (ling) | new method loses; idea published (Canby 2024) | no | – | low | medium | leaning not |
| – | Learned GCM edge weights (local) | −0.08 held out on simulated data, worse on proteins | no | free | new | low | not promising |
| – | MAGUS-lite / small backbones (local) | +1.5–2 pts or worse | no | faster | – | medium | not promising |

## Measured end-to-end runtime and accuracy (same 4-core machine type, everything run from unaligned sequences)

Ten ROSE 1000-sequence datasets (R0 of every condition):

| method | mean wall-clock | mean error (SPFN+SPFP)/2 | paired vs MAGUS |
|---|---|---|---|
| PASTA 1.8.3 (the paper's version) | 25.1 min | 9.32% | MAGUS vs PASTA −2.67 pts (10/0/0, p = 0.002) at 0.88× the time |
| MAGUS | 21.9 min | 6.64% | – |
| MAGUS(Slow) | 23.9 min | 6.39% | −0.26 (6/0/4, p = 0.38), 1.09× |
| self-soft MAGUS | 23.8 min | 5.87% | −0.77 (9/1/0, p = 0.004), 1.09× |
| slow-soft MAGUS | 24.8 min | 5.92% | −0.72 (8/1/1, p = 0.014), 1.13× |

Other datasets: RNASim 1000 R0: PASTA (rerun with `-d rna`) 72.2 min / 10.10%, MAGUS 9.80%, Slow 9.05%, self-soft 9.42%, slow-soft 9.13% (the first PASTA run failed because it used `-d dna` on RNA data;
wall-clock 87–92 min; a second machine measured 71 min for the same MAGUS run, so this is real: the end-to-end benchmark runs MAGUS's original pure-Python graph builder (`--gcmx-fastgraph false`), which is slow on RNASim's long alignments; the ~15 min seen elsewhere used our vectorized builder, which builds the identical graph).
BAliBASE BBA0101 / BBA0190: PASTA 5.9 / 12.6 min at 29.68 / 24.16%; MAGUS 14.9 / 41.6 min at 27.98 / 23.22%;
self-soft 21.0 / 44.0 min at 27.13 / 23.27%. On proteins MAGUS is 2.5–3.3× slower than PASTA, which is
where cheaper backbones (row 3) would matter. 16S.M R0: PASTA (rerun with `-d dna`, `claude/cs581-pastafix`) 18.4 min / 14.06% (the first run, 13.05%, had used `-d protein` because of a type-check bug), MAGUS 24.5 / 13.01%,
Slow 23.3 / 13.19%, self-soft 26.5 / 12.98%, slow-soft 27.2 / 13.12% (all within 0.2 points).

**Reproduction check (our runs vs the paper's published alignments of the same replicate, rescored with
the same FastSP):** median difference +0.02 points for PASTA (n = 11; mean +0.61, driven by 1000L1 +3.8 and
1000M2 +2.3), −0.02 for MAGUS (n = 12; mean +0.05, range −1.5 to +1.5) and +0.02 for MAGUS(Slow) (n = 12;
mean +0.37). BAliBASE: PASTA +0.51 / +0.25, MAGUS −0.47 / +0.10 (BBA0101 / BBA0190). Single runs differ
from the published ones by up to ±1.5 points because MAGUS and PASTA are unseeded; on average the harness
reproduces the paper. PASTA reruns with the correct datatype: RNASim 10.10% vs published 10.08%; 16S.M 14.06%
vs 12.99% (+1.07, single unseeded run).

On 4 cores MAGUS is only ~12% faster than PASTA on 1000 sequences (the paper's 2.5× used 16-core nodes).
Caveat found tonight: the MCL-threading patch used by the soft merges could silently fall back to one
thread (MAGUS runs the JSON copy of a task); fixed in `run_magus.py`. Same clusters, so accuracy is
unaffected and soft-MAGUS times above are, if anything, overstated.

## Running overnight (results will update rows 3 and 4)

- `claude/cs581-bbtool-1..6`: Clustal vs L-INS-i backbones on all 8 BAliBASE RV100 sets × 3 fresh MAGUS
  draws (paired merge-only + measured end to end), MAFFT-only controls (`--auto`, L-INS-i without `--ep`,
  G-INS-i), mixed 10+10 backbones, and nucleotide controls (RNASim, 1000M2, 16S.M, 1000L1).
- `claude/cs581-bbevidence`: mechanism (false cross-subset edges, error decorrelation), when it helps vs
  hurts, more/cheaper/mixed backbones, MAFFT-only recipes (gap penalties, `--unalignlevel`, column masking).
- `claude/cs581-protbench`: does it generalize (10AA, HomFam, simulated proteins with tree accuracy).
- A literature scout for ideas beyond the slides; pilots for the best of those will follow.

**gcmvote fan-out and first held-out rows (02:55 UTC, Oct 11).**
- Split across 26 more machines (one dataset/replicate each); see SESSIONS.md. 8 of the 14 new tree helpers stalled on a self-matching
  `pgrep -f` wait and were restarted at ~02:45.
- Training, proteins (main, n = 5, Δ SP error vs magus): es4 −1.20 (3/1/1, p = 0.19), vote hard −1.73 (4/1/0, p = 0.06), soft −0.33,
  soft2 −0.56. SIMHIGH_R1: magus 23.48, es4 −5.16, hard −5.48, soft +1.01.
- Held-out so far (SP error %):
  - SIMMOD_R2 (h5): MAGUS 11.12; es4 6.02, es5 5.63; vote hard 6.24, hard-bb 6.07, soft 11.29, soft4 9.56.
  - SIMHIGH_R2 (h5): MAGUS 24.03; es4 14.86, es5 14.54; vote hard 16.75, hard-bb 16.69, soft 25.39, soft4 19.68.
  - 1000L3 (h4, DNA): MAGUS 11.86; es4 11.55; vote hard 12.38.
  - 1000M1 R0 (h10, DNA): MAGUS 9.81; es4 9.70; vote hard 10.40, hard-bb 10.40, soft 10.40 (three distinct alignments).
  - BBA0154 (h3): MAGUS 21.23; es3 21.16, es4 20.86, es5 20.43 (vote rows pending).
- Early read: soft weighting does not work; the hard vote cutoff is close to es4 on proteins but worse on DNA (it keeps too
  little evidence there?). Checks: vote `magus` reproduces MAGUS exactly (h2 R9, h5, h10).

**gcmtrees FINAL (02:50 UTC; `claude/cs581-gcmtrees`, 16/16 datasets).**
- Recipe vs MAGUS, FastTree nRF: −0.19 (8/2/6, p = 0.63) although SP error drops on 16/16 (−5.7; ΔSPFN −10.8, ΔSPFP −0.6).
  SIMHIGH −0.49 (p = 0.46; true alignment −3.51, p = 0.008, so ~14% of the room). IQ-TREE SIMHIGH R1–R3 −0.47 (p = 0.25).
- Why: noise (needs ~50 SIMHIGH sets to resolve −0.5 RF; ΔRF tracks ΔSPFN only weakly, ρ = 0.32); recovered pairs are in
  informative columns; oracle split shows MAGUS's tree loss is mostly over-splitting, and the recipe's remaining wrong pairs cost
  more RF (+1.85 vs +1.10) than MAGUS's, so gains and losses nearly cancel. Route to tree gains: FP control on merged clusters + power.
- basemeth (02:51 interim, 14 protein sets): no aligner beats L-INS-i; FAMSA worse (Holm p = 0.003). Final ~04:15 UTC.
- HomFam (h3d, real proteins, Homstrad-seed scoring, SP error %): blmb MAGUS 22.81; es4 21.45; vote hard 21.59, hard-bb 20.84,
  soft 23.03. aat MAGUS 27.54; es4 27.83 (es3 28.18); vote hard 25.51, hard-bb 26.42, soft 24.65, soft4 24.42. First sign the
  vote model can beat es4 on real proteins (aat −2.0 vs +0.3); vote `magus` reproduces gg `linsi` on both.
