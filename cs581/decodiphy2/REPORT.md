# Identifiability in phylogenetic distance deconvolution — round 2 (theory + identifiability-aware output)

*CS581 project exploration, Oct 2026. Builds on `../decodiphy/REPORT.md` (round 1). Full proofs:
`theory.md`. Code: `code/`. Tables: `results/`. Source paper: Arasti, Şapcı, Rachtman, El-Kebir &
Mirarab, "Deconvolving Phylogenetic Distance Mixtures", RECOMB 2026 (bioRxiv 10.64898/2026.01.18.700179).*

## Headline numbers

* **Proofs:** k = 1, k = 2 (both non-adjacent and adjacent placements), and k = 3 are now **fully
  proved**. For k = 3 non-adjacent placements: non-identifiable ⇔ continuum ⇔ closed claw.
* **General k:** we stated a Hall-type criterion. Non-adjacent placements S are non-identifiable iff
  some set U of internal nodes has |N[U] ∖ V(S)| < |U|. One direction is proved. For the other:
  * the continuum ⇔ Hall step held in all 145,543 (shape, S) pairs checked (n ≤ 10, k ≤ 5, 4 length
    draws);
  * alternative ⇔ continuum held in **≈ 140k MILP-certified continuum-free matchings with k = 4–5,
    n ≤ 12, with 0 counterexamples**.
* **Corollary:** more than n/2 non-adjacent placements are never identifiable.
* **Supplement checked** (SB.3, bioRxiv media-1.pdf). The paper proves the adjacency continuum
  (Claim 3), boundary ambiguity x ∈ {0, 1} and ȳ-only identifiability (Claim 2), and k = 2 *edge*
  identifiability only for p = (½, ½) (Claim 5). Our k = 2 theorem (any p, x, ȳ, alternatives with
  k' ≤ 2, and the adjacent case) and everything about claws, Hall deficiency and k = 3 are new.
  The paper's Fig. 1d counterexample uses **adjacent** placements. With random lengths the same
  configuration still has exact alternatives, so symmetry is not the cause.
* **On the authors' own stored E1 runs** (4,171 runs, 12 trees, k ∈ {2, 3, 5, 7, 10}):
  * 12.1% of the *true* placement sets are adjacent and 4.0% have a Laplacian kernel (claw-type).
  * 4.7% of DecoDiPhy's *final fits* have a non-trivial equivalence class (11.4% noise-free; 25% for
    noise-free k = 10), and 8.9% carry at least one flag (also counting x at 0 or 1).
  * Where the class is non-trivial, the abundance range inside it has median 0.10 (90th percentile
    0.23).
* **This matters for accuracy.** On runs where DecoDiPhy found exactly the true edges, the abundance L1
  error is **0.0003 when the class is trivial vs 0.119 when it is not** (noise-free). Flagged runs
  carry **97%** of all abundance error noise-free and 30–38% with noise. The class contains the truth
  (truth inside the reported intervals in 86% of flagged runs, vs 0% for DecoDiPhy's point estimate).
  No point estimate can repair this: a central "canonical" representative is not significantly
  better (p ≥ 0.2). **The fix is to report the class (intervals plus a flag), not a better point.**

## 1. Theory (details and proofs in `theory.md`)

Round 1 showed d = F m + ȳ 1, where m is the node measure, and that d determines (m, ȳ) up to
Laplacian moves g = −Lα at internal nodes. So an alternative edge set S' exists for some parameters
iff a sign pattern on g is feasible: g < 0 only on V(S), g > 0 only on V(S'), and 0 elsewhere.
Round 2 adds a *private-neighbour lemma*: each leaf of the subtree spanned by supp α pushes a same-sign,
non-adjacent sibling pair into V(S) or V(S'). Combined with a covering count (pairwise non-adjacent
nodes need one edge each) and short "far partner" arguments, this gives:

| result | status |
|---|---|
| k = 1 identifiable | proved |
| k = 2, non-adjacent: (S, p, x, ȳ) identifiable vs all k' ≤ 2 (round 1's open sub-case: two T_U-leaves of opposite sign, closed by a "far partner" argument) | **proved** |
| k = 2, adjacent: S, m, ȳ identifiable; (p, x) a 1-dim continuum for **every** parameter value (explicit fibre) | **proved** |
| k = 3, non-adjacent: (A) alternative for some parameters ⇔ (B) continuum ⇔ (C) closed claw | **proved** (four-case analysis) |
| closed claw ⇒ different edge set on an open parameter set (generic, not symmetric) | proved (round 1) |
| continuum (B) ⇔ Hall deficiency (H) | ⇐ proved; ⇒ verified on 145,543 cases |
| k > n/2 non-adjacent ⇒ never identifiable | proved |
| general k non-adjacent: (A) ⇔ (B) ⇔ (H) | **conjecture**; 0 counterexamples, MILP over all S' (n ≤ 12, k = 4–5; k = 6 only occurs with k > n/2 up to n = 11) |
| adjacent placements, k ≥ 3: edge sets change *without* a continuum | observed (e.g. 3,213 / 4,729 at n = 9, k = 4); mechanism given |
| per-query pendant lengths y_q | everything carries over verbatim (d depends on y only through ȳ = pᵀy); y_q never identifiable for k ≥ 2 |

Figures: `figs/fig_claw.png` (k = 3 closed claw, n = 5, three exact solutions on different edge sets),
`figs/fig_adjacent.png` (the paper's Fig. 1d configuration with random lengths, still non-identifiable),
`figs/fig_splitclaw.png` (k = 4 split claw: no closed claw, still a continuum). All were checked
numerically (`results/figs_verification.txt`, |d − d'| ≤ 1e−15 for the constructed ones).

MILP evidence for the conjecture (`results/milp_enum_*.md`; random lengths, unit-length rerun identical):

| n | k | matchings | with continuum | continuum-free (each MILP-checked against all S', \|S'\| ≤ k) | counterexamples |
|---|---|---|---|---|---|
| 8 | 4 | 400 | 292 | 108 | 0 |
| 9 | 4 | 1,652 | 844 | 808 | 0 |
| 10 | 4 | 6,752 | 2,556 | 4,196 | 0 |
| 11 | 4 | 21,544 | 6,148 | 15,396 | 0 |
| 12 | 4 | 78,420 | 17,636 | 60,784 | 0 |
| 10 | 5 | 5,660 | 4,880 | 780 | 0 |
| 11 | 5 | 25,120 | 17,604 | 7,516 | 0 |
| 12 | 5 | 118,116 | 67,572 | 50,544 | 0 |
| 8–11 | 5–6 with k > n/2 | 17,988 | all (Corollary 6) | 0 | — |

## 2. Overlap with the paper's supplement (question 2)

Retrieved the main PDF from PMC and the supplement (`media-1.pdf`, 8 pp.) from bioRxiv with headless
Chromium. Contents of SB.3:

* **SB.3.1 / Claim 1:** two queries on one edge merge into one.
* **SB.3.2 / Claim 2:** x ∈ {0, 1} is unidentifiable (a 4-leaf example), and **SB.3.3:** only ȳ is
  identifiable, not the individual y_q. Our §7 per-query result is this claim restated in the measure
  framework.
* **SB.3.4 / Claim 3:** adjacent edges give infinitely many (p, x), by counting independent leaf
  equations and ignoring the inequality constraints. This is our C1 / Theorem 3(2). Ours adds that the
  fibre is explicit and non-empty for every parameter value, and that S, m and ȳ remain identifiable
  at k = 2.
* **SB.3.6 / Claim 5:** for k = 2 **and p = (½, ½)** the placement *edges* are unique (four-point
  condition, 3 configuration cases). Theorem 3 strictly generalizes it.
* Nothing on claws, Laplacian moves, k ≥ 3 characterization, Hall conditions, or how often these occur.
  Main text, verbatim: "we suspect (with no proof) that conditions that break identifiability require
  extreme cases of symmetry … We leave a full characterization of identifiability to future work."
* SB.1.1 (Algorithm S1) post-processes fits with x ≈ 0 or 1 by trying neighbouring edges. That is the
  only identifiability handling in their code, besides `check_identifiability` (flags placements
  within 2 edges, used only in `--score` mode).

**Novelty verdict:** the k = 2 generalization, the k = 3 theorem, the Hall criterion and conjecture,
the k > n/2 corollary, and the refutation of the symmetry conjecture are not in the paper or its
supplement. Round 1's literature search found no follow-up work.

## 3. Identifiability-aware output (question 3)

`code/idclass.py` post-processes any PDD fit (DecoDiPhy's `all_rounds.json` + labeled tree; CLI
`code/decodiphy_idcheck.py`):

1. Converts the rooted tree to unrooted form (suppresses the degree-2 root) and computes the node measure.
2. Raises three flags: **boundary** (x within DecoDiPhy's own 1e−4 / 1e−3 thresholds of 0 or 1),
   **adjacent**, and **kernel** (rank test on L[V ∖ V(S), I], the continuum condition (B)).
3. Computes the **equivalence class on the fitted edges**: every (p', x', ȳ') with the same fitted d.
   This is a polytope in (α, per-node shares). 2k LPs (HiGHS) give the interval [lo_q, hi_q] for each
   abundance, and the mean of the 2k optimal vertices is a canonical central representative. Every
   member has exactly the same loss, so the solver's choice among them is arbitrary.
4. Runtime is negligible: 4,343 runs on 2 cores in ~25 min, dominated by birds-stiller (n = 360).

Example (`results/idcheck_example.txt`, 1kp, k = 10, noise-free): DecoDiPhy reports p = 0.0138 for one
placement. The class interval is [0.0029, 0.0154] and the truth is 0.005; the adjacent placement trades
mass with it.

**Measured on the authors' stored E1 runs** (`DecoDiPhy-Data/biotrees`, 4,343 run directories, 4,171
complete; `code/analyze_data.py`, `code/summarize_data.py` → `results/data_summary.md`):

| | noise-free | noise1 (placement noise, scale 1) | noise2 (scale 2) |
|---|---|---|---|
| true sets adjacent / with kernel | 12.1% / 4.0% | same truth | same truth |
| fits with non-trivial class | 11.4% (k = 10: 24.8%) | 3.2% (k = 10: 8.8%) | 3.9% (k = 10: 10.3%) |
| fits with any flag (incl. boundary x) | 12.6% | 9.3% | 7.3% |
| abundance L1 on exact-edge runs: class trivial | 0.0003 (n = 460) | 0.0044 (1,129) | 0.0040 (1,130) |
| abundance L1 on exact-edge runs: class non-trivial | **0.119** (n = 41) | **0.215** (10) | **0.126** (22) |
| share of all abundance error in flagged runs | **97%** | 30% | 38% |
| truth inside reported intervals (flagged, exact edges) | 88% | 90% | 82% |

* Within a flagged class the abundance range has median 0.10 and mean 0.125; 49% of flagged fits have
  a range > 0.10. On the tree EMD (weighted-UniFrac-type, pendants ignored), the spread across class
  vertices averages 0.0032, against a mean EMD error of 0.0044 for these runs.
* **Canonical vs DecoDiPhy's point:** EMD 0.00435 vs 0.00444 (111 better / 84 worse, Wilcoxon p = 0.2).
  Abundance L1 0.136 vs 0.134 (p = 0.77). Over the whole benchmark the EMD changes by −0.15%. **A
  different point estimate does not help.** The information is simply not in d.
* So the identifiability flag is a near-perfect *error predictor* for abundances on correctly placed
  fits, and interval output turns 0% exact coverage into 82–90% coverage. The remaining ~12% misses are
  within 0.002–0.018 L1 of the interval and come from the fits' residuals.
* Placement metrics: Jaccard is unaffected by the class on fixed edges. Edge-set alternatives arise only
  via boundary sliding (4.6% of fits have a boundary x, up to 13% at k = 10 with noise).

## 4. Question 4 (per-query pendant lengths)

d depends on the y_q only through ȳ = pᵀy, so the characterization is unchanged. Individual y_q are
never identifiable for k ≥ 2: the fibre is the polytope {y ≥ 0 : pᵀy = ȳ}. Constraints y_q ≥ 0 are
equivalent to ȳ ≥ 0, so they neither add nor remove ambiguity. Placement-dependent pendant heuristics
(DecoDiPhy's y_q ∝ clade height) would make the model nonlinear and are not covered (`theory.md` §7).

## 5. Verdict

**Verdict: promising.** As a 4-week CS581 project this is a crisp theory contribution with a measurable
practical consequence on the authors' own benchmark.

* Go/kill criteria: complete proofs for k ≤ 2 ✔; a proven k = 3 theorem ✔ plus a sharply stated,
  heavily tested general-k conjecture; an identifiability-aware output that matters on their data ✔.
  Flagged fits hold 97% of noise-free abundance error, intervals cover the truth 82–90% of the time
  vs 0% for points, and the flag fires on 11–25% of noise-free fits at k ≥ 5.
* Not "more accurate than the baseline" in point metrics: canonical representatives give no
  significant gain. Pitch it as **characterization + honest uncertainty output**, not accuracy.

**Weeks 1–4.**
1. Week 1: write up Theorems 2–4 cleanly (they are complete). Try to prove (B) ⇒ (H) via the
   all-minors matrix-tree theorem (tree Laplacian minors do not cancel). Start (A) ⇒ (B) for general
   k with the private-neighbour lemma (induction on |U| or on the number of T_U leaves). Fallback: k ≤ 3
   theorem plus conjecture with the MILP evidence pushed to n = 13–14.
2. Week 2: package `idclass` as a DecoDiPhy post-processor: flags, abundance intervals, and
   enumeration of alternative edge sets reached by boundary sliding. Optionally fold the Hall/adjacency
   test into the search as a tie-breaker or into the k-selection rule from round 1 (adjacency-stop).
3. Week 3: experiments: the E1 table above for all 12 trees and noise settings, plus a run of DecoDiPhy
   itself on fresh simulations, reporting interval coverage, interval width vs k and flag rates. If
   time allows, E2-style read-based d̂ on one WoL subtree.
4. Week 4: report: theory section, figures, tables; discussion of what remains open (general k,
   non-binary trees, noisy-d stability via the grounded Laplacian).

**Risks.** (i) The general-k (A) ⇔ (B) may resist proof in 4 weeks. The k ≤ 3 theorem plus evidence
is a safe fallback. (ii) The interval output is "only" uncertainty quantification. Some may see it
as incremental, and the strong numbers (97%, 0% → 88%) are noise-free; with noise the share is 30–38%.
(iii) Coverage is assessed on the authors' simulator, whose spacing rule (no query whose parent or
grandparent is taken) already *reduces* adjacency. Real data could have more or fewer flags. (iv)
Non-binary reference trees (polytomies) need the deg − 2 generalization, which the code supports but
the proofs do not cover.
