# Pre-registration: better reference-free vote models for GCM edge support (gcmvote2)

AI-assisted (Claude), exploration code for CS581 project.

Written and pushed **before any variant alignment was scored** (2026-10-11). Seen before writing:
- the `magus` control merge on BBA0039 (avgErr 0.046892) and 1000M2 (0.082250): both reproduce the cached MAGUS
  runs to the last digit; the v1 `hard-bb` re-implementation keeps exactly the edges vote.py keeps (BBA0039: 114,847);
- fit-only summaries (kept edges, fitted parameters; **no reference, no accuracy**) of every model on BBA0039 and
  1000M2 (training). They led to two definition changes before this file: (i) a continuous-k beta-binomial is not a
  proper likelihood (positive log-likelihood, concentration at its bound), so M1-gbb rounds the graded support to an
  integer; (ii) M3 as "beta-binomial on (k, n) rescaled to round(n_eff)" collapsed 1000M2 (n' ≈ 2, every edge k' = 1,
  76 % of the weight dropped), so M3 is defined as an overlap-tempered posterior instead (below);
- from the brief / gcmvote: v1 `hard-bb` ties es4 on held-out proteins, is worse on DNA/RNA, and fails on BBA0081
  (195 sequences, so every 200-sequence backbone holds all sequences: near-identical backbones). The M3 guard was
  written with that failure in mind; BBA0081's features have not been looked at.

## Evidence (code/vote2.py)

Per cross-subset GCM edge e = (column a of subset A, column b of subset B) and backbone j:
x_j = residue pairs backbone j aligns between a and b (MAGUS's weight w = Σ_j x_j); r_j(a), r_j(b) = residues of
nodes a, b in backbone j; exposed_j ⇔ r_j(a) > 0 and r_j(b) > 0; n = #exposed; k = #{j : x_j > 0};
**s_j = x_j / (r_j(a)·r_j(b))**, the share of the residue pairs backbone j could align between a and b that it does
align (in [0, 1]). (Deviation from the brief's pairs / min(res): that ratio is not bounded by 1 — 8 × 8 fully aligned
residues give 64 / 8 = 8 — so the product is the "possible matches".) Graded categories: none (exposed, s = 0),
weak (0 < s ≤ 0.5), strong (s > 0.5). Only edges with k ≥ 1 exist, so every model is zero-truncated.

## Variants (merge-only on MAGUS's own subsets and 10 backbones; `-f 4`, MCL, minclusters, as gcmvote)

hard = keep the edge (MAGUS weight) iff posterior P(true) > 0.5; soft = weight × posterior.

| name | model |
|---|---|
| `magus` | control (must reproduce MAGUS) |
| `es4` | delete edges with k < 4 (baseline) |
| `hard-bb` | v1: zero-truncated 2-component beta-binomial mixture on (k, n) (baseline) |
| `m1-gbb` | **M1** beta-binomial mixture on graded support k_s = round(Σ_j s_j) clipped to [1, n] |
| `m1-dm`, `m1-dm-soft` | **M1** zero-truncated 2-component Dirichlet-multinomial mixture on (#none, #weak, #strong) of n |
| `m2-ds`, `m2-ds-soft` | **M2** binary Dawid–Skene: per-backbone P(vote \| true, exposed), P(vote \| false, exposed); EM with the unobserved no-vote edges as missing data |
| `m2-dsg` | **M2** graded Dawid–Skene: per-backbone categorical (none / weak / strong) per class |
| `m3-ovbb` | **M3** hard-bb fit; per-edge log-likelihood ratio × t, t = (n_eff/n) / median(n_eff/n) ∈ [0, 2]; n_eff = n² / Σ_{i,j ∈ exposed} J_ij, J_ij = Jaccard overlap of the A∪B sequences held by backbones i, j. **Guard:** if the edge-median n_eff < 2, the votes are not independent → return `magus` |
| `m23-ovds` | **M2+M3** binary Dawid–Skene, log-likelihood ratio × n_eff/n (DS assumes independent votes); same guard |
| `m4-gmm`, `m4-gmm-soft` | **M4** 2-component full-covariance Gaussian mixture on (logit (k+.5)/(n+1), logit mean s over voting backbones, log w); fitted on 300k random edges (seed 0, 3 inits); true = higher mean first coordinate |

Candidates for selection: the 10 M-variants above (`hard-bb` and `es4` are baselines, not candidates).

## Data and split

- **Training (selection only):** the gcmvote/gcmgen split — BBA0101, BBA0134, BBA0067, BBA0039, 1000M2, 1000L1,
  1000L2, 16S.M (cached MAGUS draws, `cs581/experiments/runs/*_R0/inputs.tar.xz`, the draws gcmvote trained on),
  SIMMOD_R1, SIMHIGH_R1 (bank h6b, backbones 1–10 = MAGUS's 10).
- **Held-out, primary** (rep bank, the helpers' MAGUS draws that check-in 12 used):
  - proteins (n = 9): SIMMOD_R2, SIMHIGH_R2, BBA0154, BBA0190, BBA0117, HomFam blmb, aat, Acetyltransf, PDZ;
  - DNA/RNA (n = 15): 1000L3, 1000M3, 1000S1, 1000S2, 1000M4, 1000S3, 1000M2_R1, 1000L1_R1, RNASim R0, RNASim R1,
    1000M1 R0–R3, 16S.3;
  - **BBA0081** reported on its own (and in a "proteins incl. BBA0081" line).
- **Held-out, secondary:** SIMHIGH bank R3–R12, R14–R20 (17 reps, true trees) with `magus`, `es4`, `hard-bb` and the
  selected variant (SP error and trees). Every held-out rep is scored with all variants where time allows, but
  held-out rows are never used to choose.

## Selection (training only)

Score = mean Δ error over the 10 training reps (Δ = variant − `magus`, error = (SPFN + SPFP)/2 × 100).
Eligible: candidates whose training DNA/RNA mean Δ ≤ +0.10 (v1's failure mode). Pick the eligible candidate with the
lowest score; ties within 0.05 → the smaller worst-case Δ. If none is eligible, the lowest score overall.
Selection is recorded in this file (appended) before any held-out variant row is looked at.

## Endpoints (held-out)

- **Primary:** selected variant vs `es4`, separately for proteins (n = 9) and DNA/RNA (n = 15), and pooled (n = 24):
  mean and median Δ, W/T/L (tie band ±0.1), two-sided Wilcoxon signed-rank (scipy default). Same vs `magus` and vs
  `hard-bb`.
- BBA0081: Δ vs `magus`, `es4`, `hard-bb`.
- Secondary: every candidate on the primary held-out reps (exploratory, labelled as such); SIMHIGH bank rows.
- Trees (if time allows): FastTree 2 `-lg -gamma` nRF vs the true tree on the SIMHIGH bank reps for true alignment,
  `magus`, `es4`, `hard-bb`, selected; Δ nRF vs `magus` and vs `es4`, paired Wilcoxon.

## Verdict rules (fixed now)

- **"Better than es4"**: pooled mean Δ vs es4 < 0 with Wilcoxon p < 0.05, and neither data type's mean Δ vs es4 > +0.10.
- **"Principled replacement for es4"** (gcmvote's rule): pooled mean Δ vs es4 ≤ +0.10 and no data type worse than es4
  by > 0.25.
- **"Fixes BBA0081"**: selected − `magus` ≤ +0.5 on BBA0081 (v1: +4.22).
- Otherwise the selected model is reported as not an improvement over es4.

## Addendum 1 (04:45 UTC, orchestrator request; written before any variant score was looked at)

The orchestrator's M5/M6 request (sent 04:17) arrived after this file was pushed (04:16) and after training merges
had started (04:17). No training or held-out variant score had been read when this addendum was written (only the
`magus` controls quoted above).

**M5, dataset-level gate (reference-free).** Per replicate, gate = filter with the selected M-variant if a signal is
on the filtering side of a threshold θ, else `magus`. Candidate signals, all from the 10 backbones' evidence, no
reference: (1) π of the hard-bb fit; (2) p1 − p0 of the hard-bb fit; (3) share of MAGUS's cross-edge weight on edges
with k < 4; (4) share of weight the selected variant removes; (5) median n_eff; baselines (6) **AOS** (MUMSA-style
average overlap score: mean over backbone pairs i < j of Σ_e min(s_i, s_j) / Σ_e max(s_i, s_j) over edges exposed in
both) and (7) **dispersion** (Muscle5-style: weight-weighted mean over edges of the s_j variance across exposed
backbones). For each signal and direction, θ = midpoint between adjacent training values maximising the training mean
Δ vs magus. The gate's signal is chosen by **leave-one-out CV** over the 10 training reps (θ and direction refitted per
fold; score = mean held-out-fold Δ); ties within 0.05 → AOS, then dispersion (the simple baselines win ties). The
gate counts as useful only if its LOO-CV score beats or matches both AOS and dispersion and the ungated selected
variant. Strictness (a third arm) is not fitted: 10 training points are too few for two thresholds. Held-out
endpoint: gated vs es4 and vs magus, as the primary endpoints above (secondary to the pre-registered selected model).

**M6, clustering combo (secondary).** Selected variant + Leiden-CPM (γ = 0.02, gcmclust's implementation) vs the same
graph + MCL, and es4 + CPM, on the primary held-out reps as time allows; Δ, W/T/L, Wilcoxon per data type.

**Extra held-out nucleotide set:** 16S.T (bank h12) joins the DNA/RNA held-out list if it can be fetched (n = 16).

Framing (novelty audit): WITCH already uses a weighted GCM (HMM-probability weights per input alignment); M1–M4 are
statistically defined edge / backbone weights. Support thresholding is known (Lassmann & Sonnhammer 2007; GUIDANCE;
TCS); M-Coffee found aligner weighting does not beat no weighting, so for M2 the spread of fitted backbone
reliabilities is reported and a near-uniform result is a finding.
