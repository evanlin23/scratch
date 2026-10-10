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
