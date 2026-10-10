# Ranked comparison of explored project ideas (final, 2026-10-10 06:40 UTC)

Supersedes `ranking_v1.md`. Sources: every exploration session's `cs581/<dir>/REPORT.md` (branch
`claude/cs581-<dir>`), the novelty audits (`novelty_audits.md`), the MAGUS-variant fan-out (54 replicates,
`../experiments/variants/SUMMARY.md`) and the measured end-to-end benchmark (13 of 14 datasets,
`claude/cs581-e2e-*`). **Nothing is chosen here**; the overnight runs (bottom) may move rows 3 and 4.

"Course" = how easily the choice is justified from CS581 material: **high** = an open question printed
verbatim in the slides or on the instructor's project list, about methods taught in lecture; **medium** =
methods taught, question ours; **low** = needs material not covered.

| rank | idea (branch) | strongest result | beats strongest baseline? | runtime | novelty (audit) | course | verdict |
|---|---|---|---|---|---|---|---|
| 1 | **ASTRAL-Pro is inconsistent under GDL because of its own rooting/tagging** (gdlcons) | Exact 4-taxon formula; with true gene trees stock ASTRAL-Pro3 returns the wrong species tree in 40/40, 20/20, 4/4 datasets (500 / 2k / 10k families); correct with true tags. Holds **even with constant rates** (λ = μ on every branch), but only at extreme turnover (λ = μ = 8: wrong 20/20 and 4/4; never at λ ≤ 4). Tag-error threshold q* = 0.080 ± 0.005 observed vs 0.079 predicted | n/a (theory + simulation) | cheap (simulation) | core result novel: counterexamples to Zhang et al.'s consistency conjecture; partly known context | **high**: slide 35 asks verbatim "Is ASTRAL-Pro consistent for GDL under ... error for rooting and tagging?"; GDL lecture "Phylogenomics, part 2" | promising (theory-led). Caveats: numerical, not proofs; constant-rate failures need extreme rates; wQFM-GDL probably inherits the issue (unverified) |
| 2 | **ASTRID-Pro: a GDL-corrected internode distance** (gdl) | Theorem: ortholog-only internode distance is a tree metric with the species-tree topology iff no supercritical branch; counterexample matches to 3 decimals; held-out DISCO test on high-GDL data vs ASTRID-multi −0.0096 FN (16/12/5, Holm p = 0.021); ties ASTRID-DISCO / ASTRAL-Pro | only on high-GDL data; ties elsewhere | cheap | NOVEL (narrowly); do not claim ASTRID-multi is answered | **high**: slide 35 verbatim ("distance correction for GDL so ASTRID/NJst are consistent?"); ASTRID taught in ILS lectures | promising (theory-led) |
| 3 | **MAGUS with Clustal Omega backbones (proteins)** (sota2026 → bbtool-1..6, bbevidence, protbench) | Pilot: aligning only GCM's 10 backbones with Clustal Omega instead of MAFFT L-INS-i: −1.6 to −1.9 pts on BAliBASE (5/5 sets, up to −3.0), via SPFP; backbones much cheaper. **First confirmation run went the other way** (BBA0117 draw 0: +1.3 paired, +0.3 end to end). No effect on RNASim | open (confirmation running: 8 sets × 3 draws) | expected faster on proteins (MAGUS is slow there: BBA0190 41.6 min vs PASTA 12.6); measured end-to-end timings pending | NOVEL (no published non-MAFFT GCM backbones); mechanism claim partly known (M-Coffee: consensus needs decorrelated errors) | **high**: instructor's project list says verbatim "study MAGUS with different base methods" | lead, unconfirmed. Conflicts with "keep MAFFT a black box" unless a MAFFT-only recipe matches it (being tested) |
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

Other datasets: RNASim 1000 R0: MAGUS 9.80%, Slow 9.05%, self-soft 9.42%, slow-soft 9.13% (PASTA failed;
wall-clock 87–92 min on that machine, about 6× our other RNASim timings, so treat it as unreliable).
BAliBASE BBA0101 / BBA0190: PASTA 5.9 / 12.6 min at 29.68 / 24.16%; MAGUS 14.9 / 41.6 min at 27.98 / 23.22%;
self-soft 21.0 / 44.0 min at 27.13 / 23.27%. On proteins MAGUS is 2.5–3.3× slower than PASTA, which is
where cheaper backbones (row 3) would matter. 16S.M pending.

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
