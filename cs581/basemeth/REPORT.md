# MAGUS (and PASTA) with different subset aligners: pilot report

Projects addressed (instructor list, `cs581/notes/CS581-project-suggestions.txt`): (a) "study MAGUS with different
base methods", (b) "Prank or ProbCons within PASTA or MAGUS", (c) "Regressive vs PASTA and MAGUS". Pre-registration:
`PREREG.md` (committed before any result). Raw rows: `results/*.jsonl`; all tables: `results/tables.md`
(`code/summarize.py`). One ~6 h session on a 4-core VM.

## Verdict (numbers first)

15 protein datasets (5 BAliBASE RV100, 6 HomFam at 2,000 sequences, 4 AliSim simulations), one MAGUS draw each,
merge-only so that only the subset aligner changes. d = variant minus MAGUS's own L-INS-i subsets, in SP-error
points (negative = better).

| subset aligner | all 15 protein | BAliBASE (5) | HomFam (6) | AliSim (4) | Holm p (15 sets) | subset CPU vs L-INS-i |
|---|---|---|---|---|---|---|
| MUSCLE5 `-align` | +0.93 (5/0/10) | **-0.47 (4/0/1)** | +1.92 (0/0/6) | +1.19 (1/0/3) | 0.055 | 4-7x |
| ProbCons | +2.23 (4/1/10) | **-0.37 (4/1/0)** | +5.04 (0/0/6) | +1.27 (0/0/4) | 0.043 (worse) | 7-10x |
| Clustal Omega | +3.27 (4/0/11) | **-0.31 (4/0/1)** | +3.63 (0/0/6) | +7.22 (0/0/4) | 0.037 (worse) | 0.2-0.6x |
| FAMSA2 | +2.65 (2/0/13) | +0.26 (1/0/4) | +3.52 (1/0/5) | +4.33 (0/0/4) | 0.0017 (worse) | 0.03-0.3x |
| G-INS-i (not pre-registered) | +0.29 (8/2/5) | +0.05 (2/1/2) | +0.97 (3/0/3) | **-0.42 (3/1/0)** | (p = 0.42 raw) | ~1x |
| Kalign 3 (not pre-registered) | +3.19 (3/0/12) | +0.26 | +5.79 | +2.95 | (p = 0.03 raw) | 0.03-0.2x |

- **Pre-registered primary test fails.** The best candidate by mean (MUSCLE5) is +0.93 points *worse* than L-INS-i
  inside MAGUS (Holm p = 0.055); ProbCons, Clustal Omega and FAMSA2 are significantly worse (Holm p = 0.043,
  0.037, 0.0017). No newer base method makes MAGUS more accurate overall.
- **BAliBASE is the exception, and it is small.** MUSCLE5, ProbCons and Clustal Omega lower MAGUS's error on 4 of
  5 RV100 sets by 0.3-0.5 points, mostly through SPFP. But re-aligning the same subsets with L-INS-i itself moves
  MAGUS by 0.03-0.31 points on these sets (median 0.11 over all 15; the `orig` vs `linsi` control), so the
  BAliBASE gain is only 2-4x the run-to-run noise and n = 5 cannot reach p < 0.05 (smallest p = 0.0625).
- **It reverses on simulated proteins (true alignment known) and on HomFam.** On AliSim every pre-registered
  candidate is worse on all 4 replicates, except MUSCLE5 on SIMHIGH_R2 (-0.06); Clustal Omega costs ~7 points. On HomFam (seed-scored, noisy) all three "BAliBASE winners" lose on 6/6 families (mean +1.9 to +5.0 points). This
  mirrors the earlier Clustal-backbone result (bbtool / protbench): BAliBASE rewards conservative, low-SPFP
  aligners; indel-rich simulated data penalizes them.
- **Mechanism is simple:** GCM preserves the subset alignments, so MAGUS's final error tracks the subset
  aligner's own accuracy on the 40-80-sequence subsets (subset-level table below). L-INS-i (and G-INS-i) are the
  most accurate subset aligners on everything but BAliBASE.
- **Only G-INS-i is interesting** (not pre-registered, exploratory): never clearly worse on simulated data,
  -0.76 on both SIMHIGH replicates (and -0.27 on ROSE 1000M2). Its trees are mixed: FastTree nRF 6.8% vs 10.1%
  for L-INS-i on SIMHIGH_R1 but 8.5% vs 6.7% on SIMHIGH_R2. HomFam is mixed (3 wins, 3 losses, mean +0.97 driven by
  Acetyltransf +8, a 6-seed family).
- **Trees barely follow these SP differences**: over 4 simulated sets, MUSCLE5 / ProbCons subsets change FastTree
  nRF by -1.3 to +1.6 points vs L-INS-i in both directions (table below); no tree-accuracy case for any swap.
- **Cost.** MUSCLE5 and ProbCons make the subset stage 4-10x more expensive (e.g. ~1,100-1,600 vs ~165 CPU s
  per 1,000-sequence AliSim dataset); FAMSA, Kalign and Clustal make it nearly free but are the least accurate.
  The subset stage is a small part of MAGUS's total time either way (MAGUS's 10 L-INS-i backbones dominate).
- **Whole-dataset tools are worse than MAGUS on BAliBASE and AliSim**: FAMSA2 alone is 0.8-13 points worse there
  (but better than MAGUS on 2 of 6 HomFam families, PDZ and zf-CCHH); regressive T-Coffee (NJ tree + Clustal children; see caveat) 18.3% vs MAGUS 4.7% on BBA0039. See
  the baselines section for what finished.

**Verdict for a 4-week CS581 project: not promising** as "find a better base method for MAGUS's subsets" (the
pre-registered question has a clear negative answer on 15 datasets, and the only positive signal, on BAliBASE,
is within 2-4x noise and contradicted by simulations). PASTA with ProbCons/Prank (project b) was impractical here
(1 h cap hit on 3 of 4 runs; ProbCons on BBA0039 +0.08 worse than MAFFT at 3.2x the time), and regressive T-Coffee
(project c) in the configuration that ran is far worse than MAGUS. It is a perfectly good *negative-result*
project (straightforward, cheap, clear story: GCM inherits subset accuracy; BAliBASE vs simulation disagree);
G-INS-i subsets on high-divergence data are the one thread worth a week for a positive angle.

## Main table: MAGUS final SP error (%) by subset aligner

Merge-only: the same 25 MAGUS subsets and the same 10 L-INS-i backbones in every column; only the subset aligner
changes. Bold = best in the row. `linsi` = MAGUS's own subset command re-run (the control).

| dataset | linsi | ginsi | muscle5 | famsa | probcons | clustalo | kalign | prank |
|---|---|---|---|---|---|---|---|---|
| BBA0039 | 4.66 | 4.64 | 4.48 | 4.92 | 4.61 | **4.45** | 5.03 |  |
| BBA0067 | 26.40 | 26.20 | 25.37 | 26.53 | **25.33** | 25.97 | 27.27 |  |
| BBA0101 | 28.88 | 29.32 | **28.11** | 29.63 | 28.65 | 29.00 | 29.21 |  |
| BBA0154 | 21.76 | 21.65 | 21.86 | 21.56 | 21.72 | **21.16** | 21.31 |  |
| BBA0190 | 23.43 | 23.55 | 22.98 | 23.77 | 22.98 | **22.98** | 23.63 |  |
| HF_Acetyltransf | **26.03** | 34.07 | 27.90 | 40.67 | 44.79 | 36.75 | 41.41 |  |
| HF_PDZ | 16.61 | **15.39** | 17.79 | 17.53 | 18.05 | 16.84 | 26.17 |  |
| HF_aat | 15.07 | **14.23** | 15.46 | 15.67 | 15.56 | 16.41 | 17.49 |  |
| HF_rrm | 18.41 | 19.61 | 20.81 | 20.09 | 20.51 | 20.48 | **16.98** |  |
| HF_sdr | 23.35 | **21.93** | 26.48 | 26.73 | 27.55 | 27.50 | 35.08 |  |
| HF_zf-CCHH | 12.25 | 12.33 | 14.83 | 12.13 | 15.51 | 15.51 | **9.32** |  |
| SIMMOD_R1 | 12.79 | 12.83 | 13.82 | 16.41 | 13.73 | 20.38 | 14.99 |  |
| SIMMOD_R2 | 11.06 | **10.86** | 11.41 | 13.88 | 11.46 | 17.06 | 12.07 |  |
| SIMHIGH_R1 | 25.00 | **24.19** | 28.45 | 30.21 | 28.06 | 34.53 | 30.68 |  |
| SIMHIGH_R2 | 23.63 | **22.92** | 23.57 | 29.32 | 24.29 | 29.36 | 26.51 |  |
| 1000M2 | 8.24 | **7.97** |  | 23.01 |  | 13.26 | 11.52 |  |

## Paired tests vs L-INS-i subsets (SP-error points; negative = better than MAGUS's default)

| method | set | n | mean d | median d | W/T/L | p (Wilcoxon) | d SPFN | d SPFP | d TC |
|---|---|---|---|---|---|---|---|---|---|
| ginsi | all protein | 15 | +0.29 | -0.11 | 8/2/5 | 0.421 | +0.22 | +0.37 | +0.1 |
| ginsi | BAliBASE | 5 | +0.05 | -0.02 | 2/1/2 | 1 | +0.22 | -0.13 | +0.3 |
| ginsi | HomFam | 6 | +0.97 | -0.39 | 3/0/3 | 1 | +0.75 | +1.19 | -0.7 |
| ginsi | SIMMOD | 2 | -0.08 | -0.08 | 1/1/0 | 1 | +0.09 | -0.24 | +1.0 |
| ginsi | SIMHIGH | 2 | -0.76 | -0.76 | 2/0/0 | 0.5 | -1.27 | -0.25 | +1.3 |
| ginsi | nucleotide | 1 | -0.26 | -0.26 | 1/0/0 | nan | -0.25 | -0.27 | +0.7 |
| muscle5 | all protein | 15 | +0.93 | +0.39 | 5/0/10 | 0.0554 | +1.22 | +0.64 | -0.8 |
| muscle5 | BAliBASE | 5 | -0.47 | -0.45 | 4/0/1 | 0.125 | -0.25 | -0.68 | +1.4 |
| muscle5 | HomFam | 6 | +1.92 | +2.13 | 0/0/6 | 0.0312 | +2.58 | +1.26 | -0.4 |
| muscle5 | SIMMOD | 2 | +0.69 | +0.69 | 0/0/2 | 0.5 | +0.54 | +0.84 | -5.3 |
| muscle5 | SIMHIGH | 2 | +1.69 | +1.69 | 1/0/1 | 1 | +1.48 | +1.90 | -3.0 |
| famsa | all protein | 15 | +2.65 | +0.93 | 2/0/13 | 0.000427 | +2.64 | +2.66 | -4.1 |
| famsa | BAliBASE | 5 | +0.26 | +0.26 | 1/0/4 | 0.188 | +0.44 | +0.07 | -0.5 |
| famsa | HomFam | 6 | +3.52 | +1.30 | 1/0/5 | 0.0625 | +3.22 | +3.82 | -3.0 |
| famsa | SIMMOD | 2 | +3.22 | +3.22 | 0/0/2 | 0.5 | +3.66 | +2.77 | -10.6 |
| famsa | SIMHIGH | 2 | +5.45 | +5.45 | 0/0/2 | 0.5 | +5.39 | +5.51 | -9.8 |
| famsa | nucleotide | 1 | +14.78 | +14.78 | 0/0/1 | nan | +15.33 | +14.22 | -12.4 |
| probcons | all protein | 15 | +2.23 | +0.67 | 4/1/10 | 0.0215 | +2.53 | +1.94 | -1.7 |
| probcons | BAliBASE | 5 | -0.37 | -0.23 | 4/1/0 | 0.0625 | -0.01 | -0.72 | +0.8 |
| probcons | HomFam | 6 | +5.04 | +2.68 | 0/0/6 | 0.0312 | +5.55 | +4.53 | -0.6 |
| probcons | SIMMOD | 2 | +0.67 | +0.67 | 0/0/2 | 0.5 | +0.46 | +0.88 | -7.6 |
| probcons | SIMHIGH | 2 | +1.86 | +1.86 | 0/0/2 | 0.5 | +1.86 | +1.87 | -4.9 |
| clustalo | all protein | 15 | +3.27 | +2.07 | 4/0/11 | 0.0125 | +3.13 | +3.41 | -4.6 |
| clustalo | BAliBASE | 5 | -0.31 | -0.43 | 4/0/1 | 0.125 | -0.15 | -0.48 | +1.3 |
| clustalo | HomFam | 6 | +3.63 | +2.66 | 0/0/6 | 0.0312 | +3.36 | +3.89 | -1.9 |
| clustalo | SIMMOD | 2 | +6.80 | +6.80 | 0/0/2 | 0.5 | +7.20 | +6.40 | -16.0 |
| clustalo | SIMHIGH | 2 | +7.63 | +7.63 | 0/0/2 | 0.5 | +6.55 | +8.71 | -15.9 |
| clustalo | nucleotide | 1 | +5.03 | +5.03 | 0/0/1 | nan | +5.38 | +4.67 | -13.5 |
| kalign | all protein | 15 | +3.19 | +1.01 | 3/0/12 | 0.0302 | +6.61 | -0.23 | -5.5 |
| kalign | BAliBASE | 5 | +0.26 | +0.33 | 1/0/4 | 0.438 | +3.39 | -2.86 | +0.6 |
| kalign | HomFam | 6 | +5.79 | +5.99 | 2/0/4 | 0.219 | +11.09 | +0.48 | -11.6 |
| kalign | SIMMOD | 2 | +1.61 | +1.61 | 0/0/2 | 0.5 | +2.23 | +0.98 | -1.9 |
| kalign | SIMHIGH | 2 | +4.28 | +4.28 | 0/0/2 | 0.5 | +5.59 | +2.97 | -6.1 |
| kalign | nucleotide | 1 | +3.28 | +3.28 | 0/0/1 | nan | +7.32 | -0.76 | +3.3 |

## Pre-registered primary test

| candidate | n | mean d vs L-INS-i | W/T/L | raw p | Holm p |
|---|---|---|---|---|---|
| muscle5 | 15 | +0.93 | 5/0/10 | 0.0554 | 0.0554 |
| probcons | 15 | +2.23 | 4/1/10 | 0.0215 | 0.0431 |
| famsa | 15 | +2.65 | 2/0/13 | 0.000427 | 0.00171 |
| clustalo | 15 | +3.27 | 4/0/11 | 0.0125 | 0.0374 |

## Why: subset accuracy carries through to the final alignment

Subset-level SP error (%), the 25 subset alignments of each method scored against the reference restricted to the subset, pooled:

| set | linsi | ginsi | muscle5 | famsa | probcons | clustalo | kalign | prank |
|---|---|---|---|---|---|---|---|---|
| BAliBASE (n=5) | 10.43 | 10.48 | 10.14 | 10.77 | 10.21 | 10.06 | 10.91 |  |
| HomFam (n=6) | 8.86 | 9.40 | 12.30 | 11.19 | 17.60 | 11.84 | 14.17 |  |
| SIMMOD (n=2) | 3.35 | 3.12 | 4.26 | 6.01 | 4.30 | 10.18 | 5.02 |  |
| SIMHIGH (n=2) | 7.69 | 7.43 | 9.65 | 12.88 | 9.87 | 17.77 | 12.23 |  |
| nucleotide (n=1) | 3.91 | 3.66 |  | 17.40 |  | 8.37 | 6.94 |  |

GCM keeps each subset alignment intact (the final alignment restricted to a subset is that subset's alignment), so
a better subset aligner should give a better MAGUS alignment, and it does in direction almost everywhere: on
BAliBASE MUSCLE5, ProbCons and Clustal Omega align MAGUS's subsets 0.2-0.4 points better than L-INS-i and the
final alignment improves by 0.3-0.5; on the simulated sets the same tools align the subsets 1-10 points worse and
the final alignment is 1-9.5 points worse. G-INS-i is the only method that is never clearly worse on subsets.

## Runtime

Total CPU seconds of the 25 subset alignments per dataset (mean over the set); ratio to L-INS-i.

| set | linsi | ginsi | muscle5 | famsa | probcons | clustalo | kalign | prank |
|---|---|---|---|---|---|---|---|---|
| BAliBASE | 81 (1x) | 77 (0.95x) | 317 (3.9x) | 24 (0.3x) | 627 (7.7x) | 45 (0.55x) | 17 (0.21x) |  |
| HomFam | 215 (1x) | 202 (0.94x) | 1478 (6.9x) | 6 (0.03x) | 1995 (9.3x) | 38 (0.17x) | 6 (0.028x) |  |
| SIMMOD | 164 (1x) | 178 (1.1x) | 1100 (6.7x) | 27 (0.17x) | 1560 (9.5x) | 53 (0.32x) | 12 (0.071x) |  |
| SIMHIGH | 234 (1x) | 251 (1.1x) | 1267 (5.4x) | 26 (0.11x) | 1632 (7x) | 64 (0.27x) | 12 (0.053x) |  |
| nucleotide | 1184 (1x) | 1306 (1.1x) |  | 71 (0.06x) |  | 313 (0.26x) | 568 (0.48x) |  |

GCM merge wall time (4 threads), L-INS-i subsets: median 27 s, range 1-54 s.

Control (merge of MAGUS's own subsets vs re-aligned L-INS-i subsets, SP error %):

| dataset | MAGUS (original) | merge-linsi |
|---|---|---|
| BBA0039 | 4.69 | 4.66 |
| BBA0067 | 26.28 | 26.40 |
| BBA0101 | 29.18 | 28.88 |
| BBA0154 | 21.66 | 21.76 |
| BBA0190 | 23.40 | 23.43 |
| HF_Acetyltransf | 26.76 | 26.03 |
| HF_PDZ | 16.14 | 16.61 |
| HF_aat | 15.16 | 15.07 |
| HF_rrm | 18.57 | 18.41 |
| HF_sdr | 23.65 | 23.35 |
| HF_zf-CCHH | 12.25 | 12.25 |
| SIMMOD_R1 | 12.71 | 12.79 |
| SIMMOD_R2 | 10.99 | 11.06 |
| SIMHIGH_R1 | 25.04 | 25.00 |
| SIMHIGH_R2 | 23.51 | 23.63 |
| 1000M2 | 8.23 | 8.24 |

All jobs shared the 4-core VM with other jobs of this pilot (MAGUS, PASTA), so wall times are inflated; CPU
seconds are the fair comparison. In end-to-end MAGUS the subsets are a minor part of the cost (MAGUS's own
backbones are 10 L-INS-i runs on 200 sequences each; the merge takes 1-54 s). Swapping in MUSCLE5 or ProbCons
multiplies the subset stage by 4-10x; FAMSA, Kalign and Clustal Omega make it nearly free, but they are the
least accurate.

## Trees (simulated proteins)

FastTree 2.1.11 `-lg -gamma`; normalized RF (%) vs the true 1,000-taxon tree (`code/trees.py`).

| dataset | merge_ginsi | merge_linsi | merge_muscle5 | merge_probcons | true |
|---|---|---|---|---|---|
| SIMHIGH_R1 | 6.82 | 10.13 | 10.63 | 8.83 | 6.12 |
| SIMHIGH_R2 | 8.53 | 6.72 | 8.32 | 7.22 | 5.22 |
| SIMMOD_R1 |  | 6.32 | 6.12 | 6.22 | 6.52 |
| SIMMOD_R2 | 7.02 | 7.52 | 7.82 | 7.82 | 6.52 |

## PASTA 1.8.3 with its built-in aligners, and whole-dataset baselines

PASTA rows are one PASTA iteration (`--iter-limit 1`) so that the aligner variants fit the time budget; PASTA
re-estimates its tree between iterations, and the MAGUS paper's default is 3. ProbCons on PASTA's 200-sequence
subproblems and the PASTA-bundled 2010 PRANK are very slow; runs that hit the 1 h cap are marked fail/timeout.
Whole-dataset MUSCLE5 `-align` hit its 20-min cap on BBA0039, BBA0067 and SIMMOD_R1 (machine shared with other jobs) and was
not run elsewhere for lack of time; on HF_rrm (`-super5`, 2,000 sequences) it finished. Regressive T-Coffee ran
only in a non-standard configuration (NJ guide tree, Clustal Omega children), so its very poor numbers here should
not be read as the published method's accuracy. MUSCLE5 vs MAGUS on whole datasets therefore remains open.

| dataset | MAGUS (L-INS-i subsets) | best MAGUS variant | famsa | muscle5 | pasta-mafft-it1 | pasta-prank-it1 | pasta-probcons-it1 | regressive |
|---|---|---|---|---|---|---|---|---|
| BBA0039 | 4.66 | 4.45 (clustalo) | 5.49 (268 s) | fail/timeout | 4.48 (542 s) | fail/timeout | 4.56 (1752 s) | 18.30 (33 s) |
| BBA0067 | 26.40 | 25.33 (probcons) | 32.04 (343 s) | fail/timeout |  |  |  | 36.14 (22 s) |
| BBA0101 | 28.88 | 28.11 (muscle5) | 35.24 (284 s) |  |  |  |  | 38.30 (25 s) |
| BBA0154 | 21.76 | 21.16 (clustalo) | 25.26 (255 s) |  |  |  |  | 28.19 (25 s) |
| BBA0190 | 23.43 | 22.98 (clustalo) | 33.67 (436 s) |  |  |  |  | 31.52 (77 s) |
| HF_Acetyltransf | 26.03 | 26.03 (linsi) | 51.32 (20 s) |  |  |  |  |  |
| HF_PDZ | 16.61 | 15.39 (ginsi) | 12.32 (21 s) |  |  |  |  |  |
| HF_aat | 15.07 | 14.23 (ginsi) | 15.89 (124 s) |  |  |  |  |  |
| HF_rrm | 18.41 | 16.98 (kalign) | 18.60 (24 s) | 24.71 (155 s) |  |  |  | 24.71 (41 s) |
| HF_sdr | 23.35 | 21.93 (ginsi) | 24.51 (128 s) |  |  |  |  |  |
| HF_zf-CCHH | 12.25 | 9.32 (kalign) | 10.06 (24 s) |  |  |  |  | 15.27 (32 s) |
| SIMMOD_R1 | 12.79 | 12.71 (orig) | 19.99 (282 s) | fail/timeout |  |  |  | 54.00 (60 s) |
| SIMMOD_R2 | 11.06 | 10.86 (ginsi) | 16.09 (295 s) |  |  |  |  |  |
| SIMHIGH_R1 | 25.00 | 24.19 (ginsi) | 38.22 (374 s) |  |  |  |  | 71.50 (68 s) |
| SIMHIGH_R2 | 23.63 | 22.92 (ginsi) | 29.55 (372 s) |  |  |  |  |  |
| 1000L1 |  |  | 33.08 (1558 s) |  |  |  |  |  |
| 1000M2 | 8.24 | 7.97 (ginsi) | 43.59 (922 s) |  | 12.96 (2029 s) | fail/timeout | fail/timeout | 81.16 (275 s) |
| 16S.M |  |  | 16.50 (555 s) |  |  |  |  |  |
| RNASim |  |  | fail/timeout |  |  |  |  |  |

## Nucleotide side (secondary)

Incomplete (time): merge-only results only for ROSE 1000M2; subset-level results for 1000M2, 1000L1 and part of
RNASim. 16S.M was not reached. MUSCLE5 was dropped on nucleotides (~400-500 CPU s per 1,000-bp subset).

| | L-INS-i | G-INS-i | FAMSA | Clustal Omega | Kalign | PRANK |
|---|---|---|---|---|---|---|
| 1000M2 MAGUS final SP error % (merge-only) | 8.24 | **7.97** | 23.01 | 13.26 | 11.52 | merge not finished |
| 1000M2 subset SP error % | 3.91 | 3.66 | 17.4 | 8.37 | 6.94 | 45.3 |
| 1000L1 subset SP error % | 2.30 | 2.48 | 17.5 | 7.45 | 6.98 | 57.1 |
| RNASim subset SP error % | 4.39 | 4.23 | 37.0 | 5.13 | 6.48 | 5.62 (6 of 25 subsets) |
| subset CPU s per dataset (1000M2) | 1,184 | 1,306 | 71 | 313 | 569 | 5,874 |

PRANK (v.250331, default options; `+F` and `-DNA` gave the same) over-gaps ROSE subsets badly: its subset
alignments are ~3x longer than the true alignment, with 45-57% subset SP error, at 5x L-INS-i's CPU. On RNASim it is
sane but still worse than L-INS-i. So Prank is not a candidate MAGUS base method on these data; this agrees with the
protein result that MAFFT's consistency-based modes are hard to beat on subsets of ~40 sequences.

## Methods

### Design: change only the subset aligner

For each dataset, one MAGUS run with the paper's flags (25 subsets, PASTA-style decomposition with a
300-sequence skeleton, 10 backbones of 200 sequences aligned by MAFFT L-INS-i, GCM graph from the backbones, MCL,
minclusters trace, `-f 4`) supplies the 25 subsets and the 10 backbone alignments. BAliBASE and nucleotide
runs reuse earlier MAGUS runs with the same flags (`cs581/experiments/runs/*/inputs.tar.xz`); HomFam and AliSim
runs are new (`code/basemeth.py prep`). Each subset is then re-aligned from its unaligned sequences by each base
method (one thread per subset, 4 subsets at a time), and the GCM merge alone
(`gcmx.run_magus -s SUBSETS -b BACKBONES`, same backbones, same MCL / minclusters flags, 4 threads) produces the
final alignment. So the only difference between variants is the subset aligner. `linsi` re-runs MAGUS's own
subset command (`mafft --localpair --maxiterate 1000 --ep 0.123`) and is the control; `orig` merges MAGUS's
original subset alignments untouched (checks that re-running L-INS-i reproduces MAGUS).

Base methods (single thread, default options unless stated):

| key | tool | version | command |
|---|---|---|---|
| linsi | MAFFT L-INS-i (MAGUS-bundled binary) | 7.450 | `--localpair --maxiterate 1000 --ep 0.123 --anysymbol` |
| ginsi | MAFFT G-INS-i | 7.450 | `--globalpair --maxiterate 1000 --anysymbol` |
| muscle5 | MUSCLE | 5.3 | `-align` (PPP, no ensemble) |
| famsa | FAMSA | 2.4.1 | default |
| probcons | ProbCons | 1.12 | default (2 consistency passes, 100 refinement passes) |
| clustalo | Clustal Omega | 1.2.4 | default |
| kalign | Kalign | 3.6.0 | default |
| prank | PRANK | v.250331 (bioconda 251117) | default (`-f=fasta`, no `+F`), nucleotides only |

Scoring: FastSP (SPFN, SPFP, TC); SP error = (SPFN+SPFP)/2 in percent. HomFam alignments are scored on the
Homstrad seed sequences only (the estimate restricted to them, all-gap columns removed). Paired difference =
variant minus `linsi`, in SP-error points; W/T/L with a tie band of |d| < 0.05 points; two-sided Wilcoxon
signed-rank (scipy); Holm over the four pre-registered protein candidates. Subset-level accuracy = each
method's 25 subset alignments scored against the reference restricted to each subset, pooled by homology
counts.

Runtime: per-subset wall and CPU seconds (children rusage) of the aligner; merge wall/CPU. All runs shared a
4-core / 15 GB VM with other jobs of this pilot (MAGUS runs, PASTA), so **wall times are inflated by contention;
CPU seconds are the fair cost comparison**.

Whole-dataset baselines (`code/baselines.py`, 4 threads, wall cap): FAMSA 2.4.1 default; MUSCLE 5.3 `-align`
(`-super5` above 1,000 sequences); regressive T-Coffee 12.00.7fb08c2 `-reg -reg_nseq 100 -reg_tree nj
-reg_method clustalo_msa` (the mBed guide tree hung in this bioconda build and `famsa_msa` is not available in it,
so this is not the published configuration). PASTA 1.8.3 (commit 738bec5, the MAGUS paper's version) with
`--aligner mafft|probcons|prank`, its bundled tools (MAFFT 7.149b, ProbCons, PRANK v.100311 with `+F`), OPAL
merger, FastTree, `--iter-limit 1` for the aligner comparison (3 = default for the reference run).

## Data

| set | datasets | size | reference | MAGUS inputs |
|---|---|---|---|---|
| BAliBASE RV100 (length-filtered, `cs581/data/balibase_clean`) | BBA0039, 0067, 0101, 0154, 0190 | 274-732 seqs | structural, all sequences | cached MAGUS runs (`cs581/experiments/runs/*_R0`) |
| HomFam (SALMA/EMMA release IDB-2567453, 2,000-seq subsamples, seed 1, `code/make_homfam.py`) | aat, sdr, rrm, PDZ, Acetyltransf, zf-CCHH | 2,000 seqs | Homstrad seeds only (6-20 seqs) | new MAGUS runs |
| AliSim proteins (as in the protbench branch: Yule tree, LG+G4, 300 aa, indel 0.05/0.05 POW 1.7/40; seeds 100r+7 / 100r+13) | SIMMOD R1-R2 (mean branch 0.06), SIMHIGH R1-R2 (0.10) | 1,000 seqs | true alignment and tree | new MAGUS runs |
| ROSE 1000M2, 1000L1 R0; RNASim 1000 R0; 16S.M R0 (MAGUS paper data) | 4 | 740-1,000 | true / CRW | cached MAGUS runs |

AliSim R3/R4 (in PREREG scope) were dropped for time: each MAGUS run took ~20-40 min under load.
I could not verify that the regenerated AliSim files are byte-identical to protbench's (no checksums were kept);
same commands and seeds, IQ-TREE 3.1.4, and MAGUS's error is in the same range (SIMMOD_R1 12.7% here vs 13.9% there;
MAGUS draws differ).

## Plan if pursued (4 weeks)

- **Week 1**: rerun the merge-only swap on more simulated replicates (AliSim R3-R10 at two or three rates; ROSE
  1000M/L) with 2-3 MAGUS draws each, to estimate noise properly; add G-INS-i and E-INS-i as base methods.
  Re-use `code/basemeth.py` unchanged.
- **Week 2**: trees: FastTree / IQ-TREE on every variant for all simulated sets (is the G-INS-i tree gain real?);
  PASTA with ProbCons/Prank on 2-3 simulated sets with 3 iterations (needs ~1-2 h per run on 4 cores).
- **Week 3**: whole-dataset comparison the literature lacks: MAGUS vs MUSCLE5, FAMSA2, regressive T-Coffee
  (fix the mBed hang: newer T-Coffee or a container), PASTA, on the same datasets; runtime under an idle machine.
- **Week 4**: write-up: GCM inherits subset accuracy; BAliBASE-vs-simulation disagreement; recommended base
  method per data type.

## Risks / caveats

- **BAliBASE vs simulation disagree** (as in the bbtool/protbench sessions): a project must report both and not
  pick the favorable benchmark. AliSim with LG+G4 and Zipfian indels favors MAFFT-like aligners.
- **One MAGUS draw per dataset**; paired merge-only differences remove decomposition noise, but re-running
  L-INS-i on the same subsets already moves results by up to 0.3 (BAliBASE) / 0.7 (HomFam) points.
- **HomFam is scored on 6-20 seeds** (often 1-2 subsets contain >= 2 seeds): per-family numbers are noisy.
- **Timing under contention**: CPU seconds are reliable, wall seconds are not.
- **Tool builds**: regressive T-Coffee could not run its published configuration here (mBed hang, no famsa_msa);
  PASTA's bundled MUSCLE is 3.8 and its PRANK is from 2010, so "MUSCLE5/PRANK in PASTA" needs code changes.
- **Scope cuts**: AliSim R3/R4 dropped; MUSCLE5 not run on nucleotide subsets (~400-500 CPU s per 1,000-bp
  subset); PASTA only on BBA0039 and 1000M2 with 1 iteration.
