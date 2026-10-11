@@VERDICT
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
  candidate is worse on all 4 replicates, except MUSCLE5 on SIMHIGH_R2 (-0.06); Clustal Omega costs ~7 points. On HomFam (seed-scored, noisy) all three "BAliBASE winners" lose on 6/6 families by 2-5 points. This
  mirrors the earlier Clustal-backbone result (bbtool / protbench): BAliBASE rewards conservative, low-SPFP
  aligners; indel-rich simulated data penalizes them.
- **Mechanism is simple:** GCM preserves the subset alignments, so MAGUS's final error tracks the subset
  aligner's own accuracy on the 40-80-sequence subsets (subset-level table below). L-INS-i (and G-INS-i) are the
  most accurate subset aligners on everything but BAliBASE.
- **Only G-INS-i is interesting** (not pre-registered, exploratory): never clearly worse on simulated data,
  -0.76 on both SIMHIGH replicates, and on SIMHIGH_R1 its MAGUS tree is much closer to the true tree (FastTree nRF
  6.8% vs 10.1% for L-INS-i; one replicate). HomFam is mixed (3 wins, 3 losses, mean +0.97 driven by
  Acetyltransf +8, a 6-seed family).
- **Cost.** MUSCLE5 and ProbCons make the subset stage 4-10x more expensive (e.g. ~1,100-1,600 vs ~165 CPU s
  per 1,000-sequence AliSim dataset); FAMSA, Kalign and Clustal make it nearly free but are the least accurate.
  The subset stage is a small part of MAGUS's total time either way (MAGUS's 10 L-INS-i backbones dominate).
- **Whole-dataset tools are worse than MAGUS on BAliBASE and AliSim**: FAMSA2 alone is 0.8-13 points worse there
  (but better than MAGUS on 2 of 6 HomFam families, PDZ and zf-CCHH); regressive T-Coffee (NJ tree + Clustal children; see caveat) 18.3% vs MAGUS 4.7% on BBA0039. See
  the baselines section for what finished.

**Verdict for a 4-week CS581 project: not promising** as "find a better base method for MAGUS's subsets" (the
pre-registered question has a clear negative answer on 15 datasets, and the only positive signal, on BAliBASE,
is within 2-4x noise and contradicted by simulations). It is a perfectly good *negative-result* project
(straightforward, cheap, clear story: GCM inherits subset accuracy; BAliBASE vs simulation disagree), and the
G-INS-i / tree observation is the one thread worth a week if the student wants a positive angle.

@@BASE
PASTA rows are one PASTA iteration (`--iter-limit 1`) so that the aligner variants fit the time budget; PASTA
re-estimates its tree between iterations, and the MAGUS paper's default is 3. ProbCons on PASTA's 200-sequence
subproblems and the PASTA-bundled 2010 PRANK are very slow; runs that hit the 1 h cap are marked fail/timeout.

@@DNA
__DNA_TEXT__

@@PLAN
- **Week 1**: rerun the merge-only swap on more simulated replicates (AliSim R3-R10 at two or three rates; ROSE
  1000M/L) with 2-3 MAGUS draws each, to estimate noise properly; add G-INS-i and E-INS-i as base methods.
  Re-use `code/basemeth.py` unchanged.
- **Week 2**: trees: FastTree / IQ-TREE on every variant for all simulated sets (is the G-INS-i tree gain real?);
  PASTA with ProbCons/Prank on 2-3 simulated sets with 3 iterations (needs ~1-2 h per run on 4 cores).
- **Week 3**: whole-dataset comparison the literature lacks: MAGUS vs MUSCLE5, FAMSA2, regressive T-Coffee
  (fix the mBed hang: newer T-Coffee or a container), PASTA, on the same datasets; runtime under an idle machine.
- **Week 4**: write-up: GCM inherits subset accuracy; BAliBASE-vs-simulation disagreement; recommended base
  method per data type.

@@RISKS
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
