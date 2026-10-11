# Pre-registration: equal-time, scale and estimated-alignment tests of two-step pipelines (CS581 pilot `fragscale`)

Written and pushed **before any tree in this pilot was estimated or scored** (2026-10-11).
Builds on `claude/cs581-fragml2` (code copied into `code/`, same arm names). Background there:
on 1000M1-HF with the true alignment, `place_ft_0.5_fix_rxfast` and `constr_ft_0.5` tie RAxML-NG
(ΔFN +0.44 / +0.02, n.s.) at ~1/3 of its CPU. The question here is whether they **beat** the best ML
method once time is equalised, at larger scale, and on estimated alignments.

## Tools (pinned in `code/install.sh`)
RAxML-NG 2.0.3, IQ-TREE 3.1.4, FastTree 2 (Ubuntu), EPA-ng 0.3.8 stock and patched
(`epa-ng-fix.patch` from `claude/cs581-epang`), gappa, UPP (SEPP, bioconda) and/or WITCH (pip).
All tree tools run with 1 thread; CPU = user+sys of all steps; peak RSS = max over steps.

## Anytime measurement
- **RAxML-NG anytime curve**: one full default search per replicate (GTR+G, one parsimony start,
  seed 1, 1 thread). RAxML-NG rewrites `<prefix>.raxml.lastTree.TMP` (current best tree) at every
  checkpoint; a poller copies it on every change together with the process' CPU time. "RAxML-NG@T" =
  the last such tree written at CPU ≤ T (the tree a run killed at T would have reported; topology only,
  so FN is exact). If no tree exists by T (still in the parsimony/initial phase), RAxML-NG@T = its
  parsimony starting tree (`.raxml.startTree`) when available, else missing (reported).
- **IQ-TREE 3 default anytime**: one default run per replicate (GTR+G, seed 1, 1 thread, `-cptime 10`),
  CPU capped at the same replicate's full RAxML-NG CPU; best tree read from each checkpoint
  (`.ckp.gz`) on change, with CPU time. IQ-TREE@T defined as above.
- Budgets reported: T = 1/3, 1/2, 1× of that replicate's full RAxML-NG CPU, and T = the pipeline's own CPU.

## Q1 (primary) — equal CPU, 1000M1-HF, true alignment
Data: M1HF R0–R4 = exact Park, Zaharias & Warnow 2021 inputs (IDB-7008049); R5–R9 = ROSE 1000M1
R5–R9 fragmented with `make_frag.py` (same protocol as fragml2). Planned n = 8 (R0–R7), minimum 6.

**Primary pipeline (single, declared now):** `place_ft_0.5_fix_rxfast` (FastTree backbone of sequences
≥ 0.5 × median length → RAxML-NG `--evaluate` → patched EPA-ng → gappa graft → RAxML-NG fast-mode polish).

**Primary test:** FN of `place_ft_0.5_fix_rxfast` vs FN of RAxML-NG@T_pipe, where T_pipe = that
replicate's total pipeline CPU, paired over replicates, two-sided exact Wilcoxon signed-rank (scipy;
zero differences dropped, "wilcox" zero_method), W/T/L with a 0.5-point tie band.
Interpretation fixed in advance: **two-step better at equal time** = mean ΔFN < −0.5 points and p < 0.05;
**equal** = |mean ΔFN| ≤ 0.5; otherwise worse/unclear.

## Secondary (descriptive; p-values unadjusted unless stated)
- Q1: RAxML-NG@{1/3, 1/2, 1}×full and IQ-TREE@{1/3, 1/2, 1}×RAxML-NG-full (IQ-TREE on R0–R4 at least),
  IQ-TREE `--fast`, `constr_ft_0.5` vs RAxML-NG@its own CPU; anytime curve plot (mean FN vs CPU fraction).
- **Q2 — RNASim 10K HF** (RNASim 10K R0 from the MAGUS Datasets.zip, fragmented with `make_frag.py`;
  analysis on the true alignment masked to columns with < 95% gaps among full-length sequences, as in
  fragml2/epang; the Smirnov & Warnow Dryad deposit is not reachable from this VM). Arms: FastTree,
  IQ-TREE `--fast`, RAxML-NG (default search, fixed budget: anytime curve up to a cap of 4 CPU-hours;
  the cap is the "generous budget"), `place_ft_0.5_{fix,stock}_{graft,ft,rxfast}`; BSCAMPP(e) placement
  only if time allows. FN, CPU, wall, peak RSS. R1 if time allows (n ≤ 2: descriptive only).
- **Q3 — estimated alignments.** 1000M1-HF R0–R4 (n = 5) and RNASim 1K HF R0 (one RNASim size), aligned
  with UPP (default, as in Smirnov & Warnow 2021) — WITCH instead if UPP fails to install/run. Insertion
  (lowercase) columns are removed (UPP's masked alignment) and all methods use the same alignment. Arms:
  RAxML-NG anytime, `place_ft_0.5_fix_rxfast`, `constr_ft_0.5`, FastTree. Same primary-style comparison
  (pipeline vs RAxML-NG@T_pipe), reported as secondary because n = 5 cannot reach p < 0.05 with Wilcoxon.
  Alignment time is reported separately (common to all arms).

## Priorities (4 cores, ~5 h of compute)
1. Q1 R0–R4: RAxML-NG anytime, primary pipeline, constr_ft_0.5; then IQ-TREE anytime R0–R4; Q1 R5–R7.
2. Q2 R0 (RAxML-NG with the 4 CPU-h cap runs in parallel with everything else).
3. Q3 M1HF R0–R4, then RNASim1K R0.
Arms not finished are reported as not run.

## Amendments
1. (2026-10-11 ~01:20 UTC, before any tree was scored; compute only.) UPP (SEPP 4.5.6 from bioconda)
   crashes in its own backbone step (`KeyError: 'user_options'` in `sepp/jobs.py`), so Q3 uses the
   pre-declared fallback, WITCH (MAGUS backbone, `code/align_witch.py`). WITCH's default backbone rule
   (±25% of the median of ALL lengths) is replaced by UPP's `-M 0.75` rule (±25% of the 3rd quartile),
   because on HF data the overall median falls between fragments and full-length sequences.
   Budget cuts (4 cores; MAGUS takes ~30 CPU-min per replicate): Q1 n = 6 (R0–R5; R6–R7 only if time);
   `constr_ft_0.5` and `base_iqfast` on Q1 R0–R2 only; IQ-TREE anytime on R0–R2; Q3 RNASim 1K and
   Q2 `base_iqfast` / `place_ft_0.5_fix_rxfast` (10K) moved to "if time allows". Nothing else changes.
2. (2026-10-11 ~02:45 UTC, AFTER seeing Q1 R0–R2 primary numbers; **post hoc, secondary only**.) RAxML-NG's
   anytime trajectory is very coarse: the first checkpoint after the parsimony start tree comes only after
   the first FAST SPR round (~11–15 CPU-min on M1HF), so "RAxML-NG@T_pipe" is often the parsimony tree.
   A user with a 1/3 budget would rather pick a faster RAxML-NG setting, so a second comparator is added:
   `base_rxfastmode` = RAxML-NG's own fast mode (1 parsimony start, `--opt-topology simplified
   --stop-rule kh-mult`), compared with the pipeline on FN and CPU (paired, Wilcoxon). The primary
   test and its interpretation are unchanged; the report states the result of both comparators.
3. (2026-10-11 ~02:55 UTC, before any Q3 alignment existed; compute only.) MAGUS needs > 75 min per replicate
   at the CPU share it gets here, so the WITCH backbone is aligned with MAFFT `--auto` (FFT-NS-i, 1 thread)
   instead of MAGUS. Fragments are still added with WITCH. Everything else in Q3 is unchanged.
