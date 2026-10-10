# Pilot: merging two alignments (exact MWT DP vs profile–profile mergers)

Pilot for one candidate CS581 project, run 2026-10-09/10 on a 4-core, 15 GB machine. The pilot ran about 4.5 h of wall clock; the oracle and 1000-taxon runs share the 4 cores with other jobs. Code is in `code/`, raw results in `results/` (`merge.jsonl`, `pasta_e2e.jsonl`) and tables in `results/merge_summary.md`. Prior-art notes with DOIs are in `prior_art.md`.

**Verdict: not promising as a 4-week method project. Borderline as an evaluation study.** Details are in §7.

## 1. Question

The question comes from the first-day lecture, slide "Open problems (and possible course projects)", item *"Multiple Sequence Alignment: Merging two alignments"*. It is item 4 in `literature/slide_open_problems.md`.

PASTA merges pairs of adjacent subset alignments with OPAL (DNA) or MUSCLE (protein), then applies a transitivity merge. MAGUS replaced this step with GCM, which merges all subsets at once. This repo has an exact two-alignment dynamic program for the MWT objective (`code/gcmx/progressive.py`, trace `progdp`). Its edge weights count how many backbone alignments put two columns together, as in MAGUS and WITCH-NG. For many-way merges it was no better than GCM.

The pilot asks two questions:
- When we merge exactly two alignments, is progdp (or GCM) more accurate than profile–profile aligners?
- If so, does that carry over to PASTA?

## 2. Prior art (full list with DOIs: `prior_art.md`)

**Exact merging of two alignments is old.**
- Kececioglu & Zhang, *Aligning alignments*, CPM 1998 (doi:10.1007/BFb0030790).
- Ma, Wang & Zhang 2003 (doi:10.1007/3-540-44888-8_19).
- Kececioglu & Starrett, RECOMB 2004 (doi:10.1145/974614.974626). Aligning two alignments is NP-complete under sum-of-pairs scoring with exact affine gap counts, and polynomial without exact gap counts.
- OPAL (Wheeler & Kececioglu 2007, doi:10.1093/bioinformatics/btm226) solves it exactly in practice. OPAL is **PASTA's default merger for DNA**.

**The other profile mergers are heuristics.**
- PASTA (doi:10.1089/cmb.2014.0156) and SATé-II (doi:10.1093/sysbio/syr095).
- MUSCLE 3 `-profile` (doi:10.1093/nar/gkh340). MUSCLE5 (doi:10.1038/s41467-022-34630-w) has no user-facing profile–profile command.
- MAFFT `--merge`, which is documented only in the manual; the MAFFT 7 paper is doi:10.1093/molbev/mst010.
- T-Coffee/M-Coffee (doi:10.1093/nar/gkl091), Clustal Omega `--profile1/2`, FAMSA2 and TWILIGHT (doi:10.1093/bioinformatics/btaf212).

**The MWT objective with backbone evidence is also known.**
- MAGUS (doi:10.1093/bioinformatics/btaa992).
- Zaharias, Smirnov & Warnow, TCBB 2023 (doi:10.1109/TCBB.2022.3191848). They define MWT-AM and note in Theorem 1 that it is polynomial for a fixed number of alignments B (O(2^B L^B B^2)), so B = 2 is a quadratic DP. They did not study that case.
- WITCH-NG (doi:10.1093/bioadv/vbad024) uses exactly this DP for one sequence against an alignment, with zero gap cost.
- EMMA (doi:10.1186/s13015-023-00247-x) merges by transitivity through a fixed constraint alignment.

**Novelty.** The algorithm is not new: it is the B = 2 case of a published theorem and the m × n generalisation of WITCH-NG's DP. What is untested is the combination, i.e. MWT with backbone evidence used as PASTA's pairwise merger and compared with OPAL and MUSCLE. I found no 2022–2026 paper that does this.

## 3. Data

The data are the MAGUS paper's simulated sets (Illinois Data Bank, `/opt/data/Datasets`):
- ROSE 1000M2, 1000M4, 1000L2 and 1000S2, replicates R0–R1;
- RNASim 1000, R0–R1.

That gives 10 replicates of 1000 sequences, about 1000 bp each (RNASim about 1550 bp).

Each replicate is split into two halves at the centroid edge of the **true** tree (446/554 up to 500/500 taxa). Only the true-tree split was used; there was no time to add a guide-tree split.

## 4. Method

Each merger gets the same two input alignments, and every merger keeps the input columns intact. This was checked on every run (`constraintsKept`).

| merger | what is run |
|---|---|
| `opal` | `java -jar opal.jar --in A --in2 B --align_method profile`, exactly as PASTA calls it |
| `muscle3` | MUSCLE 3.8.31 `-profile -in1 A -in2 B`, as PASTA calls it |
| `mafft-merge` | MAFFT 7.505 `--merge` (default progressive strategy) |
| `gcm` | MAGUS GCM (MCL + minclusters trace) on the two subsets, with MAGUS(Fast) flags |
| `progdp` | exact two-alignment MWT DP on the same MAGUS graph, built from the same backbones as `gcm` |

**Evidence.** There are 10 backbones per replicate, each a MAFFT L-INS-i alignment of a random sample of 100 sequences, built by MAGUS once per replicate. The paper uses 200-sequence backbones; I used 100 because L-INS-i on 200 sequences of about 1000 bp costs about 6 core-minutes per backbone, which did not fit the budget. `gcm` and `progdp` always see identical evidence.

**Input conditions.**
- `oracle`: the true alignment restricted to each 1000-taxon half. All error is then merger error.
- `fftnsi`: MAFFT FFT-NS-i (`--retree 2 --maxiterate 2`) on each half. L-INS-i on 500 sequences of about 1000 bp takes more than 30 core-minutes per half, which was over budget.
  - **Caveat:** FFT-NS-i fails on these data. The halves themselves have 75–95% SPFN on 1000M2, 1000L2 and 1000S2, so this condition mostly measures input error. It is reported but should not be relied on.
- `oracle200` / `linsi200`: **PASTA-scale merges**. These use 200 random taxa from each half, which is PASTA's maximum subproblem size, with either the true sub-alignment or MAFFT L-INS-i (PASTA's subset aligner). The backbones are the replicate's backbones restricted to the 400 taxa, so each holds about 40 sequences. This is the most realistic condition.

**Metrics.**
- FastSP SPFN and SPFP against the true alignment of the merged taxa, through `gcmx/score.py`; *error* = (SPFN + SPFP)/2.
- SPFN/FP restricted to **cross-half** homology pairs. These are the only pairs a merger decides, so this isolates merger error.
- Statistics: paired by replicate, two-sided Wilcoxon signed-rank test, mean difference, W/T/L with a tie band of 0.1 percentage points.

## 5. Results (error in %)

### 5a. Oracle: merger error only (1000-taxon halves, n = 10, except n = 8 for gcm/progdp)

| merger | 1000L2 | 1000M2 | 1000M4 | 1000S2 | RNASim | all | cross-half FN / FP | merge time |
|---|---|---|---|---|---|---|---|---|
| opal | 0.15 | 0.28 | 0.05 | 0.38 | 1.25 | **0.42** | 0.85 / 0.87 | 2 s |
| progdp | 0.57 | 0.80 | 0.07 | 0.77 | 2.14 | 0.73 | 1.50 / 1.44 | 2–3 s + evidence |
| gcm | 0.63 | 0.93 | 0.19 | 0.76 | 2.54 | 0.85 | 2.14 / 1.33 | 2–60 s + evidence |
| mafft-merge | 4.88 | 10.47 | 1.93 | 8.68 | 10.10 | 7.21 | 19.05 / 11.00 | 6 s |
| muscle3 | 6.66 | 19.07 | 3.17 | 7.65 | 28.24 | 12.96 | 31.99 / 24.87 | 5 s |

Paired comparisons against OPAL:

| merger | mean diff (points) | W/T/L | p |
|---|---|---|---|
| progdp | +0.36 | 1/2/5 | 0.078 |
| gcm | +0.49 | 1/1/6 | 0.055 |
| MAFFT `--merge` | +6.8 | 0/0/10 | 0.002 |
| MUSCLE | +12.5 | 0/0/10 | 0.002 |

progdp vs GCM on the same graph: −0.13 points, W/T/L 4/4/0, p = 0.11.

At the PASTA scale (`oracle200`, 200 + 200 true sub-alignments, n = 10):
- OPAL 0.44, progdp 0.77, GCM 0.92, MAFFT `--merge` 5.87 and MUSCLE 10.40.
- progdp vs OPAL: +0.33 points, W/T/L 2/2/6, p = 0.11.
- GCM vs OPAL: +0.48 points, W/T/L 1/2/7, p = 0.014.
- progdp vs GCM: −0.15 points, W/T/L 7/2/1, p = 0.037.

The ranking is the same as for the 1000-taxon halves.

**Reading.** Given correct inputs, OPAL's exact sum-of-pairs merge is essentially perfect, with less than 1% cross-pair error. The MWT DP is no better: it is slightly worse, though not significantly. MUSCLE 3 and MAFFT `--merge` are 7–13 points worse, so the large accuracy differences between "profile–profile" mergers come from those two tools, not from OPAL.

### 5b. PASTA-scale merges with estimated inputs (`linsi200`, 200 + 200 taxa, L-INS-i halves)

n = 10 replicates (2 per model condition).

| merger | 1000L2 | 1000M2 | 1000M4 | 1000S2 | RNASim | all |
|---|---|---|---|---|---|---|
| opal | 8.31 | 24.05 | 1.18 | 10.50 | 10.21 | 10.85 |
| progdp | 8.07 | 20.89 | 1.20 | 10.44 | 9.88 | **10.10** |
| gcm | 8.07 | 20.92 | 1.25 | 10.45 | 9.92 | 10.12 |
| mafft-merge | 10.43 | 26.75 | 2.21 | 13.97 | 13.99 | 13.47 |
| muscle3 | 9.50 | 30.35 | 1.61 | 13.23 | 29.51 | 16.84 |

Paired comparisons against OPAL:

| merger | mean diff (points) | W/T/L | p |
|---|---|---|---|
| progdp | −0.75 | 6/4/0 | 0.037 |
| gcm | −0.73 | 4/6/0 | 0.16 |
| MAFFT `--merge` | +2.6 | 0/0/10 | 0.002 |
| MUSCLE | +6.0 | 0/0/10 | 0.002 |

**Reading.** With realistic, erroneous inputs, the evidence-based mergers are never worse than OPAL and sometimes better. Almost all of the gain comes from 1000M2, about −3 points on both replicates. On RNASim the gain is −0.3 points, and on M4, L2 and S2 the difference is under 0.25 points. The gain is significant at n = 10, but only just (p = 0.037).

The likely mechanism is that the backbones add information the inputs lack: 10 independent L-INS-i alignments that also cover the merged taxa. Where the input sub-alignments are wrong, the backbone votes can still pair the right columns, whereas OPAL can only score the columns it is given. That information is not free (§6). In the oracle condition there is nothing to correct, and OPAL wins.

### 5c. FFT-NS-i halves (`fftnsi`, 1000-taxon halves)

GCM and progdp are better than OPAL here: −2.9 points (W/T/L 6/1/1, p = 0.039) and −1.6 points (6/2/0, p = 0.023). However, the inputs are 75–95% wrong on M2, L2 and S2, so absolute errors are 70–90% for every merger. Only 1000M4 (all mergers at 2.4–4.3%) and RNASim (19–44%) are interpretable. I do not draw conclusions from this condition.

### 5d. End to end in PASTA (sanity check only)

PASTA 1.9.x (bioconda) was run with:
- 1 iteration;
- maximum subproblem size 200, centroid decomposition, L-INS-i on the subsets, FastTree;
- a shared starting tree (FastTree on a MAFFT FFT-NS-2 alignment). PASTA's own HMMER-based initial alignment crashes in this build for more than about 200 sequences.

The progdp arm plugs in through PASTA's `muscle` merger slot (`code/pasta_shim.sh` → `pasta_merger.py`). Each pairwise merge builds 10 L-INS-i backbones of 50 sequences from the two sides, runs progdp, and restores the input letters. Inputs are 500-taxon random subsamples of 1000M2 replicates, kept small for budget.

| replicate | PASTA (OPAL) error | PASTA (progdp) error | diff | runtime OPAL / progdp |
|---|---|---|---|---|
| 1000M2_R10 (500 taxa) | 46.61 | 46.45 | −0.16 | 806 s / 2537 s |

Only one pair finished in the budget; the second pair's progdp arm was still running when this was written. With n = 1 no test is possible. This is a tie within noise, at 3.1× the runtime. The very high absolute error comes from one iteration started from a poor FFT-NS-2 tree. Strictly, the precondition of task 3 ("one merger clearly best") was not met, because OPAL is already PASTA's default and nothing beat it significantly. The run only checks that the plug-in works and does not show an obvious gain.

## 6. Runtime

| step | cost |
|---|---|
| OPAL / MUSCLE / MAFFT merge of two halves (1000 taxa) | 1–8 s |
| progdp DP given the graph | 1–3 s (more on RNASim) |
| GCM given the graph (MCL + trace) | 1–60 s |
| **backbone evidence (10 × L-INS-i on 100 sequences)** | **15–60 single-core minutes per replicate** (some measurements under 2× CPU oversubscription) |
| PASTA, 500 taxa, 1 iteration: OPAL vs progdp merger | 806 s vs 2537 s (3.1×) |

The evidence-based mergers cost roughly 1000× more than OPAL per merge, because they align new backbones. Inside PASTA this cost recurs for every spanning-tree edge and every iteration.

## 7. Verdict, possible 4-week project, risks

**Verdict: not promising** as a "better merger" project.

- With correct inputs, OPAL, PASTA's default for DNA, is already essentially optimal (0.42% error, under 1% cross-pair error). There is no merger error left to remove.
- With realistic estimated inputs, exact MWT with backbone evidence is 0.75 points better on average (W/T/L 6/4/0, p = 0.037, n = 10). Nearly all of that comes from one model condition (1000M2, about −3 points). It costs about 1000× OPAL's compute per merge.
- The end-to-end PASTA check showed a tie (n = 1).
- The algorithm is the B = 2 case of a published theorem (Zaharias et al. 2023) and an obvious extension of WITCH-NG.

The genuine finding is about which tools PASTA should use. MUSCLE 3 and MAFFT `--merge` are much worse than OPAL as mergers (5–13 points with oracle inputs, 2.6–6.0 points with L-INS-i inputs, p = 0.002). PASTA uses MUSCLE for protein data.

**If someone still wants a 4-week project here**, the version with a chance is an evaluation study, not a new merger:
1. Week 1: evaluate OPAL vs MUSCLE vs MAFFT `--merge` vs GCM/progdp as PASTA mergers on **protein** data (BAliBASE / HomFam). PASTA defaults to MUSCLE there, and our DNA results suggest MUSCLE loses 5–13 points of merge accuracy. Test whether OPAL with a protein scoring matrix closes that gap.
2. Week 2: a cheap evidence source for progdp, so it does not need new backbones. Options are reusing the alignments PASTA already has (the previous iteration's full alignment restricted to the two subsets), or HMM posteriors as in WITCH-NG.
3. Week 3: hybrid objective. Run OPAL's sum-of-pairs merge, but let backbone or HMM votes override it only where the input sub-alignments look unreliable. The M2 result says the gain comes from correcting input errors.
4. Week 4: PASTA end to end with 3 iterations on 10+ replicates per condition, tree error included.

**Risks.**
- The headroom is small (OPAL error is already near 0 on correct inputs).
- PASTA end to end needs about 1–2 h per run on 4 cores at 1000 taxa, so 10 replicates × 2 arms × 3 iterations does not fit a laptop budget.
- The bioconda PASTA build has an outdated dendropy and a broken HMMER initial alignment (workarounds are in `code/`).
- On protein data, OPAL may be slow or unavailable in PASTA's bundle.

## 8. Deviations, bugs found and caveats

- Backbones use 100 sequences instead of MAGUS's 200, for budget. The estimated 1000-taxon halves use FFT-NS-i instead of L-INS-i, for budget, and this turned out to be uninformative (§5c). Two replicates per model condition; the true-tree split only.
- Two pipeline bugs were found and fixed before the reported runs:
  1. MAFFT writes lower case and MAGUS matches letters case-sensitively, so GCM/progdp on MAFFT inputs or backbones silently found no matches and returned unmerged alignments. Inputs and backbones are now upper-cased.
  2. MAGUS silently skips a merge if its output file already exists, so outputs are now deleted before each run.
  
  Rows affected by either bug were deleted and recomputed.
- The GCM time in the oracle condition includes building the backbones once per replicate; the other conditions reuse them.
- n is small: 2 replicates per model condition, 8–10 replicates per comparison. p-values are descriptive.
