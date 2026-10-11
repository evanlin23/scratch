# WITCH-lite pilot: cheaper HMM selection for adding reads with WITCH

Overnight pilot, 2026-10-10. Code: `code/`. Results: `results/` (`tables.md` holds every table, the
per-query data is in `*.json`, and `acc_vs_time.png` is the summary figure). Everything below was run on a
4-core / 15 GB machine with nothing else running during timed runs. The two runs that overlapped another job
were re-timed, and the clean timings are used.

## Question and source open problem

WITCH (Shen, Park & Warnow, JCB 2022) adds each query to a backbone alignment as follows:

1. It scores the query against **every** HMM in the UPP-style ensemble with `hmmsearch --max`.
2. It turns those scores into weights (adjusted bit-scores).
3. It aligns the query to the top k = 10 HMMs and merges the results.

TIPP3 (Shen, Wedell, Pop & Warnow, PLOS CB 2025) states the open problem: *"The step in TIPP3 where reads are
added into the marker gene alignment using WITCH is the biggest contribution to runtime"*, and *"developing new
methods for this step that are substantially faster but not much less accurate than WITCH"* is *"a promising
direction"*. BSCAMPP (2025) adds that *"research into speeding up these alignments while maintaining high accuracy
is also merited"*. TIPP-SD (2026) defaults to BLASTN because WITCH *"increased the runtime substantially"*.

**Question:** can WITCH skip scoring most of the HMMs and stay within about 1 SPFN point of WITCH at under 20% of
its runtime?

## Prior art and novelty check

Full notes are in `REPORT_priorart.md`. In brief:

- **UPP2** (Park et al., *Bioinformatics* 2023) already does idea (a). It descends the decomposition tree,
  scoring two children per level, plus an "EarlyStop" variant. It does this for **top-1** HMM selection in UPP.
- **J-bandit** (Mazooji & Shomorony, *Bioinformatics* 2024) uses k-mer bandit sampling, also for top-1 in UPP.
  It does not mention WITCH or reads.
- **WITCH-NG** and **HMMerge** speed up or alter the *merge* step. Both still score every HMM.
- witch-msa 1.0.10 (read directly) has no option to restrict scoring.
- TIPP3 0.5 calls WITCH on the full marker alignment, so the decomposition is rebuilt on every run.
- TIPP3 always runs BLASTN to bin reads, so **every read already has a BLAST top hit**.
- **Not found anywhere:**
  - selection for WITCH's top-k weighting;
  - a beam search over the decomposition tree;
  - BLAST-guided HMM selection;
  - adaptive k.
- **Novelty is therefore modest.** This is an extension and benchmark of UPP2's search for WITCH and for reads,
  not a new algorithmic idea.

## Methods

**Data.** Five instances. Each uses the reference (true) alignment restricted to a random backbone, as TIPP3
uses curated reference alignments. Backbone sequences are drawn from those within 25% of the median length.
Queries are one random 150-bp window ("read") from each held-out sequence.

| instance | source | backbone | queries | HMMs in ensemble |
|---|---|---|---|---|
| 16S_b1000 | 16S.B.ALL (CRW, real) | 1000 | 1000 | 277 |
| 16S_b5000 | 16S.B.ALL | 5000 | 1000 | 1387 |
| rna_b1000 | RNASim 10k R2 | 1000 | 1000 | 283 |
| roseM1_R0 / R1 | ROSE 1000M1 (p-distance about 0.69) | 700 | 300 | 203 / 199 |

The backbone tree is FastTree `-nt -gtr`. The WITCH ensemble is the default
(`-A 10`, hierarchical centroid decomposition). The flat on-disk subsets are re-nested by set containment.

**Accuracy.**
- SPFN and SPFP are computed exactly on **query–backbone homology pairs**: the backbone is fixed and correct,
  so each query residue either lands in its true backbone column or loses all of its true pairs.
- The scorer was checked against FastSP on single-query alignments: shared, reference and estimated homology
  counts match exactly (e.g. 115684 / 136931 / 136046).
- Paired per-query differences against WITCH are reported as W/T/L and a Wilcoxon signed-rank p-value.
- **Placement:** EPA-ng places each read on the backbone tree (best-LWR edge). Delta error is the number of true
  bipartitions the tree loses when the read is added, compared against the true tree. This needs a true tree,
  so it is only computed for RNASim and ROSE.

**Methods compared.** All variants use WITCH's own alignment and merge stage (WITCH-NG mode). Each selection rule
writes a `weights.txt`, which WITCH reads in place of its own search.

- **WITCH:** the default run.
- **all/kN:** score all HMMs, then keep the top k (k = 1 is UPP-like) or use **adaptive k** (keep the fewest
  top HMMs whose weights sum to at least τ = 0.99 or 0.95).
- **hier / hier_es:** UPP2's descent and its EarlyStop variant, followed by WITCH weights over the scored HMMs.
- **beam-b:** keep the best b children per level (b = 2, 3, 4), scoring about 2b HMMs per level.
- **blastpath (+sib):** score only the root-to-leaf decomposition path of the read's top BLASTN hit (plus the
  siblings of the path nodes). The guide uses `-task blastn`, because megablast misses most RNASim/ROSE reads.
- **hybrid:** use the BLAST paths and siblings of the top 3 hits if the top hit scores at least 100 bits,
  otherwise use beam-3. Its selection was replayed from cached scores, and its selection time is **estimated**
  as BLAST time plus (scores per read) × (measured per-score cost of online beam-3).
- **BLASTN-only:** TIPP3-fast's alignment (megablast) and a sensitive `-task blastn` version.

**Scoring backend.** Online selectors call the same `hmmsearch --cpu 1 --max -E 99999999` binary WITCH uses, on
4 threads. My all-HMM scorer takes 243 / 277 / 1213 s against WITCH's 243 / 280 / 1217 s, so per-score cost
matches. pyhmmer 0.12 was 4.7× slower per score here and was dropped.

**Runtime.** Lite time = WITCH decomposition (same as WITCH) + selection + WITCH alignment stage. A "% excl.
decomposition" column is also reported, because a TIPP3 reference package could ship the decomposition
precomputed.

## Results

### Summary

Each cell is ΔSPFN versus WITCH in percentage points / runtime as a % of WITCH wall time.

| instance | HMMs | WITCH SPFN % | WITCH wall (search share) | beam-3: HMMs scored per read | selection speed-up vs WITCH search | beam-3 k=10 | hybrid | BLAST path+sib | hier (UPP2) | EarlyStop (UPP2) | BLASTN -task blastn |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 16S_b1000 | 277 | 0.54 | 322 s (75%) | 38 | 6x | +0.01 / 35% | -0.00 / 30% | -0.00 / 27% | +0.19 / 27% | +1.36 / 19% | +0.03 / 1% |
| 16S_b5000 | 1387 | 1.29 | 1419 s (86%) | 52 | 20x | +0.03 / 17% | +0.01 / 16% | +0.02 / 15% | +0.30 / 14% | +1.05 / 12% | **-0.81** / 1% |
| rna_b1000 | 283 | 2.72 | 368 s (76%) | 38 | 5x | +0.02 / 37% | +0.01 / 32% | +0.08 / 29% | +0.51 / 29% | +1.19 / 24% | +8.40 / 1% |
| roseM1_R0 | 203 | 7.85 | 61 s (65%) | 36 | 4x | +0.34 / 49% | +0.34 / 50% | +6.35 / 47% | +3.11 / 42% | +11.29 / 37% | +59.8 / 1% |
| roseM1_R1 | 199 | 10.03 | 58 s (66%) | 36 | 4x | +1.49 / 50% | +1.49 / 48% | +6.67 / 47% | +6.35 / 42% | +14.90 / 37% | +51.4 / 1% |

The paired statistics for these rows:

| instance | beam-3 W/T/L vs WITCH | Wilcoxon p |
|---|---|---|
| 16S_b1000 | 4/983/13 | 0.10 |
| 16S_b5000 | 2/985/13 | 0.009 |
| rna_b1000 | 0/993/7 | 0.018 |
| roseM1_R0 | 0/290/10 | 0.005 |
| roseM1_R1 | 0/291/9 | 0.008 |

- Losses are rare: 2–3% of reads. But a lost read is badly misaligned, because the descent went down the wrong
  subtree.
- Beam-4 narrows the gap on ROSE: +0.19 pp on R0 and +1.15 pp on R1.

### Placement delta error

Mean delta error per read; in brackets, the difference from WITCH and the Wilcoxon p.

| instance | WITCH | all k=1 | hier (UPP2) | beam-3 | BLAST path+sib | BLASTN -task blastn |
|---|---|---|---|---|---|---|
| rna_b1000 | 1.117 | 1.130 (+0.01, p=.46) | 1.157 (+0.04, p=.6) | 1.118 (+0.00, p=.32) | 1.114 (-0.00) | 1.903 (+0.79, p=3e-15) |
| roseM1_R0 | 0.830 | 0.790 | 1.123 (+0.29, p=.09) | 0.823 (-0.01) | 1.193 (+0.36) | 6.80 (99 reads unplaced) |
| roseM1_R1 | 0.797 | 0.760 | 1.393 (+0.60, p=3e-4) | 0.947 (+0.15, p=.06) | 1.437 (+0.64) | 4.71 (94 unplaced) |

### Offline selection replay

From `results/selection_replay.txt`: the share of WITCH's full-ensemble top-10 weight that each rule's selected
HMMs cover.

- **beam-3:** 0.98 (16S_b1000), 0.99 (RNASim), 0.99 / 0.97 (ROSE), using 36–52 scores per read.
- **BLAST paths:** cover only 0.55–0.70 on ROSE, because BLAST hits there are poor.
- **Fallback:** scoring all HMMs when beam-3's best bit-score is under 20 gives 1.000 coverage on ROSE, with
  25–28% of reads falling back (`results/fallback_sim.txt`). The same trigger catches nothing on 16S or RNASim,
  where the few misses have high bit-scores. A better uncertainty signal is still needed.

### Where the time goes

**HMM scoring is the bottleneck.** It is 65–86% of WITCH's wall time, and its share grows with backbone size
(75% at 1k on 16S, 86% at 5k).

**Selection cuts the scoring step by 4–20x and scales like log(#HMMs).** Beam-3 scores 36–52 HMMs per read
whatever the ensemble size (199–1387). At the 5k backbone it took 60 s of selection against 1213 s.

**After selection, WITCH's own alignment stage is the floor.** It costs 50–90 s per 1000 reads for k = 10, and
still 22–38 s with k = 1. About half of that is per-read Python overhead building constraint backbones, not
hmmalign. Adaptive k cuts mean k to 2–5 for +0.03–0.3 pp, but saves only 10–20% of the stage.

The outcome against the < 20% criterion therefore depends on backbone size:

| backbone | lite runtime vs WITCH | meets < 20%? |
|---|---|---|
| 16S_b5000 | 15–17% (10–11% excluding decomposition) | yes |
| 1k backbones | 27–37% | no |
| tiny ROSE instances | about 50% (fixed costs dominate) | no |

### The uncomfortable baseline result

**Sensitive BLASTN alone matches or beats WITCH on 16S.**
- 16S_b1000: +0.03 pp.
- 16S_b5000: **−0.81 pp** (0.48% vs 1.29%; 448 reads better, 57 worse), at about 1% of the runtime.
- With a dense real 16S backbone the nearest reference sequence is very close, and pairwise alignment to it
  beats the HMM ensemble.
- TIPP3's own default (megablast) is also better than WITCH on 16S_b5000 (0.94%).

**BLAST collapses on divergent data**, where WITCH and the lite variants win by a wide margin:

| instance | sensitive BLASTN SPFN | megablast SPFN |
|---|---|---|
| RNASim | 11% | 79% |
| ROSE | 61–68% | 89% |

So WITCH-lite only matters in the regime where reads are far from every reference sequence. These pilots did
not test whether TIPP3's 38 protein-coding marker genes are in that regime.

## Go/kill

The criterion was: within about 1 SPFN point of WITCH, at under 20% of its runtime.

| case | accuracy | runtime | verdict |
|---|---|---|---|
| 16S, 5k backbone | met (+0.01 to +0.03 pp) | met (15–17%) | GO |
| 16S, 1k backbone | met | 27–37% | not met on time |
| RNASim, 1k backbone | met (+0.01 to +0.08 pp) | 27–37% | not met on time |
| ROSE, high divergence | +0.3 to +1.5 pp (beam-3/4) | about 50% | borderline |

Other findings:
- UPP2's own search transfers poorly to WITCH: hier is +0.2 to +6.4 pp and EarlyStop +1 to +15 pp, because
  WITCH needs the top 10 HMMs, not just the best one.
- BLAST-guided selection is good when BLAST is good, which is exactly when WITCH isn't needed.

## Verdict: **unclear, leaning not promising as a "new method" project; promising only as a narrow engineering + benchmark project.**

**Why not stronger:**
- The core idea (hierarchical descent) is UPP2's.
- The speed-up is capped by WITCH's alignment stage unless that stage is also re-engineered.
- On 16S, the strongest baseline for TIPP is BLAST, not WITCH, so "near-WITCH" is not the right bar there.

**What would make it worthwhile:**
- Showing the WITCH > BLAST gap on TIPP3's real marker genes (protein-coding, larger backbones).
- In that setting, delivering an "eHMM beam search + uncertainty fallback + cached decomposition + faster
  stage" that gets WITCH to at most 15% of runtime at at most 0.2 pp SPFN loss.

### If pursued: 4 weeks

1. **Week 1:** download the TIPP3 reference package (Illinois Data Bank IDB-4931852) and pick 3–4 marker genes.
   Simulate reads from held-out reference sequences. Measure WITCH vs megablast vs blastn per marker. **This is
   the go/kill gate: is WITCH clearly better than BLAST on real markers?**
2. **Week 2:** beam-b selection with an uncertainty fallback (child-score margin, or a per-read entropy of the
   weights). Cache the decomposition. Run it inside TIPP3 (the BLAST hit is free there).
3. **Week 3:** cut the alignment-stage floor: cache per-HMM constraint backbones, and use adaptive k. Measure
   end-to-end TIPP3 profiling accuracy (Hellinger distance) and runtime against TIPP3 and TIPP3-fast.
4. **Week 4:** scaling curves (backbone 1k → 50k, ensemble size), write-up.

### Main risks

- WITCH may not beat BLAST on real markers. That kills the motivation, and the 16S evidence here points that way.
- The gains are engineering. The novelty over UPP2 is thin (beam + top-k + fallback).
- The remaining time sits in WITCH's alignment stage, which is somebody else's code.
- The TIPP3 reference package is large; download and per-marker WITCH baselines take days of CPU on 4 cores.

## Caveats

- Five instances, a single read length (150 bp), no sequencing errors.
- ROSE instances are small (300 reads), so fixed costs dominate their runtime ratio.
- The hybrid's selection time is estimated, not measured.
- Placement delta error was computed against a FastTree backbone tree; its reference is the true tree.
