# epangdown: what fixing EPA-ng's large-tree bug buys downstream (BSCAMPP, PICRUSt2)

Pilot, one overnight session (4 cores, 15 GB RAM, no GPU), 2026-10-10. Follows the root cause
found on branch `claude/cs581-epang` (cs581/epang/). With `--rate-scalers auto` (the default) and
**more than 2,000 tips**, EPA-ng 0.3.8 turns on per-rate scalers, and then:
- **bug 2**: `shift_partition_focus` moves the scale buffers by `offset` sites instead of
  `offset*rate_cats`, so they are misaligned in the pre-masked thorough phase (accuracy);
- **bug 1**: `make_partition` sets `attributes = PLL_ATTRIB_RATE_SCALERS;`, which drops the SIMD
  bits, so large trees run on scalar kernels (speed).

Before starting I checked that branch's latest results (commit ec9e39c, RNASim 10K nested
subtrees). The fix removes the jump there: for fragments at k=3000, stock gives mean delta 1.74
and fixed gives 0.78. So I went ahead.

## Bottom line

0. **(Added after the orchestrator's follow-up.) Whole-tree EPA-ng, BSCAMPP paper Exp. 5.**
   EPA-ng on the whole backbone does not fit for 77K/49K leaves, so I used a random
   10,000-leaf nt78 sub-backbone and an 8,000-leaf RNASim sub-backbone, with the same 1,000
   fragments. **Stock whole-tree EPA-ng is much worse than BSCAMPP, and fixed whole-tree EPA-ng
   is as good as BSCAMPP:**
   - nt78 10K: stock whole-tree 2.55, fixed whole-tree 0.78, BSCAMPP(e) stock b=2000 0.79.
     Fixed vs stock whole: 565/381/54, p=2e-82.
   - RNASim 8K: stock whole-tree 1.72, fixed whole-tree 0.78, BSCAMPP b=2000 0.81. Fixed vs
     stock: 470/436/94, p=6e-52.

   So the paper's "EPA-ng on the full tree is much more error-prone" (Exp. 5) is, on these
   data, the EPA-ng bug and not a weakness of whole-tree placement. With the fix, whole-tree
   EPA-ng is also fast: 24 s vs 32-39 s stock, and BSCAMPP takes 43-69 s. Its limit is memory
   (12-13 GB at 8-10K leaves).
1. **The fix explains BSCAMPP's open problem end to end.** On three benchmarks with 26K-77K
   leaves, BSCAMPP with stock EPA-ng gets worse delta error when the subtree size goes from 2,000
   to 5,000 (and 10,000): 2.8x on nt78 (1.72 → 4.77; 5.09 at 10,000), 2.7x on RNASim 50K
   (0.53 → 1.43), and +11% on the noisy 16S.B.ALL (25.3 → 28.0). With the patched EPA-ng, delta error at 5,000 and 10,000 is
   the same as at 2,000. **All of the accuracy effect comes from bug 2.** The bug-2-only build
   matches the full fix on 98% of queries. The bug-1-only build matches stock on 98.6%.
2. **The fix does not make BSCAMPP more accurate than the published default.** With the fix,
   larger subtrees are no longer harmful, but they are not measurably better either. Fixed
   b=5,000/10,000 vs stock b=2,000 differs by -0.02/-0.03 on nt78 (n.s.), -0.02 on RNASim 50K
   (p=0.09), and +1.0 to +1.6 on 16S.B.ALL (n.s., noisy). The default b=2,000 that the BSCAMPP
   paper chose sits just below the bug threshold, so the default pipeline was never affected.
3. **Speed: the fix is a real, free 1.6-3.3x speedup of EPA-ng on trees over 2,000 tips**
   (single thread, nt78 subtrees, same placements as `--rate-scalers off`). Both bugs cost time
   on fragments. Inside BSCAMPP at b=5,000 the EPA-ng time drops by 18-41%. At b=2,000, BSCAMPP
   runs entirely below the threshold and gains nothing.
4. **Downstream tool (PICRUSt2, 12,000-tip subset of its bacterial reference, 749 V4 ASVs):**
   2-2.7% of the ASVs change placement (two random reference subsets; a stock re-run changes
   0). On one subset, five ASVs are grossly misplaced by stock (NSTI 32 vs 1.3) and dropped by
   the NSTI filter, and one sample's predicted metagenome moves by 16%. On the other subset the
   worst sample moves by 3%. The median sample barely changes (≤0.2%). The real 26,868-tip
   reference needs at least 16 GB for EPA-ng and did not fit here.

**Verdict for a 4-week CS581 project built on this bug fix: unclear, but better than it first
looked.** Within BSCAMPP it is not promising as a "downstream accuracy" project (point 2). But
the whole-tree result (point 0) gives the project a publishable-style claim: a published
comparison (EPA-ng vs BSCAMPP/SCAMPP) is confounded by the bug, and fixed EPA-ng matches
BSCAMPP whenever it fits in memory. BSCAMPP's remaining value is then memory and scaling,
not accuracy. The diagnosis is a solid, novel result (explanation + fix +
attribution). But the payoff in BSCAMPP is "the larger-subtree setting stops being broken", not
"BSCAMPP gets more accurate". The default was already safe, and memory caps EPA-ng subtrees at
about 10K leaves on a laptop. The speed result is clean but modest. The most promising angle is
the tools that run EPA-ng *directly* on references with more than 2,000 tips and fragmentary
queries (PICRUSt2, TIPP3, PEWO-style pipelines), which have no b=2,000 escape hatch. But the
PICRUSt2 pilot shows a smaller effect than BSCAMPP: 2-3% of the placements change. All the V4
amplicons crop the alignment to the same window, which seems to limit the damage. Queries with
varying fragment positions (metagenomic reads, TIPP3) are the case to test (see weeks 1-4
below).

## 1. Methods

**Builds** (`code/build_epang.sh`). EPA-ng v0.3.8 (github.com/pierrebarbera/epa-ng tag v0.3.8,
commit cff6a47) is built once with cmake, Release `-O3`, gcc 13.3, OpenMP. The other variants are
rebuilt incrementally in the *same* build tree after `git apply`, so the compiler flags are
identical and only the patched translation units change:
`epa-ng-stock`, `epa-ng-fix` (both hunks of `code/epa-ng-fix.patch`), `epa-ng-bug2only`
(`pll_util.cpp` hunk: scaler shift times `rate_cats`), and `epa-ng-bug1only` (`file_io.cpp`
hunk: `|=`, keep the SIMD bits).

**BSCAMPP** (github.com/ewedell/BSCAMPP, v1.0.8, commit ea75e98, pip-installed). Its
`epang_path` points to `code/epa_sel.sh`, which runs the build named by `$EPA_BIN` and logs every
EPA-ng call (tips, seconds, return code). Settings: `-V 5 --threads 4 --cpus-per-job 4`, which
means one EPA-ng job at a time with 4 threads (two 2-thread jobs ran out of memory at b>=5,000).
Wall-clock is measured with `/usr/bin/time -v`.

**Data.** The BSCAMPP paper's RNASim-VS backbones are on Dryad (doi 10.5061/dryad.78nf7dq), and
its API now requires a bearer token, so they were not reachable. Instead I used:
- **nt78**: SCAMPP benchmark (Illinois Data Bank IDB-9257957; FastTree-2 16S-like ROSE
  simulation, 78,132 seqs, 1,287 sites). The true tree is known. The backbone tree is the
  provided FastTree topology with RAxML-NG branch lengths and model (`alignment.phylip.raxml.*`).
- **16S.B.ALL**: same release (CRW, 27,643 seqs, 6,857 sites masked to 1,603). The "true" tree
  is the published RAxML reference tree (as in SCAMPP/BSCAMPP), and the backbone uses the same
  topology, so delta error = placement error relative to that tree.
- **RNASim 50K**: Recursive-MAGUS release (IDB-1048258, `RNASim/50000_R_1/R0`, true alignment
  and true tree; 12,458 sites masked to 1,624). The backbone tree is FastTree-2 `-nt -gtr
  -gamma` on 49,000 sequences (56 min, 1 thread). A RAxML-NG `--evaluate` on all 49K tips does
  not fit in 15 GB, so the GTR+G parameters come from RAxML-NG `--evaluate` on a random
  5,000-leaf pruned subtree (`code/model_from_subtree.py`), and the branch lengths stay
  FastTree's.

On each dataset, 1,000 random sequences are removed from the backbone all at once (batch
placement as in BSCAMPP; the tree is pruned and kept unrooted). Columns with more than 95% gaps
among the backbone sequences are masked. The queries are fragments: a random start and length ~
N(10% of the ungapped length, sd 10 nt), as on the epang branch (`code/prep_scampp.py`,
`code/prep.py`). Mean fragment length after masking is 128 nt (nt78), 135 (16S) and 154
(RNASim).

**Delta error** uses the epang branch's scorer (`code/deltalib.py`, unchanged): FN(T+q, T*|B+q) −
FN(T, T*|B) on the full backbone. Subtree placements are mapped back to the backbone. I made a
memory-lean driver (`code/score_multi.py`) because the original per-query cache does not fit
for 77K leaves. Pairing is by query. The tables report mean delta, W/T/L (W = lower delta), and
the Wilcoxon signed-rank p (`code/analyze.py`).

**Noise floor.** BSCAMPP + EPA-ng is not exactly repeatable. Exactly tied placements (identical
sequences, zero-length branches) are broken differently from run to run, and this is not fixed
by `PYTHONHASHSEED`. Identical-code reruns at b=2,000 differ as follows:
- nt78: on 10-11 of 1,000 queries (mean ±0.02);
- RNASim: on 0 of 1,000;
- 16S.B.ALL: on about 600 of 1,000 (mean ±0.3-0.6, never significant). The differing queries
  all have best LWR ≈ 0.003 with identical LWR in both runs, so they are many-way ties.

16S is therefore only useful for large effects.

## 2. BSCAMPP accuracy and wall-clock (1,000 fragmentary queries each)

### nt78 (77,132-leaf backbone, simulated, true tree known)

| run | mean delta | median | % delta=0 | vs stock b=2000: diff | W/T/L | Wilcoxon p | BSCAMPP wall | EPA-ng time (sum) | max RSS |
|---|---|---|---|---|---|---|---|---|---|
| stock b=2000 (published default) | 1.720 | 1 | 37 | | | | 5:06 | 273 s | 2.4 GB |
| fix b=2000 | 1.721 | 1 | 37 | +0.001 | 10/981/9 | 0.81 | 5:36* | 284 s | 2.4 GB |
| **stock b=5000** | **4.769** | 3 | 16 | **+3.049** | 103/309/588 | 3e-79 | 9:16* | 516 s | 6.5 GB |
| fix b=5000 | 1.697 | 1 | 37 | -0.023 | 39/924/37 | 0.56 | 7:36* | 424 s | 6.4 GB |
| bug2only b=5000 | 1.686 | 1 | 37 | -0.034 | 38/928/34 | 0.37 | 9:00* | 507 s | 6.4 GB |
| bug1only b=5000 | 4.756 | 3 | 16 | +3.036 | 103/310/587 | 4e-79 | 6:53* | 381 s | 6.5 GB |
| **stock b=10000** | **5.089** | 3 | 16 | **+3.369** | 101/288/611 | 3e-84 | 10:38 | 610 s | 13.1 GB |
| fix b=10000 | 1.690 | 1 | 38 | -0.030 | 46/913/41 | 0.32 | 8:15 | 461 s | 13.0 GB |

Fixed vs stock at the same size: b=5000 -3.07 (W/T/L 584/325/91, p=1e-80); b=10000 -3.40
(620/286/94, p=3e-87). bug2only vs fix: identical delta on 98.3% of queries (n.s.). bug1only vs
stock: identical on 98.6% (n.s.).
\* = run while FastTree (1 thread) was building the RNASim tree on the same machine. The wall
times are comparable within the starred group; the per-call EPA-ng sums are more reliable.

### RNASim 50K (49,000-leaf backbone, simulated, true tree known; PYTHONHASHSEED=0)

| run | mean delta | % delta=0 | vs stock b=2000: diff | W/T/L | p | BSCAMPP wall | EPA-ng time | max RSS |
|---|---|---|---|---|---|---|---|---|
| stock b=2000 | 0.530 | 68 | | | | 3:06 | 158 s | 2.8 GB |
| fix b=2000 | 0.530 | 68 | 0.000 | 0/1000/0 | 1 | 2:57 | 152 s | 2.8 GB |
| **stock b=5000** | **1.434** | 41 | **+0.904** | 78/488/434 | 4e-53 | 6:30 | 366 s | 7.4 GB |
| fix b=5000 | 0.507 | 68 | -0.023 | 22/963/15 | 0.089 | 5:11 | 287 s | 7.4 GB |

Fixed vs stock at b=5000: -0.927, 439/487/74, p=1e-55. b=10,000 does not fit (about 16 GB for
EPA-ng at 1,624 sites).

### 16S.B.ALL (26,643-leaf backbone, biological; reference tree as truth)

| run | mean delta | % delta=0 | vs stock b=2000: diff | W/T/L | p | BSCAMPP wall | EPA-ng time |
|---|---|---|---|---|---|---|---|
| stock b=2000 | 25.27 | 8 | | | | 2:44 | 158 s |
| stock b=2000 (rerun) | 24.95 | 8 | -0.32 | 317/409/274 | 0.25 | 2:20 | |
| stock b=2000 (seed 0) | 25.25 | 8 | -0.02 | 294/412/294 | 0.98 | 2:10 | |
| fix b=2000 / seed 0 | 25.44 / 24.67 | 8 | +0.17 / -0.60 | ~300/410/290 | 0.56 / 0.49 | 2:01 | 115 s |
| **stock b=5000 / seed 0** | **28.03 / 27.93** | 4 | **+2.75 / +2.66** | 389/125/486 | 6e-4 / 4e-4 | 4:12 / 3:33 | 247 s |
| fix b=5000 / seed 0 | 26.85 / 26.23 | 8 | +1.58 / +0.96 | 323/355/322 | 0.45 / 0.43 | 2:33 / 2:08 | 146 s |
| bug2only b=5000 | 25.33 | 8 | +0.06 | 335/351/314 | 0.58 | 3:15 | 191 s |
| bug1only b=5000 | 28.20 | 4 | +2.93 | 391/115/494 | 3e-4 | 2:24 | 140 s |

The unseeded 16S runs (first five rows without "seed 0", plus bug*only) ran while FastTree was running, so their times are noisy. For example, the EPA-ng sums at b=2,000 (158 s vs 115 s) come from identical code. Fixed vs stock at b=5000: -1.18 (p=0.009) and -1.69 (seed 0, p=0.003). The large absolute delta
error comes from biological 16S fragments with many identical or near-identical sequences in a
27K-leaf tree. Ties make this dataset noisy (see noise floor). Note "% delta=0": stock at 5000
halves the exactly-correct placements on all three datasets (37→16%, 68→41%, 8→4%).

### BSCAMPP(p) (pplacer) vs BSCAMPP(e) at matched subtree sizes (nt78, 1,000 fragments)

BSCAMPP's bundled pplacer v1.1.alpha19 with taxit refpkgs, on the same backbone tree. Its model
comes from the dataset's `RAxML_info.REF` (RAxML 7 GTR+Γ), while EPA-ng uses the RAxML-NG
model. I ran one pplacer job at a time with 4 threads (BSCAMPP passes `-j` threads).

| run | mean delta | % delta=0 | BSCAMPP wall | max RSS |
|---|---|---|---|---|
| BSCAMPP(p) b=2000 | 1.705 | 41 | 13:30 | 2.0 GB |
| BSCAMPP(p) b=5000 | 1.633 | 42 | 23:31 | 5.1 GB |
| BSCAMPP(e) stock b=2000 | 1.720 | 37 | 5:06 | 2.4 GB |
| BSCAMPP(e) fixed b=2000 | 1.721 | 37 | 5:36* | 2.4 GB |
| BSCAMPP(e) stock b=5000 | 4.769 | 16 | 9:16* | 6.5 GB |
| BSCAMPP(e) fixed b=5000 | 1.697 | 37 | 7:36* | 6.4 GB |
| BSCAMPP(e) fixed b=10000 | 1.690 | 38 | 8:15 | 13.0 GB |

Paired comparisons:
- (p) vs fixed (e) at b=2000: -0.016, W/T/L 105/818/77, p=0.27.
- (p) vs fixed (e) at b=5000: -0.064, 112/821/67, p=0.006.
- (p) b=5000 vs (p) b=2000: -0.072, p=0.04.
- Stock (e) at b=5000 is 3x worse than (p) at b=5000.

**pplacer is about as accurate as fixed EPA-ng: slightly better at b=5000, and with more exact
placements (41-42% vs 37%). It is 2.6-3.1x slower and uses about 20% less memory.** Unlike stock
EPA-ng, pplacer shows no jump at b>2000, which fits a bug in EPA-ng rather than a property of
large subtrees. For "TIPP3 with BSCAMPP(p) instead of BSCAMPP(e)": at b=2000 the two are tied
on accuracy, and (e) is much faster. Above 2,000, (e) needs the fix to be usable.

### Memory is the practical cap on subtree size
EPA-ng RSS is about 1.3 GB per 1,000 tips at 1,286 sites (5,000 → 6.4 GB, 10,000 → 13 GB) and
grows with the number of sites. On a 15-16 GB machine, BSCAMPP b=10,000 only fits for nt78 when
nothing else runs. On 16S and RNASim it does not fit. So even with the fix, "use larger
subtrees" is mostly unavailable to a laptop user, and it would not improve accuracy anyway
(above).

## 3. Speed: stock vs fixed EPA-ng standalone (`code/speedtest.sh`)

nt78 backbone, k-leaf subtree (first k leaves in preorder), 200 queries, `-T 1`, each build run
alone. Wall-clock (`/usr/bin/time`):

| k | queries | stock | bug1only (SIMD back) | bug2only | **fix** | stock `--rate-scalers off` | speed-up stock→fix |
|---|---|---|---|---|---|---|---|
| 5,000 | fragments | 126.9 s | 53.1 s | 89.0 s | **38.6 s** | 34.2 s | **3.3x** |
| 5,000 | full-length | 23.2 s | 13.3 s | 23.1 s | **13.5 s** | 12.0 s | **1.7x** |
| 10,000 | fragments | 133.1 s | 64.6 s | 87.7 s | **46.0 s** | 37.2 s | **2.9x** |
| 10,000 | full-length | 45.0 s | 28.2 s | 44.2 s | **28.5 s** | 22.6 s | **1.6x** |

Max RSS: 5.43 GB (k=5,000) and 10.85 GB (k=10,000) for all four builds; 5% less with rate
scalers off.

Best-edge agreement (`results/speed/placement_agreement.txt`): for full-length queries all five
runs agree on 100% of queries. That is expected: pre-masking crops nothing from a full-length
query, so the misaligned shift is 0. For fragments, fix = bug2only = rate-scalers-off on 100%,
and stock = bug1only on 100%. Stock and fix agree on only 0-6%. (These 200 queries mostly do not
belong to the subtree, which magnifies the disagreement. In BSCAMPP, where the subtree is chosen
around the query, 32% of the queries keep the same delta.) So bug 1 is purely a speed bug, bug 2
is an accuracy bug that also costs time (the misaligned scalers make the thorough phase slower),
and the fixed build gives the same answers as `--rate-scalers off` at 1.1-1.3x its run time.

## 4. Is it a CS581-sized project? Weeks 1-4

What is already done (both pilots): root cause, a 2-line fix, attribution of accuracy to bug 2
and speed to bug 1 and bug 2, and the end-to-end BSCAMPP picture on three datasets. What is left
is (a) a stronger downstream story and (b) the obvious report items. A realistic plan:

- **Week 1:** Make the fix upstream-ready. Check the patch against the EPA-ng master branch
  (0.3.8 is the last tag). Test it with `--rate-scalers on` on small trees (the epang branch
  saw a `fix_rson` anomaly at k=500 that needs explaining), and check the binary-dump and
  `--no-pre-mask` paths. Open the issue and PR upstream.
- **Week 2:** The direct-EPA-ng users. Re-run PICRUSt2 (full 26,868-tip bacterial reference;
  needs a 32-64 GB machine, e.g. a campus cluster node) and TIPP3 or a PEWO-style
  pruning benchmark on references with more than 2,000 tips and short amplicon/fragment
  queries. Measure how many placements, NSTI values and predicted functions change.
- **Week 3:** BSCAMPP/SCAMPP sweep with replicates (several RNASim 50K/200K replicates once the
  Dryad backbones are available; b = 1,000-10,000 on a big-memory node). Show that with the fix
  the delta-vs-b curve is flat, and report the optimal b for time given accuracy.
- **Week 4:** Write-up: the bug, why the threshold is exactly 2,000 tips, figures (delta vs b,
  stock vs fixed; speed table), and a recommendation (`--rate-scalers off` as a workaround
  for older binaries).

Risk: the BSCAMPP accuracy payoff is null with the fix (fixed b>2,000 ≈ stock b=2,000). A
project whose headline needs "fixed BSCAMPP beats published BSCAMPP" will not get it. The
headline that holds is "a silent EPA-ng bug degrades every placement of a fragment into a
>2,000-tip reference, and here is what that does to tools X and Y".

## 5. PICRUSt2 (item 4)

**Setup.** PICRUSt2 2.6.3 (bioconda). Its post-link download of the reference files is blocked
here, so I took the reference from the git tag v2.6.2. `place_seqs.py` calls EPA-ng with the
defaults plus `--filter-acc-lwr 0.99 --filter-max 100`, so rate scalers are on above 2,000 tips.
On the real 26,868-tip bacterial reference (1,578 sites), EPA-ng is killed at an 8 GB address
limit within 6 s. Scaling from the 12,000-tip run below (7.3 GB), it needs at least 16 GB, so
it does not fit here. I therefore used a **random 12,000-tip subset of the bacterial reference**
(`code/picrust2_subref.py`). The tree is pruned, the alignment subset, the HMM rebuilt with
`hmmbuild`, the model kept, and the default 16S/KO/EC trait tables subset to those tips.
Steps (`code/run_picrust2_subref.sh`): `place_seqs.py` → `hsp.py` (16S + NSTI, KO, EC; max
parsimony) → `metagenome_pipeline.py`. One gotcha: hsp with r-castor 1.8.7 fails ("invalid
'ncol'"), so I downgraded it to 1.7.11.

Queries: the QIIME 2 "Moving Pictures" tutorial (770 ASVs, 120 bp V4, 34 samples;
docs.qiime2.org 2024.10 `rep-seqs.qza`/`table.qza`). 749 ASVs pass the alignment filter. The
same ASVs and reference are placed with stock and with fixed EPA-ng (`-T 4`).

| | stock | fixed |
|---|---|---|
| EPA-ng placement step, wall / max RSS | 23.7 s / 7.3 GB | 18.4 s / 7.4 GB |
| ASVs above the NSTI cut-off (2.0), dropped downstream | 9 | 4 |
| mean NSTI over 749 ASVs | 0.341 | 0.140 |

Comparison (`results/picrust2/compare.txt`):
- **Best placement edge changes for 20 of 749 ASVs (2.7%)** and NSTI for 22. Of the 22, 5 ASVs
  get an NSTI of 31.9-32.9 with stock (a blown-up pendant length, the typical sign of the
  scaler misalignment) vs 1.27-1.35 with the fix. These are the extra ASVs dropped by the
  NSTI ≤ 2 filter. For the other 17 the NSTI moves both ways: the fix gives lower NSTI for 7
  of 22.
- **Predicted gene content** (KO and EC) changes for 16 ASVs.
- **Predicted metagenome** (34 samples): per-sample relative L1 difference median 0.18%, max
  16% (KO) and 15.5% (EC). The worst sample is L3S360, where the changed ASVs make up 10.8% of
  reads. Per-sample Spearman between stock and fixed is at least 0.955. 3,681 of 7,840 KOs
  change in at least one sample.

Controls: (i) stock re-run on the same subset gives **0 changes**, so PICRUSt2 + EPA-ng is
deterministic here and all 20 changes are due to the fix. (ii) A second random 12,000-tip subset
(seed 2): best edge changes for **15 of 749 (2.0%)**, NSTI for 24 and gene content for 14 ASVs.
There are no NSTI blow-ups this time (mean NSTI 0.124 stock vs 0.137 fixed, 5 vs 5 ASVs over
the cut-off). The predicted metagenome changes by 0% at the median and up to 3.0% (KO) / 3.3%
(EC) in the worst sample; Spearman ≥ 0.974.

Reading: for full-region 16S amplicons that all cover the same V4 window, the pre-mask range is
the same for every query, and only a few percent of placements change. Those few sometimes include
grossly wrong placements (NSTI about 32 on subset 1, none on subset 2), which shift the
predictions of the samples where those ASVs are abundant by up to about 15% (3% on subset 2). That is a real but modest downstream effect on this
dataset. It needs the full 26,868-tip reference (big-memory machine), more studies and other
amplicon regions before one can claim more. Not done: TIPP3-fast (time).

## Files
- `code/`: builds (`build_epang.sh`, `epa-ng-fix.patch`), data prep (`prep_scampp.py`,
  `prep.py`, `prep_rnasim.sh`, `model_from_subtree.py`), runs (`run_bscampp.sh`, `epa_sel.sh`,
  `queue*.sh`, `speedtest.sh`, `make_subtree.py`, `picrust2_subref.py`,
  `run_picrust2_subref.sh`), scoring (`deltalib.py`, `score.py` from cs581/epang; `score_multi.py`,
  `analyze.py`, `speed_compare.py`, `compare_picrust2.py`).
- `results/<dataset>/`: per-query delta (`scores.tsv.gz`), `analysis.md`, BSCAMPP `times.tsv`
  (variant, b, label, wall, max RSS KB), and per-EPA-ng-call logs (`*.epa_calls.tsv`: tips,
  seconds, rc). The rows with b=10000 and an empty time in nt78/16S `times.tsv` are runs killed
  for running out of memory while FastTree was running; nt78 b=10000 was rerun alone.
- `results/speed/`: standalone timing and placement agreement.
- `results/nt78_10k/`, `results/rna_8k/`: whole-tree EPA-ng (stock vs fixed) vs BSCAMPP on 10K/8K
  sub-backbones (`code/make_subset_dataset.py`, `code/run_whole.sh`).
- `results/nt78/scores_p.tsv.gz`, `analysis_pplacer.md`: BSCAMPP(p) comparison.
