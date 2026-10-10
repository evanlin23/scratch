---
title: "CS581 Project Proposal: Why Phylogenetic Placement Breaks on Large Trees, and Fixing It"
author: "Evan Lin · [NetID] · October 2026 (draft, pilot numbers provisional)"
---

**Problem.** Phylogenetic placement adds query sequences (often short reads) to a fixed reference tree
and underlies taxonomic profiling (TIPP3), functional prediction (PICRUSt2) and tree updating. Maximum
likelihood placers such as EPA-ng and pplacer are accurate but do not scale to very large reference trees,
so the Warnow lab's SCAMPP and BSCAMPP place each query into a small subtree around its closest reference
sequences. Wedell, Shen & Warnow (BSCAMPP, IEEE/ACM TCBB 2025) report an unexplained anomaly that limits
this design: *"When the subtree size exceeds 2,000 leaves, BSCAMPP(e) had over twice the delta error than
it did for 1,000- and 2,000-leaf subtrees; therefore, we set the default subtree size to 2,000"*, and
*"EPA-ng may have some numeric issues when placing into very large trees ... Further research is needed to
understand whether this explanation is correct."* The same paper finds EPA-ng on the full 180K-leaf RNASim
tree much less accurate than every subtree-based pipeline. If the anomaly is an artifact, the default
subtree size, and the conclusions drawn from it, should be revisited.

**Pilot finding (done).** The anomaly is a software bug, not a property of large trees. EPA-ng
(v0.3.3-v0.3.9-pre, including the bioconda v0.3.8 that BSCAMPP bundles) switches on per-rate likelihood
scalers when the tree has more than 2,000 tips (`large_tree()` is hard-coded as `tip_nodes > 2000`). Its
premasking code (`shift_partition_focus`) then shifts the scaler buffers by `offset` sites instead of
`offset × rate_cats`, so for any query that does not start at the first alignment column (fragments,
amplicons) the likelihoods of the thorough placement phase use the wrong scale factors (logL off by
thousands of units, falsely confident LWR = 1.0). Evidence on RNASim (true query alignments, fragments of
~10% length, paired queries): stock EPA-ng's mean delta error is 0.82 at 2,000 leaves and 1.74 at 3,000;
forcing scalers on at 500-2,000 leaves reproduces the jump; disabling premasking removes it; full-length
queries are unaffected, as the mechanism predicts. A 5-line patch removes the jump (0.78 / 0.79 / 0.73 at
3k / 5k / 9k leaves) and reproduces scaler-free placements exactly. A second hunk restores the SIMD kernels
that the same code path silently drops (speed only). On two large benchmarks the patched EPA-ng is flat in subtree size where stock EPA-ng degrades
2.7-3× (77K-leaf nucleotide set: stock 1.72 / 4.77 / 5.09 at 2k / 5k / 10k leaves vs patched 1.72 / 1.70 /
1.69, p = 1e-80 at 5k; RNASim 50K: 1.43 vs 0.51 at 5k, p = 1e-55), and on fragments it is 2.9-3.3× faster
standalone. Larger subtrees do *not* make BSCAMPP more accurate than its default 2,000 on these sets (ties;
a small gain only on RNASim 10K, −0.125, p = 0.005): the default happens to sit just under the bug's
trigger. PICRUSt2 2.6.3 is on the affected path; on a 12K-tip subset of its reference, 2.7% of V4 ASVs
change placement, a few by a lot (nearest-sequenced-taxon index ~32 vs ~1.3), shifting one sample's
predicted metagenome by 16% (median 0.2%). Most importantly, the BSCAMPP paper's conclusion that whole-tree EPA-ng is "much more" error-prone
than the subtree pipelines is confounded by the bug: on 8-10K-leaf backbones, stock whole-tree EPA-ng has
mean delta 1.72-2.55, the patched one 0.78, the same as BSCAMPP(e) at its default (0.79-0.81) and faster
(24 s vs 43-69 s), with memory (12-13 GB at 10K leaves) as the remaining reason to decompose. Against
the instructor-suggested comparison, BSCAMPP(p) (pplacer) is about as accurate as BSCAMPP(e) with the
patched EPA-ng (−0.016, n.s., at 2,000; −0.064, p = 0.006, at 5,000) but 2.6-3.1× slower. No issue,
release note or paper reports the bug.

**Research questions.**
(1) *Characterization.* Which queries, data and settings are affected: fragment length and start offset,
number of rate categories, DNA vs protein, EPA-ng heuristics (`--dyn-heur`, `--baseball-heur`,
`--no-heur`, where the bug is much worse), and how the error depends on tree size above 2,000 tips.
(2) *Re-tuning scalable placement.* With a correct EPA-ng, what subtree size should BSCAMPP and SCAMPP use (pilot: larger subtrees tie the default on two large sets, so the question is the accuracy/time/memory trade-off and when larger subtrees help)?
Accuracy-runtime-memory trade-off on the BSCAMPP benchmarks (RNASim up to 200K leaves, plus a biological
dataset), compared with BSCAMPP(pplacer), APPLES-2 and EPA-ng on the whole tree; does whole-tree EPA-ng
remain worse once fixed (pilot: no, it ties BSCAMPP on 8-10K-leaf trees; the question becomes memory and scale)?
(3) *Downstream impact.* How many placements and taxonomic assignments change in TIPP3-fast (which uses
BSCAMPP) and in PICRUSt2's placement step (EPA-ng on a ~26,900-tip default reference tree), and does that
change profiling accuracy on simulated communities?
(4) *Speed.* What restoring SIMD kernels buys on large trees.

**Approach and evaluation.** Build stock and patched EPA-ng from the same source (each hunk separately and
together). Data: RNASim (Illinois Data Bank) with the BSCAMPP protocol (masked alignment, fragmentary
queries, true and estimated query alignments), plus a 16S reference set for TIPP3/PICRUSt2. Metrics:
delta error (as in SCAMPP/BSCAMPP), placement edge distance, LWR calibration, taxonomic profiling error
(TIPP3's metrics), wall-clock and peak memory. All comparisons paired by query (Wilcoxon, W/T/L), several
replicates. Deliverables: the patch with a minimal reproducing test, a report to the EPA-ng maintainers,
and recommended BSCAMPP/TIPP3 settings.

**Timeline.** Week 1: reproduce on the full BSCAMPP RNASim benchmark; finish characterization (RQ1).
Week 2: subtree-size sweep for BSCAMPP and SCAMPP with the fix (RQ2), memory profiling. Week 3: TIPP3-fast
and PICRUSt2 downstream study (RQ3), speed (RQ4). Week 4: biological dataset, write-up, upstream report.

**Risks.** (a) The core diagnosis is already done, and the pilot shows the downstream accuracy gains are
small (BSCAMPP ties its default; PICRUSt2 changes 2-3% of placements), so the contribution is the diagnosis,
correctness on large trees, the 3× speed-up, and a careful characterization rather than a big accuracy gain. (b) Memory grows with
subtree size (12.9 GB at 9,000 leaves on 4 cores), which may cap the sweep; mitigated by fewer parallel
jobs. (c) The maintainers might fix it independently; the evaluation stands regardless. (d) Downstream
effects might be small if most queries are full length; fragments and amplicons are the common case in
metagenomics, which is the setting we target.

**Course connection.** Phylogenetic placement is part of the course: the divide-and-conquer lecture's
revised GTM pipeline adds sequences "using phylogenetic placement method BATCH-SCAMPP (with EPA-ng)", the
first-day slides cite pplacer-XR and pplacerDC, and the instructor's project list includes TIPP3/BSCAMPP
projects, verbatim "Evaluate TIPP3 with BSCAMPP(p) instead of BSCAMPP(e)". This project explains why
BSCAMPP(e) was limited to 2,000-leaf subtrees (the trigger is "more than 2,000 tips", so the default itself
is just below it) and re-evaluates BSCAMPP(e) vs BSCAMPP(p) once EPA-ng is correct, using the
maximum-likelihood ideas from lecture (site likelihoods, rate categories, numerical scaling).
