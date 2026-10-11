# Pre-registration: does the GCM evidence recipe improve ML tree accuracy?

Written and pushed before any tree in this study was computed (and before any alignment of the new
replicates R5-R8 was scored).

## Question

Does `wsoft0.03:linsi&fftns2#es4` (the gcmgen recipe; L-INS-i backbone pairs unconfirmed by FFT-NS-2 at
weight 0.03, then cross-subset GCM edges with support < 4 of 10 backbones deleted) improve maximum-likelihood
tree accuracy over MAGUS on simulated proteins? Same question for `linsi#es3` and for the hard filter
`linsi&fftns2-op3`.

## Data

AliSim (IQ-TREE 3.1.4) LG+G4, 1,000 taxa, root length 300, indels 0.05/0.05 POW(1.7, 40); Yule-Harding tree
`-rlen 0.001 MEAN 0.8`, MEAN = 0.06 (SIMMOD) / 0.10 (SIMHIGH); tree seed 100r+7, sequence seed 100r+13,
for r = 1..8. R1-R4 = protbench's replicates (check: SIMMOD_R1 must have 7,348 columns); R5-R8 are new.
16 datasets. Run order: SIMHIGH R1-R8 first, then SIMMOD R1-R8. If the time budget runs out, the
datasets not reached are reported as not run; nothing is dropped after its trees are seen.
Optional extras (only if time allows, reported separately and not part of the primary test):
a second MAGUS draw on SIMHIGH; a harder level.

## Methods (paired, merge-only on one MAGUS draw per dataset)

One MAGUS run per dataset with the paper's flags (25 subsets, 10 backbones x 200, MCL, minclusters).
Alignments, all on MAGUS's own subsets and backbone sequence sets:
`true` (true alignment), `linsi` (MAGUS's merge, the baseline), `wsoft0.03:linsi&fftns2#es4` (recipe),
`linsi#es3`, `linsi&fftns2-op3` (hard filter).

Trees: FastTree 2 `-lg -gamma` (single-threaded) on every alignment; normalized RF (FN, FP) to the true
tree (DendroPy, unrooted). If time allows: IQ-TREE 3 `-m LG+G4 --fast` on the SIMHIGH alignments (at least
true, MAGUS, recipe).

## Endpoints

- **Primary:** FastTree nRF (%) to the true tree, recipe minus MAGUS, paired over all simulated datasets
  run (target 16), two-sided Wilcoxon signed-rank (exact; zero differences dropped, i.e. scipy default
  `zero_method="wilcox"`). Significance at alpha = 0.05.
- **Secondary:** the same for `linsi#es3` and for `linsi&fftns2-op3`; the recipe restricted to SIMHIGH.
  Also reported: W/T/L with a 0.1 RF-point tie band; SP error (SPFN, SPFP) deltas; room =
  RF(MAGUS) - RF(true); IQ-TREE deltas if run.
- Verdict rule: "the recipe improves trees" only if the primary test has p < 0.05 with a negative mean.
  Otherwise the why-not analysis below is reported (correlation of tree change with SPFN / SPFP change;
  whether changed residue pairs fall in variable or near-constant columns of the true alignment; size of
  the effect vs FastTree / MAGUS-draw noise).
