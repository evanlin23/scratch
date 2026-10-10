# epangdown: what fixing EPA-ng's large-tree bug buys downstream (BSCAMPP)

Pilot, one overnight session (4 cores, 15 GB RAM). Builds on the root cause found on branch
`claude/cs581-epang` (cs581/epang/): with `--rate-scalers auto` and more than 2,000 tips,
EPA-ng 0.3.8 turns on per-rate scalers, and then
**bug 2** (`shift_partition_focus` moves the scale buffers by `offset` instead of
`offset*rate_cats`) misaligns them in the pre-masked thorough phase (accuracy), and
**bug 1** (`attributes = PLL_ATTRIB_RATE_SCALERS;` drops the SIMD bits) runs large trees on
scalar kernels (speed). Patch: `code/epa-ng-fix.patch` (copied from that branch).

RESULTS_PLACEHOLDER
