# gcmvote helper h13: RNASim 10K R0

AI-assisted (Claude), exploration code for CS581 project.

Data: RNASim 10,000-sequence subset, replicate R0, from the MAGUS paper's datasets (Illinois Data Bank,
doi:10.13012/B2IDB-2643961_V1, `Datasets.zip` → `RNASim/10000/R0/true_align.txt`, `true_tree.tre`), downloaded by
`cs581/code/setup.sh`. Machine: 4 cores, ~15 GB RAM. Code: `cs581/gcmvote/code_h13/` (magus_h13.py = the MAGUS draw
and replicate; pipeline.py = baselines, vote variants, trees; memrun.py = wall/CPU/peak memory).

## MAGUS non-recursive (--recurse false; the paper's 10K runs recurse)

The MAGUS draw used here is **not** the paper's configuration for 10K sequences. What happened:

1. First attempt: the paper's flags (`--maxsubsetsize 0 --maxnumsubsets 25 --decompstrategy pastastyle
   --decompskeletonsize 300 --graphbuildmethod mafft --graphclustermethod mcl --graphtracemethod minclusters
   --graphtraceoptimize false -r 10 -m 200 -f 4`, `-np 4`, vectorized graph builder), started 02:04 UTC.
   The top-level decomposition into 25 subsets took 1449 s (done 02:28). All 25 subsets have more than 200 sequences
   (295–550), above MAGUS's default recursion threshold (200, `recurse=true`), so MAGUS started a full recursive MAGUS
   on each subset. By 03:35, 8 of the 10 top-level backbones were done. The first recursive subset (subset_7) had
   finished its own decomposition (02:28–03:05) and 3 of its 10 backbones, and no other subset had started. Peak
   memory was about 2 GB. At about 1 h per recursive subset, the projection was about 25 h on 4 cores, so the run was
   stopped at 03:37 (logs: `/opt/work/h13/d0/recursive_attempt/`, not committed).
2. On the orchestrator's decision (option a), the run was restarted at 03:38 with the same flags plus
   `--recurse false`, in the same working directory. MAGUS reused the 25 top-level subset files ("Detected existing
   subset file" × 25, no re-decomposition) and the 8 finished top-level backbones (1–6, 9, 10). Backbones 7 and 8
   were unfinished, so MAGUS drew new sequence sets for them. Each subset is aligned directly by MAFFT L-INS-i
   (MAGUS's subset aligner). The reported MAGUS wall time covers only the restarted part. Add the 1449 s
   decomposition and the backbone time from attempt 1 for a full estimate.

All comparisons below are merge-only on identical subsets and backbones, so the subset aligner setting is the same
for every variant. The deviation from the paper is in the MAGUS baseline itself.
