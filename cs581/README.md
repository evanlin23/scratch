# CS581 project scratch: improving MAGUS (divide-and-conquer MSA)

Working folder for the CS581 (Algorithmic Computational Genomics) course project.
Starting point: MAGUS (Smirnov & Warnow, *Bioinformatics* 2021), which beat PASTA by
merging all subset alignments at once with the Graph Clustering Merger (GCM).

## Layout

| path | what |
|---|---|
| `papers/` | MAGUS paper (open access, CC BY-NC) + supplement + figures |
| `literature/sota_review.md` | Is MAGUS state of the art? What came after; project ideas; dataset links |
| `proposal/` | 2-page project proposal draft (`proposal.md`; `bash build.sh` -> PDF) |
| `code/setup.sh` | installs aligners, FastSP, MAGUS (pinned commit 39041fc), downloads the paper's data |
| `code/gcmx/` | experiment code: MAGUS merge variants via runtime patches (MAGUS itself untouched) |
| `code/fanout/` | parallel workers that rerun MAGUS with the paper's settings + all variants |
| `experiments/validation/` | rescoring of the paper's published alignments; comparison report |
| `experiments/runs/` | per-replicate results from the workers (`prep.json`, `pilot.jsonl`, cached inputs) |

## Apples-to-apples protocol

* **Datasets**: the paper's own files (Illinois Data Bank doi:10.13012/B2IDB-2643961_V1):
  ROSE 1000L1-3/M2-4/S1-3, RNASim 1K (20 reps each), 8 BAliBASE sets, 16S.M; 16S.3/16S.T/
  RNASim10K/16S.B.ALL for the large-data comparison.
* **Settings**: MAGUS(Fast) exactly as in the published logs: `--maxsubsetsize 0
  --maxnumsubsets K --decompstrategy pastastyle --decompskeletonsize 300 --graphbuildmethod mafft
  --graphbuildhmmextend false --graphclustermethod mcl --graphtracemethod minclusters
  --graphtraceoptimize false -r 10 -m 200 -f 4`, K = 25 (100 for 16S.3/16S.T, 50 for
  RNASim10K, 200 for 16S.B.ALL).
* **Metric**: (SPFN + SPFP) / 2 from FastSP against the same references (length-filtered
  "clean" references for biological data).
* **Paired design**: every merge variant reuses the subsets and backbones of the MAGUS run on
  the same replicate, so only the merge differs.

## Validation layers

1. `gcmx.validate_published`: FastSP on the authors' published MAGUS/PASTA alignments
   (streamed from Results.zip) vs the numbers in the paper's figures.
2. Our MAGUS reruns (paper flags) vs the published alignment of the same replicate.
3. Merge-only rerun on cached inputs (`default` variant) reproduces our full pipeline output.

`bash code/fanout/aggregate.sh` collects worker results and writes
`experiments/validation/REPORT.md` and `experiments/variants/SUMMARY.md`.

## Candidate improvements being tested (`code/gcmx/pilot.py`)

* `soft-m2`, `soft-m3`, `soft-m2-random`: **soft constraints**. Split each subset alignment
  into m groups and let GCM merge all groups at once (m = 1 is MAGUS).
* `progdp`, `progdp+opt`: exact progressive pairwise MWT merging on GCM's graph
  (+ leave-one-out refinement option).
* `default+opt`, `fm+opt`: published GCM search variants (Zaharias et al. 2023).
* `weight-frac`: normalized edge weights.
* `oracle:*`: true backbones / true subset alignments, to measure where the error comes from.
