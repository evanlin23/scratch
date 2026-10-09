## 1. Paper figures vs. our rescoring of the published alignments

Average error (SPFN+SPFP)/2 in %. *Paper* = read off the paper's figures (±0.5); *rescored* = FastSP on the alignments the authors published (Illinois Data Bank), mean over replicates.

| dataset | reps | MAGUS(Fast) paper | MAGUS(Fast) rescored | PASTA paper | PASTA rescored | MAGUS(Slow) rescored | PASTA(3)+GCM rescored |
|---|---|---|---|---|---|---|---|
| 1000L3 | 20 | 11.0 | 11.1 | 17.5 | 17.5 | 10.4 | 12.4 |
| 1000S1 | 20 | 11.1 | 11.1 | 15.7 | 15.7 | 10.7 | 12.3 |
| 1000M2 | 20 | 9.8 | 9.8 | 13.7 | 13.7 | 9.1 | 10.8 |
| RNASim | 20 | 9.6 | 9.6 | 9.8 | 9.8 | 9.3 | 9.2 |
| 1000L1 | 20 | 6.2 | 6.2 | 8.0 | 8.0 | 5.8 | 7.2 |
| 1000S2 | 20 | 5.4 | 5.4 | 7.6 | 7.6 | 5.2 | 6.0 |
| 1000S3 | 20 | 4.5 | 4.5 | 6.2 | 6.2 | 4.4 | 4.9 |
| 1000M3 | 20 | 3.5 | 3.5 | 4.8 | 4.8 | 3.4 | 4.0 |
| 1000L2 | 20 | 2.9 | 2.9 | 3.7 | 3.7 | 2.9 | 3.2 |
| 1000M4 | 20 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.2 |
| 1000M1 | 19 | 12.2 | 12.1 | – | 16.5 | 11.5 | 13.9 |
| 16S.3 | 1 | 10.3 | 10.3 | 10.3 | 10.3 | 10.2 | 10.3 |
| 16S.T | 1 | 10.0 | 9.9 | 12.9 | 12.9 | 9.8 | 10.0 |
| RNASim10K | 10 | 8.3 | 8.3 | 10.8 | 10.7 | – | 8.4 |
| 16S.B.ALL | 1 | 4.2 | 4.1 | 5.4 | 5.4 | – | 4.1 |
| BBA0081 | 1 | 57.0 | 57.0 | 74.1 | 74.1 | 55.3 | 67.1 |
| BBA0067 | 1 | 25.6 | 25.6 | 25.7 | 25.7 | 25.8 | 24.8 |
| BBA0154 | 1 | 21.1 | 21.1 | 22.4 | 22.5 | 20.9 | 21.3 |
| BBA0101 | 1 | 28.5 | 28.5 | 29.2 | 29.2 | 29.3 | 29.2 |
| BBA0190 | 1 | 23.1 | 23.1 | 23.9 | 23.9 | 23.3 | 23.6 |
| BBA0134 | 1 | 18.7 | 18.7 | 20.3 | 20.4 | 18.6 | 18.4 |
| BBA0117 | 1 | 12.1 | 12.1 | 12.1 | 12.1 | 12.8 | 12.0 |
| BBA0039 | 1 | 4.5 | 4.5 | 4.2 | 4.3 | 4.5 | 4.4 |
| 16S.M | 1 | 13.0 | 13.1 | 13.0 | 13.0 | 12.9 | 12.8 |

## 2. Our MAGUS reruns vs. the published MAGUS alignment of the same replicate

Our runs: current MAGUS (commit 39041fc) with the paper's MAGUS(Fast) flags (`--maxsubsetsize 0 --maxnumsubsets 25 -r 10 -m 200 --graphbuildhmmextend false`), 4 cores. Backbone sampling is random, so equal-in-expectation, not identical, results are expected.

| replicate | published MAGUS(Fast) | ours | diff (pts) | published PASTA(3) | our minutes (4 cores) | published minutes |
|---|---|---|---|---|---|---|

## 3. Harness check: merge-only rerun on cached inputs reproduces the full pipeline

| replicate | full pipeline | merge-only rerun (`default`) | identical error? |
|---|---|---|---|
