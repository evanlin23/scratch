# Local single-replicate findings (1000M2_R0), 2026-10-09

Exploratory runs on one replicate (1000M2 R0, MAGUS paper settings, same 25
subset alignments unless stated). Error = (SPFN+SPFP)/2 in %. MAGUS's own
merge on these inputs: **8.23**. Single replicate: these numbers pick
directions; significance comes from the fan-out workers.

## Runtime-oriented ideas

| idea | result | verdict |
|---|---|---|
| MAGUS `--graphbuildmethod initial` (guide-tree alignment as the only evidence, no MAFFT backbones), end to end | 15.20 | MAGUS bug: its evidence file holds only the 700 HMM-added queries, not the 300 skeleton sequences |
| same, evidence fixed (skeleton alignment + HMM-extended rest), same subsets | 15.13 | a single HMM-added alignment is weak evidence |
| HMM-extended subset alignments as evidence ("subsethmm") | 38.10 (SPFN 69.9) | fails |
| smaller backbones: 10 x 100 seqs (MAFFT CPU 23.5 vs 68.5 min for 10 x 200) | hard 10.70, Slow 9.92, Slow+soft m3 10.00 | accuracy loss ~1.5-2 pts; soft does not recover it |
| smaller backbones: 10 x 50 seqs (CPU 8.7 min) | hard 11.93, Slow 11.33, Slow+soft m3 12.30 | worse |
| fresh 10 x 200 backbones (different random draw, same subsets) | hard 8.82, Slow 8.94, Slow+soft m3 8.43 | **backbone draw alone moves error by ~0.6 pts** |
| MCL with `-te 3` (soft-merge graph, 160M edges) | identical clusters; 317 s vs 699 s | safe 2.2x speedup of the merge bottleneck (`--gcmx-mclthreads`) |

MAGUS stage costs (published logs, 1000M2 R0): initial tree 308 s, subset
alignments + backbones until ~921 s, merge ~47 s of 968 s. The soft-constraint
merge spends most of its time in MCL (210 of 292 s single-threaded).

## Accuracy-oriented ideas

| idea | result |
|---|---|
| self-derived evidence (MAGUS output restricted to 8 seqs/subset, HMM-extended) + soft m3, one round | **7.34** (234 s) |
| ... rounds 2, 3, 4 | 7.39, 7.38, 7.40 (no further gain) |
| adaptive splitting (tiers / tiers-hard / proportional) with PP-weighted Slow evidence | 7.74-7.99 vs 7.79 uniform m3: no gain |
| 24 PP-weighted extended backbones + soft m3 / m4 | 7.19 / 7.12 (more backbones = more runtime) |
| learned edge weights, trained and tested on this replicate (in-sample, optimistic) | 7.79 (same merge time); held-out test across 37 replicates running |
