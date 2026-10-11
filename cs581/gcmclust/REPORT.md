# Replacing MCL in MAGUS's Graph Clustering Merger (GCM)

Exploration pilot, CS581 (UIUC), Fall 2026. Instructor's boldfaced project: "Modify the MAGUS design by using a
different clustering method instead of the Markov Clustering algorithm within GCM." Merge-only: each replicate's
25 L-INS-i subset alignments and 10 L-INS-i backbones (one MAGUS draw) are fixed; only GCM's clustering step
changes, and MAGUS's own violation purge + `minclusters` trace + writer run unchanged on the new clusters.
Code: `code/` (`gc.py`, `build_graph.py`, `summarize.py`); raw rows: `results/*.results.jsonl`; protocol fixed
before results: `PREREG.md`; all tables: `results/tables_train.md`, `results/tables_held.md` (both generated before the exploratory runs), `results/tables_explore.md`.

## TL;DR (numbers first)

22 replicates (10 protein: 6 BAliBASE RV100 + 4 AliSim; 12 DNA/RNA: 9 ROSE 1000, 2 RNASim 1000, 16S.M), one
MAGUS draw each, train 10 / held-out 12 (gcmgen's split). Δ = variant − MAGUS in SP-error points
((SPFN+SPFP)/2 × 100), W/T/L at a 0.05 tie band, two-sided Wilcoxon. 31 clustering variants were run on training;
the 13 picked by the pre-registered rule were then run once on held-out.

| held-out (n = 12: 4 protein, 8 DNA/RNA) | all | proteins | DNA/RNA |
|---|---|---|---|
| **es4 + Leiden-CPM γ=0.02** (training pick) vs MAGUS | **−1.38** (7/1/4, p = 0.23) | **−4.24** (4/0/0) | +0.05 (3/1/4) |
| es4 + MCL I=4 (gcmgen's filter, MAGUS's MCL) vs MAGUS | −1.30 (7/4/1, p = 0.15) | −3.82 (3/0/1) | −0.04 (4/4/0) |
| es4 + label propagation vs MAGUS | −1.34 (7/2/3) | −4.04 (4/0/0) | +0.01 (3/2/3) |
| **es4 + Leiden-CPM vs es4 + MCL** (the clustering swap itself) | **−0.08** (6/2/4, p = 0.62) | **−0.42** (4/0/0) | +0.09 (2/2/4) |
| best raw-graph swap with no support filter: raw + LPA / MCL I=2 / Leiden-mod γ=100 | −0.53 / −0.10 / +0.01 | −2.39 / −0.39 / −1.45 | +0.40 / +0.05 / +0.75 |

Over all 10 protein replicates (training + held-out, γ was chosen on 6 of them), es4 + Leiden-CPM beats es4 + MCL
by **−0.52 (8/1/1, p = 0.010)**. Over the 12 DNA/RNA replicates it loses by +0.12 (3/2/7, p = 0.38). The gain is
fewer false positives (ΔSPFP −1.0) bought with a few more false negatives.

1. **Replacing MCL on MAGUS's raw graph does not help.** No community-detection method beats MCL I=4 on held-out
   (best: LPA −0.53, 5/1/6, p = 0.97, which helps proteins but costs DNA/RNA +0.40). Modularity methods (Leiden,
   Louvain) at low resolution build huge clusters and wreck DNA alignments (raw Leiden γ=10: +13.9 on training;
   held-out Louvain γ=30: +5.6 on DNA). Connected components (support ≥ 6, 8 or 10) are catastrophic
   (+3.5 to +22). Tuning MCL's inflation alone gives nothing: I = 2, 3, 6 → −0.11, +0.72, +0.17 on training;
   I=2 gives −0.10 on held-out.
2. **The low-support damage is in the graph, not specific to MCL.** Deleting edges supported by < 4 backbones
   helps on proteins whatever clusters the graph afterwards. Proteins, es4 − raw at the same parameter: MCL −1.8
   (train) / −3.8 (held), Leiden-CPM −2.2 / −4.9, Leiden-modularity −1.3 / −2.2, LPA −0.2 / −1.7, constrained
   agglomeration −7.0 (train). Once the graph is filtered, the clustering method matters ≤ 0.5 points
   (es4 + MCL, LPA and CPM all lie within 0.4 of each other on held-out).
3. **What does differ is cluster shape.** 48–73% of MCL's clusters break the one-column-per-subset rule (MAGUS's
   purge and trace then repair them). es4 drops that to 29–44%, and CPM to 2–20%. Protein alignments come out
   1.2–1.7× the reference length under every method. Clustering time is not a reason to switch (held-out mean,
   es4 graph: MCL 3 s, Leiden-CPM 5 s, LPA 1 s; raw graph: MCL 9 s, Leiden-CPM 19 s), but the es4 filter makes the
   trace 2–5× faster (held-out clustering + trace, graph already built: MAGUS 27 s; es4 + MCL 7 s; es4 + CPM 11 s).
4. **Constrained greedy agglomeration** (a Region-Growing relative) removes all violations but is slow in the
   trace (up to 400 s per protein set) and, on held-out, loses to es4 + MCL (+0.51, 1/1/10, p = 0.02).
5. **Verdict: unclear, leaning not promising as a stand-alone project.** The instructor's swap, done literally,
   gives no gain on the raw graph. It turns into a small protein-only gain (≈0.4–0.5 points over the support
   filter) only when combined with gcmgen's edge-support filter, and that gain is not significant on the 4
   held-out protein sets alone. It is a solid, cheap, *negative-plus-mechanism* project for 4 weeks; see §7.


## 1. Prior art

- **MAGUS** (Smirnov & Warnow, *Bioinformatics* 37:1666, 2021). GCM = alignment graph from backbones → MCL
  (inflation 4; the paper says the inflation sweep "follows a standard recommendation" and puts it in the
  supplement: 1000M1 only, values from van Dongen's 1.4/2/4/6) → trace by A* `minclusters`. MCL is motivated
  as extending the consistency principle to longer paths.
- **MWT-AM** (Zaharias, Smirnov & Warnow, *IEEE/ACM TCBB* 20:1700, 2023; AlCoB 2021). Formalises merging as Maximum
  Weight Trace (NP-hard). Tried other *clustering* steps — MLR-MCL and Region Growing (a constraint-aware greedy
  agglomeration close to our `agglo`) — and other *trace* steps (FM, MWT-greedy/search, RG-fast, plus an optimizer).
  Its headline win is in the trace (GCM(fm+opt) beat default GCM on every HomFam set); MCL stayed the default
  clustering afterwards (our notes do not record the clustering-variant numbers; to be checked in week 1). That work used MAGUS's raw graph and HomFam/ROSE, not a support-filtered graph, and
  did not test modularity/CPM community detection.
- **Recursive MAGUS** (Smirnov, *PLOS Comput Biol* 17:e1008950, 2021): recursion, guide trees, scale; GCM's
  clustering is unchanged (MCL); recursion buys scalability, not accuracy.
- **Community detection.** Louvain (Blondel et al. 2008), Leiden (Traag, Waltman & van Eck 2019; fixes
  Louvain's badly-connected communities; supports CPM, which has no resolution limit), label propagation
  (Raghavan et al. 2007). Our literature search (repo notes `cs581/literature/sota_review.md` §P5 and a web search)
  found **no alignment-graph merger that uses Leiden/Louvain/CPM**; residue-clustering MSA (MSARC, Modzelewski &
  Dojer 2014) uses hierarchical clustering of residues, not of subset-alignment columns.
- **In this repo.** `cs581/gcmgen` (branch claude/cs581-gcmgen): cross-subset GCM edges supported by < 4 of 10
  backbones are 2–25% correct; deleting them (`linsi#es4`) gives proteins −2.27 (7/1/2), DNA/RNA −0.02.

## 2. Methods

**Harness (`code/gc.py`).** For every replicate the GCM graph is built once, exactly as MAGUS builds it (a copy of
`run_edgesup.py`'s `buildMatrix`, `code/build_graph.py`), and saved with each edge's weight (residue-pair count
summed over backbones) and support (number of the 10 backbones that contribute to it). That run is also the
MAGUS control; it reproduces the cached MAGUS score (BBA0101 29.19, 1000M2 8.23, 1000L1 7.47, BBA0134 18.34; the
AliSim controls equal `bbtool_bench`'s merge-only control). Each variant then writes `graph/graph.txt` (raw or
filtered) and its own `graph/clusters.txt` into a fresh MAGUS working directory; MAGUS reads both instead of running
MCL, and runs its unchanged `purgeDuplicateClusters` → `purgeClusterViolations` → `minclusters` (A*) trace and
writer (`--graphtraceoptimize false`). Re-running MCL I=4 through this path gives MAGUS's score exactly, and
`es4:mcl:4` reproduces gcmgen's `linsi#es4` exactly on the shared cached replicates.

**Graphs.** `raw` = MAGUS's graph. `es4` = every non-self edge supported by < 4 of the 10 backbones deleted
(gcmgen's `linsi#es4`); the trace also sees the filtered graph.

**Clusterers** (all single-threaded; Leiden/Louvain/LPA seeded):

| name | method |
|---|---|
| `mcl:I` | MCL (MAGUS's bundled binary, `--abc -I I`, same input file as MAGUS) |
| `leidmod:γ` | Leiden (leidenalg 0.12), modularity (RBConfiguration) at resolution γ, raw edge weights |
| `leidcpm:γ` | Leiden, Constant Potts Model at resolution γ on degree-normalised weights w_ab/√(s_a s_b) (CPM has no resolution limit; the normalisation mimics MCL's column-stochastic scaling) |
| `louvain:γ` | Louvain (igraph multilevel), modularity at resolution γ |
| `lpa` | label propagation (igraph), weighted |
| `cc:t` | connected components of cross-subset edges with support ≥ t |
| `agglo:t` | greedy constrained agglomeration: cross-subset edges with support ≥ t in decreasing order of w_ab/√(s_a s_b); two clusters are joined only if they share no subset (union-find), so clusters never break the one-column-per-subset rule (no order constraint; the trace still orders) |

Infomap was not run (time). **Protocol:** `PREREG.md` (grid, gcmgen's train/held-out split, the selection rule:
per family × graph the parameter with the lowest mean training Δ; 30-min merge limit, no variant hit it).
**Data:** cached MAGUS draws (`cs581/experiments/runs`, worker branches) for BAliBASE, ROSE, RNASim and 16S.M;
AliSim SIMMOD/SIMHIGH R1 (train) and R2 (held-out) regenerated with cs581/protbench's seeds (`code/sim.sh`) plus one
fresh MAGUS draw each (paper flags, ~30 min each on the 4-core box). BAliBASE BBA0081/0117 dropped, as in gcmgen.

## 3. Training results (10 replicates, every grid point)

|---|---|---|---|---|---|---|---|
| `es4:agglo:1` | 10 | -0.64 | 6/1/3 | 0.131 | -0.16 / -1.12 | 0 + 55 (44) | 919 / 0.154 |
| `raw:agglo:1` | 10 | +3.81 | 2/2/6 | 0.084 | +9.17 / -1.55 | 1 + 278 (44) | 1172 / 0.044 |
| `raw:agglo:2` | 10 | +0.09 | 5/2/3 | 0.625 | +1.67 / -1.49 | 1 + 163 (44) | 1023 / 0.102 |
| `raw:agglo:4` | 10 | -1.04 | 8/0/2 | 0.037 | -1.66 / -0.42 | 0 + 48 (44) | 885 / 0.162 |
| `raw:cc:10` | 10 | +3.50 | 1/0/9 | 0.004 | +13.07 / -6.07 | 0 + 4 (44) | 859 / 0.256 |
| `raw:cc:6` | 10 | +22.01 | 0/0/10 | 0.002 | +49.82 / -5.80 | 0 + 8 (44) | 273 / 0.183 |
| `raw:cc:8` | 10 | +12.15 | 0/0/10 | 0.002 | +29.33 / -5.04 | 0 + 5 (44) | 486 / 0.214 |
| `es4:leidcpm:0.003` | 10 | -0.58 | 3/2/5 | 0.922 | -1.59 / +0.43 | 2 + 12 (44) | 754 / 0.148 |
| `es4:leidcpm:0.01` | 10 | -1.24 | 8/2/0 | 0.006 | -2.48 / +0.01 | 2 + 12 (44) | 882 / 0.159 |
| `es4:leidcpm:0.02` | 10 | -1.36 | 6/2/2 | 0.084 | -1.66 / -1.06 | 3 + 8 (44) | 990 / 0.186 |
| `raw:leidcpm:0.003` | 10 | +1.76 | 2/2/6 | 0.084 | +3.89 / -0.37 | 11 + 55 (44) | 908 / 0.043 |
| `raw:leidcpm:0.01` | 10 | +0.66 | 4/2/4 | 0.846 | +1.68 / -0.36 | 13 + 26 (44) | 1129 / 0.120 |
| `raw:leidcpm:0.02` | 10 | +0.48 | 4/1/5 | 0.492 | +3.26 / -2.29 | 16 + 16 (44) | 1189 / 0.183 |
| `es4:leidmod:10` | 10 | +5.25 | 5/0/5 | 0.375 | +11.46 / -0.96 | 2 + 3 (44) | 516 / 0.139 |
| `es4:leidmod:100` | 10 | -0.76 | 6/1/3 | 0.492 | -1.89 / +0.38 | 2 + 5 (44) | 684 / 0.139 |
| `es4:leidmod:30` | 10 | +0.90 | 6/0/4 | 0.922 | +2.18 / -0.37 | 2 + 3 (44) | 613 / 0.139 |
| `raw:leidmod:10` | 10 | +13.86 | 1/0/9 | 0.014 | +24.44 / +3.27 | 9 + 6 (44) | 374 / 0.025 |
| `raw:leidmod:100` | 10 | -0.15 | 3/2/5 | 0.922 | -0.61 / +0.30 | 10 + 30 (44) | 712 / 0.027 |
| `raw:leidmod:30` | 10 | +4.52 | 5/0/5 | 0.557 | +7.89 / +1.16 | 9 + 8 (44) | 542 / 0.025 |
| `es4:louvain:30` | 10 | +0.95 | 6/0/4 | 0.922 | +2.16 / -0.26 | 0 + 3 (44) | 612 / 0.139 |
| `raw:louvain:30` | 10 | +4.58 | 5/0/5 | 0.492 | +7.87 / +1.30 | 3 + 9 (44) | 544 / 0.025 |
| `es4:lpa` | 10 | -1.01 | 8/0/2 | 0.160 | -2.29 / +0.26 | 1 + 8 (44) | 778 / 0.139 |
| `raw:lpa` | 10 | -0.38 | 7/0/3 | 0.492 | +0.34 / -1.09 | 4 + 11 (44) | 650 / 0.025 |
| `es4:mcl:2` | 10 | -1.01 | 6/2/2 | 0.232 | -2.36 / +0.34 | 3 + 11 (44) | 817 / 0.133 |
| `es4:mcl:3` | 10 | -1.07 | 7/1/2 | 0.049 | -2.36 / +0.22 | 3 + 10 (44) | 857 / 0.133 |
| `es4:mcl:4` | 10 | -1.08 | 6/2/2 | 0.105 | -2.21 / +0.05 | 3 + 10 (44) | 879 / 0.134 |
| `es4:mcl:6` | 10 | -1.07 | 5/3/2 | 0.232 | -2.02 / -0.13 | 2 + 10 (44) | 904 / 0.135 |
| `raw:mcl:2` | 10 | -0.11 | 4/4/2 | 0.275 | -0.51 / +0.28 | 15 + 25 (44) | 771 / 0.019 |
| `raw:mcl:3` | 10 | +0.72 | 3/4/3 | 0.770 | +1.29 / +0.15 | 12 + 30 (44) | 812 / 0.020 |
| `raw:mcl:4` | 10 | +0.00 | 0/10/0 | – | +0.00 / +0.00 | 11 + 33 (44) | 840 / 0.021 |
| `raw:mcl:6` | 10 | +0.17 | 0/2/8 | 0.004 | +0.50 / -0.15 | 11 + 34 (44) | 874 / 0.023 |

#### train: proteins


Per replicate (Δ vs MAGUS) for the main variants:

| replicate | MAGUS err | `es4:mcl:4` | `es4:leidcpm:0.02` | `es4:lpa` | `raw:agglo:4` | `raw:lpa` | `raw:mcl:2` | `raw:leidmod:100` |
|---|---|---|---|---|---|---|---|---|
| BBA0101 | 29.19 | -0.39 | -1.42 | -0.13 | -1.33 | -2.56 | -0.80 | +0.37 |
| BBA0134 | 18.34 | -0.68 | -1.21 | -0.66 | -0.29 | -1.39 | +0.08 | -0.01 |
| BBA0067 | 26.28 | +0.11 | -1.72 | -0.09 | -1.02 | -1.90 | +0.04 | +0.44 |
| BBA0039 | 4.69 | -0.01 | -0.04 | -0.08 | +0.08 | -0.10 | -0.03 | -0.10 |
| SIMMOD_R1 | 13.78 | -4.53 | -4.32 | -4.27 | -3.80 | -0.97 | -0.11 | -1.62 |
| SIMHIGH_R1 | 23.97 | -5.45 | -5.81 | -5.09 | -3.65 | -2.20 | -0.24 | -1.56 |
| 1000M2 | 8.22 | +0.35 | +0.70 | +0.18 | +0.15 | +3.23 | +0.10 | +0.70 |
| 1000L1 | 7.47 | -0.10 | +0.01 | -0.09 | -0.06 | +0.18 | -0.04 | +0.09 |
| 1000L2 | 4.49 | -0.03 | +0.51 | +0.16 | -0.42 | +2.01 | -0.10 | +0.21 |
| 16S.M | 12.86 | -0.06 | -0.35 | -0.07 | -0.06 | -0.07 | +0.00 | -0.05 |



**Selection (PREREG rule).** es4 graph: MCL I=4 (I=3 at −1.07 ties I=4 at −1.08 within 0.01 → MAGUS's setting;
I=3 was also run on held-out and is identical), Leiden-CPM γ=0.02 (−1.36), Leiden-mod γ=100 (−0.76), Louvain γ=30
(+0.95), LPA (−1.01), agglo t=1 (−0.64). Raw graph: MCL I=2 (−0.11), Leiden-CPM γ=0.02 (+0.48), Leiden-mod γ=100
(−0.15), Louvain γ=30 (+4.58), LPA (−0.38), cc t=10 (+3.50), agglo t=4 (−1.04; note that t=4 already discards the
low-support edges, so it is a filtered-graph method). The overall training winner is es4 + Leiden-CPM γ=0.02.
γ=0.02 is the edge of the pre-registered CPM grid (see §6, exploratory).

## 4. Held-out results (12 replicates, evaluated once)

#### held: all replicates (n = 12)

| variant | n | mean Δ err vs `raw:mcl:4` | W/T/L | p | ΔSPFN / ΔSPFP | cluster s + trace s (base total s) | clusters / singleton frac |
|---|---|---|---|---|---|---|---|
| `es4:agglo:1` | 12 | -0.79 | 5/1/6 | 1.000 | -1.63 / +0.05 | 0 + 19 (27) | 1144 / 0.130 |
| `raw:agglo:4` | 12 | -1.15 | 5/3/4 | 0.791 | -2.48 / +0.19 | 1 + 18 (27) | 1114 / 0.133 |
| `es4:leidcpm:0.02` | 12 | -1.38 | 7/1/4 | 0.233 | -2.13 / -0.64 | 5 + 6 (27) | 1119 / 0.143 |
| `raw:leidcpm:0.02` | 12 | +0.89 | 3/1/8 | 0.110 | +3.25 / -1.48 | 19 + 27 (27) | 1413 / 0.110 |
| `es4:leidmod:100` | 12 | -0.77 | 3/1/8 | 0.339 | -1.82 / +0.28 | 4 + 6 (27) | 968 / 0.118 |
| `raw:leidmod:100` | 12 | +0.01 | 2/2/8 | 0.233 | -0.17 / +0.20 | 14 + 25 (27) | 973 / 0.024 |
| `es4:louvain:30` | 12 | +2.61 | 3/0/9 | 0.151 | +4.59 / +0.62 | 1 + 7 (27) | 869 / 0.118 |
| `es4:lpa` | 12 | -1.34 | 7/2/3 | 0.110 | -2.64 / -0.04 | 1 + 7 (27) | 1018 / 0.118 |
| `raw:lpa` | 12 | -0.53 | 5/1/6 | 0.970 | -0.95 / -0.11 | 4 + 16 (27) | 1009 / 0.024 |
| `es4:mcl:3` | 12 | -1.30 | 7/4/1 | 0.092 | -2.73 / +0.12 | 6 + 6 (27) | 1052 / 0.111 |
| `es4:mcl:4` | 12 | -1.30 | 7/4/1 | 0.151 | -2.69 / +0.08 | 3 + 4 (27) | 1059 / 0.111 |
| `raw:mcl:2` | 12 | -0.10 | 4/4/4 | 0.970 | -0.34 / +0.14 | 18 + 33 (27) | 1086 / 0.017 |

#### held: proteins

| variant | n | mean Δ err vs `raw:mcl:4` | W/T/L | p | ΔSPFN / ΔSPFP | cluster s + trace s (base total s) | clusters / singleton frac |
|---|---|---|---|---|---|---|---|
| `es4:agglo:1` | 4 | -2.98 | 4/0/0 | – | -5.86 / -0.09 | 0 + 40 (49) | 897 / 0.258 |
| `raw:agglo:4` | 4 | -3.48 | 2/1/1 | – | -7.17 / +0.21 | 0 + 42 (49) | 888 / 0.263 |
| `es4:leidcpm:0.02` | 4 | -4.24 | 4/0/0 | – | -8.07 / -0.41 | 3 + 12 (49) | 903 / 0.275 |
| `raw:leidcpm:0.02` | 4 | +0.66 | 2/0/2 | – | +2.85 / -1.52 | 8 + 70 (49) | 1506 / 0.159 |
| `es4:leidmod:100` | 4 | -3.62 | 2/0/2 | – | -7.63 / +0.40 | 3 + 12 (49) | 654 / 0.242 |
| `raw:leidmod:100` | 4 | -1.45 | 2/0/2 | – | -3.08 / +0.18 | 6 + 62 (49) | 667 / 0.055 |
| `es4:louvain:30` | 4 | -3.47 | 3/0/1 | – | -6.42 / -0.51 | 0 + 12 (49) | 588 / 0.242 |
| `es4:lpa` | 4 | -4.04 | 4/0/0 | – | -7.97 / -0.12 | 1 + 14 (49) | 742 / 0.242 |
| `raw:lpa` | 4 | -2.39 | 3/0/1 | – | -4.56 / -0.23 | 1 + 37 (49) | 731 / 0.054 |
| `es4:mcl:3` | 4 | -3.80 | 3/0/1 | – | -7.92 / +0.31 | 3 + 12 (49) | 821 / 0.231 |
| `es4:mcl:4` | 4 | -3.82 | 3/0/1 | – | -7.88 / +0.24 | 2 + 8 (49) | 836 / 0.231 |
| `raw:mcl:2` | 4 | -0.39 | 2/1/1 | – | -0.94 / +0.16 | 9 + 85 (49) | 914 / 0.042 |

#### held: DNA/RNA

| variant | n | mean Δ err vs `raw:mcl:4` | W/T/L | p | ΔSPFN / ΔSPFP | cluster s + trace s (base total s) | clusters / singleton frac |
|---|---|---|---|---|---|---|---|
| `es4:agglo:1` | 8 | +0.30 | 1/1/6 | 0.039 | +0.48 / +0.12 | 1 + 9 (16) | 1267 / 0.066 |
| `raw:agglo:4` | 8 | +0.02 | 3/2/3 | 0.844 | -0.14 / +0.18 | 1 + 7 (16) | 1226 / 0.069 |
| `es4:leidcpm:0.02` | 8 | +0.05 | 3/1/4 | 0.547 | +0.84 / -0.75 | 6 + 3 (16) | 1226 / 0.078 |
| `raw:leidcpm:0.02` | 8 | +1.00 | 1/1/6 | 0.078 | +3.45 / -1.45 | 24 + 6 (16) | 1366 / 0.086 |
| `es4:leidmod:100` | 8 | +0.65 | 1/1/6 | 0.039 | +1.09 / +0.21 | 5 + 4 (16) | 1124 / 0.056 |
| `raw:leidmod:100` | 8 | +0.75 | 0/2/6 | 0.023 | +1.28 / +0.21 | 18 + 6 (16) | 1126 / 0.009 |
| `es4:louvain:30` | 8 | +5.64 | 0/0/8 | 0.008 | +10.09 / +1.19 | 1 + 4 (16) | 1010 / 0.056 |
| `es4:lpa` | 8 | +0.01 | 3/2/3 | 1.000 | +0.02 / -0.01 | 1 + 3 (16) | 1156 / 0.056 |
| `raw:lpa` | 8 | +0.40 | 2/1/5 | 0.195 | +0.86 / -0.05 | 5 + 6 (16) | 1149 / 0.009 |
| `es4:mcl:3` | 8 | -0.05 | 4/4/0 | 0.148 | -0.13 / +0.03 | 7 + 4 (16) | 1168 / 0.051 |
| `es4:mcl:4` | 8 | -0.04 | 4/4/0 | 0.312 | -0.09 / -0.00 | 4 + 2 (16) | 1171 / 0.051 |
| `raw:mcl:2` | 8 | +0.05 | 2/3/3 | 0.547 | -0.03 / +0.13 | 23 + 7 (16) | 1172 / 0.004 |

Per replicate (Δ vs MAGUS):

| replicate | MAGUS err | `es4:mcl:4` | `es4:leidcpm:0.02` | `es4:lpa` | `raw:agglo:4` | `raw:lpa` | `raw:mcl:2` | `raw:leidmod:100` | `es4:agglo:1` |
|---|---|---|---|---|---|---|---|---|---|
| BBA0154 | 21.66 | +0.41 | -0.52 | -0.45 | +0.24 | -0.56 | +0.02 | +0.32 | -0.27 |
| BBA0190 | 23.40 | -0.30 | -0.70 | -0.29 | +0.02 | +0.07 | +0.25 | +0.72 | -0.24 |
| SIMMOD_R2 | 10.75 | -5.36 | -5.45 | -5.35 | -5.28 | -2.30 | -0.94 | -1.86 | -4.69 |
| SIMHIGH_R2 | 24.88 | -10.03 | -10.30 | -10.09 | -8.89 | -6.78 | -0.89 | -4.97 | -6.72 |
| 1000S1 | 9.78 | +0.01 | +0.32 | +0.08 | +0.02 | +0.10 | -0.02 | +0.10 | +0.14 |
| 1000L3 | 11.28 | -0.07 | +0.33 | +0.26 | +0.28 | +2.68 | +0.35 | +0.87 | +1.39 |
| 1000M3 | 4.52 | +0.01 | +0.06 | +0.07 | +0.24 | +0.17 | -0.07 | +0.01 | +0.39 |
| 1000M4 | 1.18 | -0.15 | -0.15 | -0.17 | -0.17 | -0.06 | +0.06 | -0.02 | -0.11 |
| 1000S2 | 4.74 | -0.11 | -0.23 | -0.13 | -0.20 | -0.16 | +0.03 | +0.14 | +0.01 |
| 1000S3 | 4.51 | +0.03 | +0.31 | -0.04 | +0.13 | +0.39 | +0.05 | +0.42 | +0.29 |
| RNASim | 9.92 | -0.08 | -0.23 | -0.07 | -0.12 | -0.00 | -0.05 | +2.42 | +0.07 |
| RNASim_R1 | 8.93 | +0.00 | -0.03 | +0.04 | -0.03 | +0.11 | +0.01 | +2.02 | +0.19 |

**The clustering swap itself (held-out, paired against es4 + MCL I=4):**

#### held: all replicates (n = 12)

| variant | n | mean Δ err vs `es4:mcl:4` | W/T/L | p | ΔSPFN / ΔSPFP | cluster s + trace s (base total s) | clusters / singleton frac |
|---|---|---|---|---|---|---|---|
| `es4:agglo:1` | 12 | +0.51 | 1/1/10 | 0.021 | +1.05 / -0.03 | 0 + 19 (7) | 1144 / 0.130 |
| `raw:agglo:4` | 12 | +0.16 | 2/4/6 | 0.266 | +0.20 / +0.11 | 1 + 18 (7) | 1114 / 0.133 |
| `es4:leidcpm:0.02` | 12 | -0.08 | 6/2/4 | 0.622 | +0.56 / -0.71 | 5 + 6 (7) | 1119 / 0.143 |
| `raw:leidcpm:0.02` | 12 | +2.19 | 3/1/8 | 0.092 | +5.94 / -1.56 | 19 + 27 (7) | 1413 / 0.110 |
| `es4:leidmod:100` | 12 | +0.53 | 0/3/9 | 0.001 | +0.87 / +0.20 | 4 + 6 (7) | 968 / 0.118 |
| `raw:leidmod:100` | 12 | +1.32 | 1/1/10 | 0.002 | +2.52 / +0.12 | 14 + 25 (7) | 973 / 0.024 |
| `es4:louvain:30` | 12 | +3.91 | 1/0/11 | 0.003 | +7.28 / +0.54 | 1 + 7 (7) | 869 / 0.118 |
| `es4:lpa` | 12 | -0.04 | 3/6/3 | 0.970 | +0.04 / -0.12 | 1 + 7 (7) | 1018 / 0.118 |
| `raw:lpa` | 12 | +0.77 | 2/0/10 | 0.021 | +1.74 / -0.19 | 4 + 16 (7) | 1009 / 0.024 |
| `es4:mcl:3` | 12 | +0.00 | 0/12/0 | 0.945 | -0.04 / +0.05 | 6 + 6 (7) | 1052 / 0.111 |
| `es4:mcl:4` | 12 | +0.00 | 0/12/0 | – | +0.00 / +0.00 | 3 + 4 (7) | 1059 / 0.111 |
| `raw:mcl:2` | 12 | +1.20 | 2/4/6 | 0.077 | +2.35 / +0.06 | 18 + 33 (7) | 1086 / 0.017 |
| `raw:mcl:4` | 12 | +1.30 | 1/4/7 | 0.151 | +2.69 / -0.08 | 9 + 18 (7) | 1139 / 0.018 |

#### held: proteins

| variant | n | mean Δ err vs `es4:mcl:4` | W/T/L | p | ΔSPFN / ΔSPFP | cluster s + trace s (base total s) | clusters / singleton frac |
|---|---|---|---|---|---|---|---|
| `es4:agglo:1` | 4 | +0.84 | 1/0/3 | – | +2.02 / -0.34 | 0 + 40 (10) | 897 / 0.258 |
| `raw:agglo:4` | 4 | +0.34 | 1/0/3 | – | +0.72 / -0.03 | 0 + 42 (10) | 888 / 0.263 |
| `es4:leidcpm:0.02` | 4 | -0.42 | 4/0/0 | – | -0.19 / -0.65 | 3 + 12 (10) | 903 / 0.275 |
| `raw:leidcpm:0.02` | 4 | +4.48 | 2/0/2 | – | +10.73 / -1.77 | 8 + 70 (10) | 1506 / 0.159 |
| `es4:leidmod:100` | 4 | +0.21 | 0/2/2 | – | +0.25 / +0.16 | 3 + 12 (10) | 654 / 0.242 |
| `raw:leidmod:100` | 4 | +2.37 | 1/0/3 | – | +4.81 / -0.06 | 6 + 62 (10) | 667 / 0.055 |
| `es4:louvain:30` | 4 | +0.35 | 1/0/3 | – | +1.46 / -0.75 | 0 + 12 (10) | 588 / 0.242 |
| `es4:lpa` | 4 | -0.22 | 2/2/0 | – | -0.08 / -0.36 | 1 + 14 (10) | 742 / 0.242 |
| `raw:lpa` | 4 | +1.43 | 1/0/3 | – | +3.33 / -0.47 | 1 + 37 (10) | 731 / 0.054 |
| `es4:mcl:3` | 4 | +0.02 | 0/4/0 | – | -0.03 / +0.07 | 3 + 12 (10) | 821 / 0.231 |
| `es4:mcl:4` | 4 | +0.00 | 0/4/0 | – | +0.00 / +0.00 | 2 + 8 (10) | 836 / 0.231 |
| `raw:mcl:2` | 4 | +3.43 | 1/0/3 | – | +6.94 / -0.08 | 9 + 85 (10) | 914 / 0.042 |
| `raw:mcl:4` | 4 | +3.82 | 1/0/3 | – | +7.88 / -0.24 | 3 + 46 (10) | 1033 / 0.044 |

#### held: DNA/RNA

| variant | n | mean Δ err vs `es4:mcl:4` | W/T/L | p | ΔSPFN / ΔSPFP | cluster s + trace s (base total s) | clusters / singleton frac |
|---|---|---|---|---|---|---|---|
| `es4:agglo:1` | 8 | +0.34 | 0/1/7 | 0.008 | +0.57 / +0.12 | 1 + 9 (6) | 1267 / 0.066 |
| `raw:agglo:4` | 8 | +0.06 | 1/4/3 | 0.742 | -0.06 / +0.18 | 1 + 7 (6) | 1226 / 0.069 |
| `es4:leidcpm:0.02` | 8 | +0.09 | 2/2/4 | 0.383 | +0.93 / -0.74 | 6 + 3 (6) | 1226 / 0.078 |
| `raw:leidcpm:0.02` | 8 | +1.04 | 1/1/6 | 0.039 | +3.54 / -1.45 | 24 + 6 (6) | 1366 / 0.086 |
| `es4:leidmod:100` | 8 | +0.70 | 0/1/7 | 0.008 | +1.18 / +0.22 | 5 + 4 (6) | 1124 / 0.056 |
| `raw:leidmod:100` | 8 | +0.79 | 0/1/7 | 0.008 | +1.37 / +0.21 | 18 + 6 (6) | 1126 / 0.009 |
| `es4:louvain:30` | 8 | +5.69 | 0/0/8 | 0.008 | +10.18 / +1.19 | 1 + 4 (6) | 1010 / 0.056 |
| `es4:lpa` | 8 | +0.05 | 1/4/3 | 0.461 | +0.11 / -0.00 | 1 + 3 (6) | 1156 / 0.056 |
| `raw:lpa` | 8 | +0.45 | 1/0/7 | 0.016 | +0.95 / -0.05 | 5 + 6 (6) | 1149 / 0.009 |
| `es4:mcl:3` | 8 | -0.01 | 0/8/0 | 0.312 | -0.05 / +0.03 | 7 + 4 (6) | 1168 / 0.051 |
| `es4:mcl:4` | 8 | +0.00 | 0/8/0 | – | +0.00 / +0.00 | 4 + 2 (6) | 1171 / 0.051 |
| `raw:mcl:2` | 8 | +0.09 | 1/4/3 | 0.195 | +0.05 / +0.13 | 23 + 7 (6) | 1172 / 0.004 |
| `raw:mcl:4` | 8 | +0.04 | 0/4/4 | 0.312 | +0.09 / +0.00 | 12 + 4 (6) | 1192 / 0.005 |



## 5. Mechanism: where does the low-support damage happen?

If the low-support edges hurt *because of MCL* (flow leaking along spurious edges, fragmenting or merging MCL
clusters), a different clustering method should be less sensitive to them, and es4 should help it less. It does not:
es4 − raw at the same clustering parameter (negative = filter helps), with each raw method's own Δ vs MAGUS:

Training:

| method:param | es4 − raw, proteins (n) | es4 − raw, DNA/RNA (n) | raw vs MAGUS, proteins | raw vs MAGUS, DNA/RNA |
|---|---|---|---|---|
| agglo:1 | -7.00 (6) | -0.62 (4) | +5.73 | +0.93 |
| agglo:2 | – | – | +0.10 | +0.08 |
| agglo:4 | – | – | -1.67 | -0.10 |
| cc:10 | – | – | +2.31 | +5.28 |
| cc:6 | – | – | +16.67 | +30.01 |
| cc:8 | – | – | +9.90 | +15.51 |
| leidcpm:0.003 | -3.75 (6) | -0.23 (4) | +2.37 | +0.85 |
| leidcpm:0.01 | -3.06 (6) | -0.14 (4) | +1.12 | -0.04 |
| leidcpm:0.02 | -2.15 (6) | -1.40 (4) | -0.27 | +1.62 |
| leidmod:10 | -2.26 (6) | -18.12 (4) | +1.65 | +32.16 |
| leidmod:100 | -1.34 (6) | +0.51 (4) | -0.41 | +0.24 |
| leidmod:30 | -1.11 (6) | -7.38 (4) | -0.75 | +12.43 |
| louvain:30 | -1.16 (6) | -7.34 (4) | -0.66 | +12.44 |
| lpa | -0.20 (6) | -1.29 (4) | -1.52 | +1.34 |
| mcl:2 | -1.49 (6) | +0.00 (4) | -0.18 | -0.01 |
| mcl:3 | -2.98 (6) | -0.00 (4) | +1.24 | -0.05 |
| mcl:4 | -1.82 (6) | +0.04 (4) | +0.00 | +0.00 |
| mcl:6 | -2.06 (6) | -0.04 (4) | +0.13 | +0.23 |

Held-out:

| method:param | es4 − raw, proteins (n) | es4 − raw, DNA/RNA (n) | raw vs MAGUS, proteins | raw vs MAGUS, DNA/RNA |
|---|---|---|---|---|
| agglo:1 | – | – | – | – |
| agglo:4 | – | – | -3.48 | +0.02 |
| cc:10 | – | – | -0.59 | +7.48 |
| leidcpm:0.02 | -4.91 (4) | -0.95 (8) | +0.66 | +1.00 |
| leidmod:100 | -2.17 (4) | -0.09 (8) | -1.45 | +0.75 |
| louvain:30 | -1.04 (4) | -5.96 (8) | -2.43 | +11.60 |
| lpa | -1.65 (4) | -0.40 (8) | -2.39 | +0.40 |
| mcl:2 | – | – | -0.39 | +0.05 |
| mcl:3 | – | – | – | – |
| mcl:4 | -3.82 (4) | -0.04 (8) | +0.00 | +0.00 |

- On proteins the filter helps **every** method, usually by as much as or more than it helps MCL (CPM −2.2/−4.9,
  agglomeration −7.0, MCL −1.8/−3.8). The one partial exception is LPA (training −0.2, held-out −1.7): on the raw
  graph LPA already gets part of the filter's gain on proteins (−1.5/−2.4 vs MAGUS) because a node adopts its
  heaviest neighbourhood label and a few light spurious edges rarely win, but LPA's giant clusters cost DNA/RNA
  (+1.3/+0.4).
- So the damage is done by the evidence (edges that are 2–25% correct, gcmgen) entering *any* clustering, and
  then the trace. The clustering step matters mainly as **how big a cluster may grow along weak edges**:
  modularity/Louvain/LPA/connected components grow large, constraint-violating clusters (max size 200–1,500 on
  proteins, 50–950 on DNA) that the trace must break apart, which inflates alignment length and false negatives;
  CPM on normalised weights and constrained agglomeration cap cluster size near the number of subsets (≤ 25–28).
- After the filter, the method barely matters on DNA/RNA (all within ±0.1 except agglo, Leiden-mod and Louvain,
  which lose) and matters ≈ 0.4–0.5 points on proteins, where CPM's smaller, cleaner clusters trade some recall for
  precision (held-out proteins ΔSPFN −0.19, ΔSPFP −0.65 vs es4 + MCL).

Cluster statistics (training; held-out in `results/tables_held.md`):

#### cluster statistics, proteins (train)

| variant | clusters (size ≥ 2) | mean size | max size | singleton nodes | clusters violating 1-col/subset | aln length / ref |
|---|---|---|---|---|---|---|
| `raw:mcl:4` | 704 | 20.3 | 81 | 0.033 | 0.60 | 1.44 |
| `raw:mcl:2` | 604 | 23.5 | 169 | 0.030 | 0.64 | 1.52 |
| `es4:mcl:4` | 750 | 16.7 | 59 | 0.189 | 0.44 | 1.43 |
| `es4:mcl:3` | 723 | 17.1 | 74 | 0.188 | 0.46 | 1.44 |
| `raw:leidcpm:0.02` | 1130 | 10.8 | 26 | 0.237 | 0.16 | 1.69 |
| `es4:leidcpm:0.02` | 889 | 13.1 | 27 | 0.255 | 0.21 | 1.54 |
| `raw:leidmod:100` | 527 | 26.2 | 228 | 0.040 | 0.66 | 1.49 |
| `es4:leidmod:100` | 482 | 23.9 | 204 | 0.195 | 0.51 | 1.72 |
| `raw:louvain:30` | 387 | 36.0 | 403 | 0.038 | 0.65 | 2.01 |
| `es4:louvain:30` | 425 | 27.5 | 345 | 0.195 | 0.49 | 2.03 |
| `raw:lpa` | 429 | 34.1 | 1475 | 0.038 | 0.58 | 2.04 |
| `es4:lpa` | 610 | 20.0 | 378 | 0.195 | 0.49 | 1.58 |
| `raw:cc:10` | 546 | 18.3 | 562 | 0.351 | 0.07 | 2.89 |
| `raw:agglo:4` | 748 | 15.2 | 25 | 0.221 | 0.00 | 1.35 |
| `es4:agglo:1` | 769 | 14.9 | 25 | 0.211 | 0.00 | 1.40 |


#### cluster statistics, DNA/RNA (train)

| variant | clusters (size ≥ 2) | mean size | max size | singleton nodes | clusters violating 1-col/subset | aln length / ref |
|---|---|---|---|---|---|---|
| `raw:mcl:4` | 1044 | 25.8 | 46 | 0.003 | 0.73 | 0.90 |
| `raw:mcl:2` | 1020 | 26.4 | 54 | 0.003 | 0.74 | 0.92 |
| `es4:mcl:4` | 1072 | 23.8 | 36 | 0.052 | 0.38 | 0.89 |
| `es4:mcl:3` | 1058 | 24.2 | 38 | 0.052 | 0.39 | 0.89 |
| `raw:leidcpm:0.02` | 1278 | 19.1 | 26 | 0.102 | 0.02 | 1.14 |
| `es4:leidcpm:0.02` | 1141 | 21.7 | 27 | 0.084 | 0.07 | 0.96 |
| `raw:leidmod:100` | 989 | 27.2 | 80 | 0.007 | 0.74 | 0.99 |
| `es4:leidmod:100` | 987 | 25.9 | 94 | 0.056 | 0.41 | 1.04 |
| `raw:louvain:30` | 778 | 35.0 | 230 | 0.007 | 0.75 | 2.36 |
| `es4:louvain:30` | 892 | 28.7 | 183 | 0.056 | 0.42 | 1.64 |
| `raw:lpa` | 980 | 27.5 | 956 | 0.007 | 0.74 | 1.16 |
| `es4:lpa` | 1030 | 24.8 | 66 | 0.056 | 0.40 | 0.92 |
| `raw:cc:10` | 1328 | 19.3 | 398 | 0.113 | 0.02 | 1.61 |
| `raw:agglo:4` | 1092 | 23.0 | 25 | 0.073 | 0.00 | 0.86 |
| `es4:agglo:1` | 1144 | 22.1 | 25 | 0.069 | 0.00 | 0.87 |


MCL's clusters are mostly *not* valid trace columns: 48–73% contain two columns from one subset (MAGUS's purge
removes the lower-scoring ones, then the trace splits crossing clusters). The es4 filter lowers that to 29–44%; CPM makes
almost all clusters valid. Despite that, CPM's held-out alignment lengths are no closer to the reference than MCL's
(proteins 1.35 vs 1.32), i.e. the over-long protein alignments (1.2–1.4× the reference) come from the graph/trace,
not from MCL.

## 6. Runtime, trees, exploratory

**Runtime** (held-out means, single thread, graph already built; MAGUS's own graph building is the same for every
variant and is not included): MAGUS (raw, MCL I=4) 9 s clustering + 18 s trace; es4 + MCL 3 + 4 s; es4 + Leiden-CPM
5 + 6 s; es4 + LPA 1 + 7 s; raw + Leiden-CPM 19 + 27 s; agglo 0–1 + 18–44 s (A* trace struggles on its many small
clusters; up to 400 s per protein set on training). The filter, not the clusterer, is what saves time.

**Trees** (FastTree `-lg -gamma`, normalised RF to the true AliSim tree; 4 simulated protein sets, R1 = training,
R2 = held-out; `code/trees.py`, `results/trees.jsonl`):

| replicate | `true` | `raw:mcl:4` | `es4:mcl:4` | `es4:leidcpm:0.02` | `es4:lpa` | `raw:agglo:4` |
|---|---|---|---|---|---|---|
| SIMMOD_R1 | 0.0652 | 0.0642 | 0.0632 | 0.0612 | 0.0642 | 0.0632 |
| SIMHIGH_R1 | 0.0612 | 0.1254 | 0.1224 | 0.1254 | 0.1214 | 0.1204 |
| SIMMOD_R2 | 0.0652 | 0.0652 | 0.0672 | 0.0722 | 0.0682 | 0.0712 |
| SIMHIGH_R2 | 0.0522 | 0.1043 | 0.0752 | 0.0772 | 0.0863 | 0.0782 |
| mean | 0.0609 | 0.0898 | 0.0820 | 0.0840 | 0.0850 | 0.0833 |

The edge-support filter lowers tree error (mean 0.090 → 0.082, almost all from SIMHIGH_R2: 0.104 → 0.075). Swapping
MCL for Leiden-CPM on top of it does not lower it further (0.084; on SIMMOD_R2 it is the worst of the variants).
The −0.4 to −0.5 SP-point protein gain of CPM does not show up in trees at n = 4.

**Exploratory (after the held-out evaluation; not used for any claim above): Leiden-CPM resolution beyond the
grid edge**, es4 graph, all 22 replicates, Δ vs es4 + MCL I=4 (cluster counts and singleton fractions in
`results/tables_explore.md`):

| γ | proteins (n = 10) | DNA/RNA (n = 12) | singleton nodes |
|---|---|---|---|
| 0.01 | −0.05 (5/2/3) | −0.07 (6/4/2) | 0.14 |
| **0.02** (selected) | **−0.52 (8/1/1, p = 0.01)** | +0.12 (3/2/7) | 0.16 |
| 0.03 | −0.69 (5/1/4) | +1.62 (1/0/11) | 0.19 |
| 0.05 | +36.3 (0/0/10) | +41.9 (0/0/12) | 0.87 |

γ = 0.02 sits just below a cliff: by 0.05 CPM leaves most columns as singletons and the alignment falls apart.
Proteins prefer γ ≈ 0.02–0.03 (precision), DNA/RNA γ ≤ 0.01 (recall). A per-dataset γ would need a
reference-free signal, and the cliff makes a fixed γ fragile. That is the main technical risk of the CPM variant.

## 7. Verdict for a 4-week CS581 project: **unclear** (leaning not promising on its own)

- *For:* boldfaced (straightforward) project; the harness exists and runs a full train/held-out sweep of 22
  replicates in ~3 h on 4 cores; the question has a clean, already-answered-once structure; a negative result plus
  the mechanism analysis above is a respectable course report.
- *Against:* on MAGUS's own graph, nothing we tried beats MCL I=4 (13 variants on held-out; best −0.53 with a
  DNA/RNA cost). The only gain is es4 + Leiden-CPM over es4 + MCL on proteins: −0.42 held-out (4/0/0, n too small
  for a test), −0.52 over all 10 protein replicates (p = 0.01, but 6 of them chose γ), and +0.09/+0.12 on DNA/RNA.
  The headline improvement in every table is gcmgen's edge-support filter, not the clustering.
- If chosen, frame it as **"Is MCL the right clusterer for GCM? A controlled comparison"**: the deliverable is the
  comparison and the mechanism (cluster validity, size caps, interaction with evidence filtering), with es4 +
  Leiden-CPM as the candidate improvement for proteins.

**Weeks 1–4 (if chosen).**
1. Week 1: read MWT-AM's clustering experiments (MLR-MCL, RG) and reproduce one; re-run this harness on 4–6 more
   protein sets (BBA0081/0117, 10AA, HomFam 2k with seed-only scoring, AliSim R3/R4) so the held-out protein n ≥ 8;
   fix γ for CPM on training (extend the grid past 0.02, §6).
2. Week 2: second MAGUS draw per dataset (draw-to-draw variance); a constrained Leiden (refuse moves that put two
   columns of one subset together) and Infomap; CPM γ adapted per data set from a reference-free statistic
   (e.g. fraction of low-support edges).
3. Week 3: end-to-end MAGUS runs with the best variant on a subset (ROSE 1000M2/L1, RNASim, 16S.M, BAliBASE), and
   FastTree/IQ-TREE nRF on AliSim.
4. Week 4: write-up; figures (Δ vs MAGUS by family, cluster-validity vs error, es4 interaction).

**Risks.** (i) The gain is mostly the support filter, which is someone else's (gcmgen) result; the clustering swap
alone may come out null — plan the report so that a null is a finding. (ii) Small protein n and one MAGUS draw per
dataset; the AliSim gains are large (−5 to −10) but come from the filter. (iii) CPM's resolution needs a scale; the
w/√(s_a s_b) normalisation is a choice, and γ at the grid edge means the optimum may be elsewhere. (iv) Leiden is
fast, but A* trace time explodes for some cluster shapes (agglo); any new clusterer must be watched for that.
(v) DNA/RNA: anything that trades recall for precision hurts there (ROSE is recall-limited, gcmgen).

## 8. Caveats

One MAGUS draw per replicate; 4 held-out protein sets; Leiden is randomised (seed 1 only); the BBA0101 debugging
probe was seen before the grid was fixed (it is a training set); tree results are on 4 AliSim sets only. All
clustering variants keep MAGUS's trace (`minclusters`, no optimizer); MWT-AM's fm+opt trace might change the ranking.

