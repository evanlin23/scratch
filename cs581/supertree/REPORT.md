# Pilot: divide-and-conquer supertrees with a Disjoint Tree Merger (DC-GTM)

CS581 Fall 2026 project pilot, 2026-10-09/10. Branch `claude/cs581-supertree`.
Code: `code/`. Raw per-run results: `results/*.jsonl`. Tables: `results/*_summary.md`. Scan: `SCAN.md`.

## 1. Question and where it comes from

- **First-day deck**, slide "Open problems (and possible course projects)": *"Supertree estimation: Need to scale to 100,000+ species with high accuracy; Approximation algorithms."*
- **Divide-and-conquer deck**: DACTAL "depends on having a scalable and accurate supertree method"; "supertree methods are generally computationally intensive". DTMs (NJMerge, TreeMerge, GTM) are presented as the scalable alternative, and GTM+ASTRAL was "faster and more accurate than ASTRAL" for species trees from gene trees.

**Pilot question.** Supertree input is a set of source trees on overlapping taxon subsets. Can a DTM pipeline give a supertree method that matches or beats the accurate-but-unscalable methods (ASTRAL-III, FastRFS) while scaling to 10^4–10^5 taxa?

The pipeline:
1. compute a cheap guide tree;
2. split the taxa into disjoint subsets of ≤ m taxa by centroid-edge decomposition of the guide;
3. restrict every source tree to each subset and run ASTRAL-III on each subset;
4. merge the subset trees with GTM using the guide.

## 2. Scan (summary; full tables with DOIs in `SCAN.md`)

| method | type | largest published test | runtime at scale | status here |
|---|---|---|---|---|
| MRP (PAUP*/TNT) | parsimony on the MRP matrix | ~5.5k (STK) | did not finish in 14 days on SMIDGenOG-5500 (cited in the SCS paper) | not run (no PAUP*) |
| MRL | ML on the MRP matrix | 2,228 (CPL) | 575 s on 558 taxa (RAxML) | **MRL-FT** (FastTree on the MRP matrix): 15–43 s at 1–2k taxa, **1,033 s at 10k** |
| SuperFine (+MRP/MRL) | SCM + polytomy refinement | 2,228 | the SCM collapses at large n | not run (python2, PAUP*) |
| BCD (Fleischauer & Böcker 2017) | min-cut + clade deletion | 10k (SCS-DCM) | ~2 h at 10k | not run |
| FastRFS (2017) | exact RFS in a constrained space | 2,228 | 3,282 s on CPL | not built (needs an old Bazel); published numbers used |
| Exact-RFS-2 / GreedyRFS (*AMB* 2021) | exact RFS of 2 trees, O(n²\|X\|); greedy pairwise | 500 taxa, 9 replicates (criterion scores only) | not reported | not run |
| ASTRAL-II/III | max quartet support | 1,000 (FastRFS paper) | **ASTRAL-III here: 21 s and 1.5 GB at 1k, 185 s and 3.7 GB at 2k; killed (exit 137) at 3.9k, 5k and 10k with an 8 GB heap under shared memory; heap-13g runs alone: see §6** | **main baseline** |
| ASTER astral4 | quartet placement + subsampling | — | 45 s at 1k | **32.1% RF vs ASTRAL-III 17.2% on SMIDGen-1000 d20** (bad on scaffold supertree input; `-R`: 28.9%) |
| TREE-QMC (2023) | weighted quartet Max-Cut | "promising as supertree" | 340 s at 1k (O(n³k)) | 22.0% RF on one SMIDGen-1000 d20 replicate (ASTRAL-III 17.4%) |
| wQFM-TREE (2025), Asteroid (2023) | quartet FM / balanced minimum evolution with missing data | species-tree data | — | not run |
| Spectral Cluster Supertree (2024) | spectral clustering of rooted trees | 10k | 20 s (paper); **65 s at 10k here** | **run: the best method on the rooted DCM data** |
| SDSR (arXiv 2026) | spectral D&C for species trees | — | up to 10× faster than ASTRAL (abstract) | not run |
| GTM / NJMerge / TreeMerge | DTMs (disjoint subsets) | 50k (GTM, ML trees) | GTM: < 1 s for 10k taxa here | **GTM used as the merger** |

**Datasets with true trees** (details in `SCAN.md`):
- SMIDGen 100/500/1,000 (doi:10.13012/B2IDB-2952208_V1);
- SMIDGenOG 100–1,000 and SMIDGenOG-5500, plus the SCS birth-death DCM-IQ data at 500–10,000 taxa (Zenodo doi:10.5281/zenodo.11118022);
- TreeMerge species-tree data (doi:10.13012/B2IDB-9570561_V1): a DTM benchmark, not a supertree benchmark;
- Exact-RFS-2: no separate deposit; it reuses SMIDGen-500.

**No public supertree benchmark with true trees beyond 10k taxa exists**, so "100,000+" would need a new simulation.

## 3. Chosen idea and why

DC-GTM (decompose by guide tree → ASTRAL-III on subsets → GTM). Reasons:
- (i) it is the D&C-lecture recipe applied to the supertree problem, and I found no published GTM-for-supertrees study;
- (ii) GTM runs in linear time and ASTRAL-III is fast on ≤ 500-taxon subsets, so it scales to 10^5 if a cheap guide exists;
- (iii) it can be fully piloted in hours with off-the-shelf tools.

Guides tried:
- the ASTRAL-III tree itself (a sanity check);
- MRL-FT (cheap, unrooted);
- SCS (cheap, needs rooted source trees);
- the true tree (an oracle, to bound the headroom).

## 4. Data and baseline reproduction

- **SMIDGen 1,000 taxa** (doi:10.13012/B2IDB-2952208_V1): 4 scaffold densities × 10 replicates, ~25 ML source trees each (one scaffold plus clade-based trees). The model trees have ~930 leaves (all taxa that appear in a source tree).
- **SCS birth-death DCM-IQ** (Zenodo 11118022): 500 / 1k / 2k / 5k / 10k taxa, replicates 0–4 (0–2 at 10k). Source trees are IQ-TREE 2 trees on Rec-I-DCM3 subsets of ≤ 100 taxa (mean 44 at 10k), rooted. There are 250 source trees at 10k, and **each taxon is in ~1.1 source trees (very little overlap)**.
- **SMIDGenOG-5500** replicates 0–4 (3,910 taxa in replicate 0; 503 RAxML source trees of 74–124 taxa; each taxon in ~13 trees; mean source-tree RF 7%).

**Reproduction (ASTRAL as a supertree method on SMIDGen-1000; RF error %; published = ASTRAL-II, FastRFS paper Table 2):**

| scaffold density | 20% | 50% | 75% | 100% |
|---|---|---|---|---|
| published ASTRAL | 16.9 | 15.7 | 13.6 | 11.6 |
| **ours, ASTRAL-III 5.7.8 (10 reps)** | **17.19** | **15.76** | **13.85** | **11.62** |
| published FastRFS-enhanced (for reference) | 16.7 | 15.1 | 13.4 | 11.8 |

All four cells are within 0.3 pp of the published numbers, so the data, the metric (unrooted RF on the common leaf set; FN = FP for binary trees) and the baseline are validated. FastRFS-enhanced is at most 0.6 pp better than ASTRAL in the published table, which makes ASTRAL-III a fair proxy for the best accurate RFS/quartet baseline.

## 5. Results

Metric: RF error % against the true tree. FN and FP are in `results/*_summary.md`; they differ by < 0.3 pp except for SCS, which outputs polytomies. Paired by replicate:
- difference = new − baseline, in percentage points;
- **W/T/L** = new wins / ties / losses, with a tie band of |Δ| ≤ 0.25 pp;
- two-sided Wilcoxon signed-rank test (with n = 5 the smallest attainable p is 0.0625).

### 5.1 SMIDGen-1000 (40 replicates)

| method | d20 | d50 | d75 | d100 | wall s (incl. guide) |
|---|---|---|---|---|---|
| ASTRAL-III | **17.19** | 15.76 | 13.85 | 11.62 | 21–24 |
| MRL-FT | 19.64 | 20.88 | 17.42 | 11.93 | 17–27 |
| ASTER astral4 (8 reps) | 32.07 | – | – | – | 45 |
| DC-GTM, ASTRAL-III guide, m = 100 | 17.71 | **15.55** | **13.74** | **11.35** | ~34 (= 22 + 12) |
| DC-GTM, MRL-FT guide, m = 100 | 20.10 | 19.92 | 16.76 | 11.69 | ~33 |
| DC-GTM, MRL-FT guide, m = 200 | 19.65 | 19.27 | 16.52 | 11.67 | ~34 |
| *DC-GTM, true tree as guide (oracle)* | *13.42* | *12.48* | – | – | 12 |

Paired, all 40 replicates:

| new vs baseline | mean Δ (pp) | W/T/L | p |
|---|---|---|---|
| DC-GTM(ASTRAL guide) vs ASTRAL-III | −0.02 | 21/5/14 | 0.65 |
| DC-GTM(MRL-FT guide, m = 100) vs ASTRAL-III | +2.51 | 5/3/32 | 4.5e-7 |
| DC-GTM(MRL-FT guide, m = 200) vs ASTRAL-III | +2.17 | 2/7/31 | 8.4e-7 |
| DC-GTM(MRL-FT guide, m = 200) vs MRL-FT | −0.69 | 24/7/9 | 0.0016 |
| MRL-FT vs ASTRAL-III | +2.86 | 4/1/35 | 2.8e-7 |
| astral4 vs ASTRAL-III (d20, n = 8) | +13.99 | 0/0/8 | 0.0078 |

### 5.2 Birth-death DCM-IQ, 500–10,000 taxa (rooted source trees with little overlap)

| method | 500 | 1k | 2k | 5k | 10k |
|---|---|---|---|---|---|
| SCS | **1.83** | **2.80** | **3.20** | **3.39** | **3.84** (3 reps) |
| ASTRAL-III | 5.71 | 8.08 | 7.23 | see §6 | see §6 |
| MRL-FT | 6.84 | 8.30 | 7.31 | – | 8.03 (1 rep) |
| DC-GTM, SCS guide, m = 100 | 5.57 | 8.11 | 7.42 | 7.85 | 8.15 |
| DC-GTM, SCS guide, m = 500 | 5.71 | 8.05 | 8.01 | 8.12 | 8.47 (2 reps) |
| DC-GTM, ASTRAL-III guide, m = 100 | 7.79 | 11.03 | 10.02 | – | – |
| DC-GTM, MRL-FT guide, m = 100 / 500 | – | – | – | – | 11.06 / 9.52 (1 rep) |
| *DC-GTM, true guide (oracle)* | – | – | *6.17* | – | – |

Paired:

| comparison | mean Δ (pp) | W/T/L | p |
|---|---|---|---|
| DC-GTM(SCS, m = 100) vs SCS, all sizes | +4.42 | 0/0/23 | 2.4e-7 |
| DC-GTM(SCS, m = 100) vs ASTRAL-III, 500–2k | +0.02 | 5/2/8 | 0.85 |
| DC-GTM(SCS, m = 500) vs ASTRAL-III, 500–2k | +0.25 | 3/9/3 | 0.58 |
| SCS vs ASTRAL-III, 500–2k | −4.40 | 15/0/0 | 6.1e-5 |
| DC-GTM(ASTRAL, m = 100) vs ASTRAL-III, 500–2k | +2.60 | 0/0/15 | 6.5e-4 |

### 5.3 SMIDGenOG-5500 (replicate 0 = 3,910 taxa)

SCS: 35.8% RF (FN 42.4%, FP 33.6%) in 281 s. DC-GTM with the SCS guide: 45.4% (m = 200) and 46.4% (m = 500). Inside the DC subsets, the ASTRAL-III subset trees miss 1,598 of 3,871 true within-subset splits, against 1,395 for the SCS tree restricted to the same subsets. The source trees themselves have only 7% RF, so the error comes from sparse deep-level information: 5 scaffold trees of 100 taxa constrain the relationships between ~500 clades. ASTRAL-III on the whole replicate: see §6. Replicates 1–4 (SCS and DC only): see `results/og5500.jsonl`.

## 6. Runtime and memory

All jobs ran single-threaded inside ProcessPools, with 3–4 jobs sharing 4 cores, so wall times are inflated by up to ~1.5×. SCS uses numpy threads. The peak is the max RSS of the tool's process.

| n taxa | ASTRAL-III | SCS | MRL-FT | DC-GTM(SCS, m = 100), DC part only | GTM merge only |
|---|---|---|---|---|---|
| 1k (SMIDGen) | 21–24 s / 1.3–1.7 GB | – | 17–27 s / 40 MB | 12 s / 120 MB | 0.2 s |
| 1k (bd) | 36 s / 1.6 GB | 9 s / 0.2 GB | 15 s | 16 s / 120 MB | < 1 s |
| 2k (bd) | 185 s / 3.7 GB | 15 s / 0.3 GB | 43 s | 43 s / 120 MB | < 1 s |
| 3.9k (SMIDGenOG) | **killed (exit 137) after 22 min at 7.7 GB RSS, 8 GB heap, memory shared with other runs** | 281 s / 1.9 GB | – | 108 s (m = 200) | 0.5 s |
| 5k (bd) | killed (exit 137) under shared memory; rerun alone pending (`results/bd5000.jsonl`) | 29 s / 0.7 GB | – | 58 s / 120 MB | < 1 s |
| 10k (bd) | killed (exit 137) after 32 s at 7.8 GB with an 8 GB heap; 13 GB heap alone: `results/bd10000.jsonl` | 65 s / 1.3 GB | 1,033 s / 1.3 GB | 112 s / 120 MB | < 1 s |

DC-GTM's memory is flat (~120 MB at m = 100, ~500 MB at m = 500) and its time is linear in n: about 11 ms per taxon in total for SCS guide + DC, or 3 min at 10k. Extrapolated, 100k taxa would take ~30 min single-threaded. That part of the idea works.

## 7. Verdict: **not promising** (as an accuracy-improving supertree method)

1. **DC-GTM never beat the best method on any dataset.** With realistic guides it lands on the accuracy of its *subset solver* (ASTRAL-III) or of its guide, whichever is worse. It is significantly worse than SCS on the rooted DCM data (+4.4 pp, 0/0/23). On SMIDGen it is significantly worse than ASTRAL-III unless ASTRAL-III itself is the guide, in which case it ties (−0.02 pp, 21/5/14, p = 0.65).
2. **The headroom is small even with a perfect guide.** With the true tree as guide, DC-GTM gets 13.4% / 12.5% on SMIDGen d20/d50 (ASTRAL-III 17.2% / 15.8%) and 6.2% on bd2000 (ASTRAL-III 7.2%). The subset trees, estimated from source trees restricted to the subset, are the bottleneck. Restricting the data loses information that the global analysis uses: on SMIDGen-1000 d20 replicate 0 (MRL-FT guide, m = 200, 7 subsets), the subset trees had 161 within-subset FN, against 147 for the global ASTRAL-III tree restricted to the same subsets.
3. **What does work:** DC-GTM gives ASTRAL-III-level accuracy at a fraction of the memory: about 120 MB versus 3.7 GB at 2k taxa. ASTRAL-III was killed (exit 137) at 3.9k, 5k and 10k with an 8 GB heap while other jobs shared the 15 GB, and its 13 GB solo runs are pending (§6). DC-GTM's cost is linear, about 3 min at 10k. But a scalable method that is better than ASTRAL-III-level accuracy already exists for rooted inputs (SCS).
4. Side findings worth a sentence in a proposal:
   - ASTER's astral4 is a poor supertree method on SMIDGen scaffold data: 32% vs 17% RF, and 29% with `-R`.
   - MRL-FT takes 17 min at 10k taxa because the MRP matrix is about 99% missing.

**If this were turned into a 4-week project anyway**, the only version I would propose changes the subset solver, not the merger:
- (a) Use *overlapping* subsets and an RFS / Exact-RFS-2 merger instead of GTM (blending, per the D&C lecture's open problem). This addresses the information loss in point 2, but each merge is O(n²|X|), so it does not obviously scale.
- (b) Iterate: guide → DC → new guide, and test whether accuracy climbs toward the oracle bound (13.4% vs 17.2% on d20). This is cheap to test: one afternoon.
- (c) Simulate a 100k-taxon SMIDGen-style benchmark, since none exists. This is useful to the field regardless of the method.

**Risks:**
- The oracle bound says the best case is −3 to −4 pp versus ASTRAL-III, and realistic guides captured none of it in this pilot.
- SCS already dominates on rooted inputs.
- FastRFS and BCD, the strongest published baselines, were not run here: one needs an old Bazel build, the other a Java build plus GSCM.
- Supertree accuracy on SMIDGenOG-5500 is limited by the data (~35–45% RF for every method tried), so differences there may be hard to show.

## 8. Reproduce

```
# tools: ASTRAL 5.7.8 (git clone smirarab/ASTRAL), GTM (vlasmirnov/GTM), micromamba env with
#   tree-qmc aster fasttree newick_utils; pip: treeswift dendropy scipy sc-supertree remotezip
cd cs581/supertree/code
python3 run_bench.py smidgen1000 3 astral3,mrlft,dc:astral3:100:astral3,dc:mrlft:100:astral3,dc:mrlft:200:astral3 cases/smidgen1000.tsv
python3 run_bench.py bd_small 3 scs,mrlft,astral3,dc:scs:100:astral3,dc:scs:500:astral3,dc:astral3:100:astral3 cases/bd_small.tsv
python3 analyze.py smidgen1000 > ../results/smidgen1000_summary.md
```

Method names have the form `dc:<guide>:<max subset size>:<subset method>`. `true` as the guide is the oracle.
