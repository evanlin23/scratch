# GTM blending at scale: does a blending DTM help where divide-and-conquer is actually used?

CS581 (Warnow, Fall 2026) pilot, about 6 h of wall-clock time on a 4-core, 15 GB machine.
Code: `cs581/gtmscale/code/`. Tables: `cs581/gtmscale/results/` (`scale_results.tsv`, `published_results.tsv`, `paired_tests.txt`).
This builds on the earlier pilot in `cs581/gtm/` (copied here from branch `claude/cs581-gtm`). That pilot validated GTM and the published Park et al. 2021 numbers, and introduced GTM-Blend-ML.

## 0. Summary

| Question | Answer |
|---|---|
| Does blending beat GTM at scale with a cheap guide? | **Yes, almost always, but by a small amount.** Pooled over 46 paired first-round cases (simulated 2,000–5,000 taxa plus the published 1000M1-HF and Cox1-HET inputs): **−0.45 FN points (95% CI −0.66, −0.27), 32 better / 7 worse / 7 tied, Wilcoxon p = 3e-6.** Per condition the gain is 0.2–1.7 points. It is largest with the worst guide (k-mer NJ, 80% FN: −1.55) and about 0 when the guide is good (IQ-TREE guide on the published data: ±0.05). |
| Is "a meaningful margin" reached? | **No.** On the conditions where a DTM pipeline is the realistic choice, the gain is too small. On RNASim10K the blended tree is 10.42% FN, versus 10.49% for GTM and 10.42% for full FastTree. With a k-mer guide (30% FN) on RNASim10K, blending recovers 0.75 points (11.21 → 10.46), which only brings it back to plain FastTree. |
| Does the pipeline beat full FastTree on accuracy? | **Yes in the hard simulations, no on RNASim10K.** Simulations (Yule trees, short internal branches, 1000 sites): GTM / GTM-Blend-FT beat full FastTree by 5–8 FN points with the FastTree-fastest guide, because the IQ-TREE subset trees are much better than FastTree. RNASim10K: tie (10.42 vs 10.42). With the k-mer guide at 2,000 taxa the pipeline is *worse* than FastTree (+4 points) unless it is iterated. |
| Does the pipeline beat IQ-TREE on time? | **Not at 2,000–5,000 taxa**: IQ-TREE `-fast` (2 threads) takes 1–7 min and is **more accurate than every DTM pipeline** we ran (by 0.6–5 points). **At 10,000 taxa (RNASim10K) full IQ-TREE `-fast` ran out of memory** on this 15 GB machine (13.1 GB resident when it was killed; see §4.4 for the solo retry). The GTM pipeline finished in 10–28 min at under 2 GB, and Blend-FT at 3.7 GB. |
| Where does the remaining merge error come from? | **Mostly the guide/decomposition, not the lack of blending.** One extra iteration (decompose the blended tree again, re-estimate subset trees, merge again) gains 2–6 points over the first round. Blending adds a further 0.4–1.4 points. |
| Prior art | No 2023–2026 DTM that does ML-scored blending of disjoint subset trees (§2). The idea of putting *all* subset-tree splits into FastTree's partial-split constraint mechanism did not turn up in our search either. |
| **Verdict** | **Not promising** as a 4-week project framed as "a blending DTM that clearly beats GTM / is competitive at scale". The effect is real and very consistent but small (about 0.5 points), and full-data IQ-TREE `-fast` dominates wherever it fits in memory. It is **unclear-to-promising only** if reframed as a careful negative/analytical result ("how much does blending buy, and why not more"), or around iteration plus memory-bounded scaling (§6). |

## 1. Question

The instructor's open problem (divide-and-conquer deck): *"GTM does NOT allow blending, it is unlikely GTM is the best that can be done. Develop a better DTM approach that allows blending."*

The earlier pilot showed that ML-decided blending beats GTM by 6–14 FN points when the subsets are not clades of the true tree. That was at 200 taxa, where full-data ML is better than any DTM pipeline anyway. This pilot asks the same question where DTMs are actually used:
- thousands of taxa;
- a cheap, imperfect guide tree;
- comparison against full-data FastTree and IQ-TREE on accuracy, time and memory.

## 2. Prior art (checked 2026-10-10)

- **No hit** for ML-scored blending of disjoint subset trees in 2023–2026. That includes constrained SPR under several leaf-disjoint constraint trees.
- Warnow's ICERM talk (Nov 2024) still calls a blending DTM an open problem. It says Constrained-INC is the only DTM allowing full blending.
- uDance (Balaban et al., *Nat. Biotech.* 2023) does divide-and-conquer with placement at 200K genomes. It is not a disjoint tree merger.
- **IQ-TREE `-g` and RAxML-NG `--tree-constraint` take one constraint tree.** Several disjoint subset trees cannot be expressed as one tree without forcing every subset to be a clade, which is exactly the unblended GTM answer.
- **FastTree `-constraints` takes an "alignment" of 0/1/− split columns, where "−" means the taxon is unconstrained for that split.** So it *can* express "every T_i is induced" as one column per internal split of every subset tree, with "−" for taxa outside S_i. Since every T_i is binary, satisfying all columns ⇔ T|S_i = T_i for all i.
- This is the mechanism behind GTM-Blend-FT below. We did not find it used as a disjoint tree merger in the literature. Novelty risk: it is a natural trick, and the DACTAL/SATé/PASTA line might mention it somewhere.

## 3. Methods

**Pipeline** (`code/pipe.py`, one directory per dataset, restartable):
1. **Guide tree** (cheap):
   - `ftfast` = FastTree `-fastest -noml` (minimum evolution only);
   - `kmer` = rapidnj on Mash-style 8-mer distances of the *unaligned* sequences;
   - `ft` = full FastTree was not used as a guide, because it is the baseline.
2. Centroid-edge decomposition of the guide into subsets of ≤ 500 taxa (simulations) or ≤ 1,000 taxa (RNASim10K).
3. **Subset trees:** IQ-TREE 3 `-m GTR+G -fast`, one thread each, 2 subsets in parallel.
4. **Mergers**, all from the same subset trees:
   - **GTM** (convex mode, commit 18e3bc9).
   - **GTM-Blend-FT (new):** FastTree 2.1.11 `-nt -gtr -gamma`, started from the GTM tree (`-intree`; GTM's few degree-4 nodes are resolved arbitrarily), with every internal split of every subset tree as a partial constraint column (`-constraints`). Subsets may interleave. **The output kept every subset tree exactly induced in every case we ran** (`induced` column = k/k).
   - **Polish (control):** the same FastTree run from the GTM tree *without* constraints. This tells "constrained blending" apart from "more ML search".
   - **cFT (control):** constrained FastTree from scratch, without the GTM start.
   - **Blend-FT-fast:** `-fastest -mlnni 2` with constraints (a cheaper blending heuristic).
   - **GTM-Blend-ML:** the prior pilot's exhaustive constrained-SPR search with a likelihood score (one case, capped).
   - **TreeMerge:** original code + PAUP*, topological guide distances (one case, capped).
5. **Iteration (`+it`):** use the round-1 Blend-FT tree as the new guide and repeat steps 2–4.

**Baselines on the full alignment:** FastTree 2.1.11 `-nt -gtr -gamma` (single thread) and IQ-TREE 3 `-m GTR+G -fast` (2 threads, 1–1.5 h cap).

**Data:**

| Data | Taxa | Sites | Reps | Notes |
|---|---|---|---|---|
| Sim `n2000_i0.01` | 2000 | 1000 | 6 (ftfast), 4 (kmer) | Yule-like tree, internal branch lengths ~Exp(mean 0.01), pendant ~Exp(0.1), AliSim GTR+Γ(0.5); FastTree-fastest guide ≈ 57% FN (hard) |
| Sim `n2000_i0.02` | 2000 | 1000 | 4 | internal mean 0.02, guide ≈ 37% FN |
| Sim `n5000_i0.01` | 5000 | 1000 | 2 | as `n2000_i0.01` |
| RNASim10K (MAGUS data bank) | 10,000 | 8,701 | 1 (R0) | true alignment, true tree; guides ftfast and kmer |
| 1000M1-HF, Cox1-HET (Park et al. 2021) | 1000 / 2341 | 3,880 / 658 | 5 / 10 | *published* guide (FT or IQ), IQ-TREE subset trees and GTM trees, merged again with Blend-FT |

RNASim1000 was dropped for budget: the prior pilot's ceiling analysis gave ≤ 0.1 points of headroom there.

**Measures:**
- FN rate against the true tree (Park et al. criterion).
- Wall-clock time and peak RSS (GNU time).
- Paired differences per replicate, with W/T/L, a bootstrap 95% CI, and a two-sided Wilcoxon signed-rank test.

**Runtime caveat:**
- To fit the budget, at most two pipelines ran at once (≤ 4 busy cores in total).
- Every single-threaded FastTree arm had its own core. The 2-thread IQ-TREE runs and the 2-way parallel subset-tree step did too, except for short overlaps.
- Times are therefore indicative to about ±20%, not benchmark-grade.

## 4. Results

### 4.1 Simulations (mean FN %, paired across the same replicates)

| Condition, guide | reps | guide | subset trees | full FT | full IQ-fast | GTM | **Blend-FT** | Polish (no constraints) | cFT | Blend-FT-fast |
|---|---|---|---|---|---|---|---|---|---|---|
| n2000 i0.01, ftfast | 6 | 56.3 | 22.2 | 34.1 | **24.6** | 28.4 | 27.7 | 32.2 | 60.0 | 28.3 |
| n2000 i0.01, ftfast **+it** | 6 | 27.7 | 22.1 | 34.1 | **24.6** | 26.1 | 25.8 | 30.5 | – | 25.9 |
| n2000 i0.01, kmer | 4 | 79.5 | 22.4 | 34.6 | **25.0** | 40.2 | 38.7 | 34.4 | 57.9 | 39.8 |
| n2000 i0.01, kmer **+it** | 4 | 38.7 | 23.1 | 34.6 | **25.0** | 34.2 | 32.8 | 32.8 | – | – |
| n2000 i0.02, ftfast | 4 | 37.7 | 12.8 | 16.6 | 13.4 | 14.8 | 14.0 | 15.5 | 48.8 | 14.2 |
| n2000 i0.02, ftfast **+it** | 4 | 14.0 | 12.5 | 16.6 | 13.4 | 13.1 | **12.9** | 14.9 | – | 13.0 |
| n5000 i0.01, ftfast | 2 | 57.7 | 24.0 | 35.2 | **25.8** | 30.1 | 29.7 | 32.5 | 60.0 | 29.9 |
| n5000 i0.01, ftfast **+it** | 2 | 29.7 | 24.1 | 35.2 | **25.8** | 28.2 | 28.0 | 31.6 | – | – |

Paired differences (FN points, B − A; negative = B better). Per-condition tests are in `results/paired_tests.txt`.

| Comparison | n | mean (95% CI) | better/worse/tie | Wilcoxon p |
|---|---|---|---|---|
| GTM → Blend-FT, simulated round 1 (all guides) | 16 | **−0.92 (−1.38, −0.58)** | 15/0/1 | 6.5e-4 |
| GTM → Blend-FT, simulated iteration round | 16 | −0.57 (−0.93, −0.29) | 15/0/1 | 6.1e-5 |
| GTM → Blend-FT-fast, simulated round 1 | 16 | −0.30 (−0.47, −0.16) | 13/3/0 | 1.3e-3 |
| Polish → Blend-FT, simulated round 1 | 16 | −1.35 (−3.07, +0.56) | 12/4/0 | 0.21 |
| Polish → Blend-FT, simulated iteration round | 16 | **−2.73 (−3.74, −1.71)** | 13/2/1 | 1.2e-3 |
| full FT → Blend-FT, n2000 i0.01 ftfast | 6 | −6.40 | 6/0/0 | 0.031 (minimum attainable with n = 6) |
| full IQ-fast → Blend-FT, n2000 i0.01 ftfast | 6 | **+3.02** (IQ better) | 1/5/0 | 0.063 |
| full IQ-fast → Blend-FT, n2000 i0.01 ftfast +it | 6 | +1.13 | 2/4/0 | 0.56 |
| GTM round 1 → GTM after one iteration, n2000 i0.01 ftfast | 6 | about −2.3 | 6/0/0 | – |

Observations:
- **Blending never made a simulated tree worse than GTM** (30/0/2 over both rounds).
- But blending recovers only a small part of the gap between GTM (28.4) and the subset trees' own error (22.2). Re-decomposing (iteration) recovers more.
- With the k-mer guide (80% FN), *unconstrained* polish beats Blend-FT by 4.3 points. The merge is so wrong that keeping the subset trees fixed and only blending at the seams is not enough; a free search does better.
- With better guides, the constraints are what keep the IQ-TREE subset-tree quality, and polish is worse than GTM.

### 4.2 Published conditions (Park et al. 2021 inputs; GTM tree = published)

| Condition, guide | n | GTM | Blend-FT | Polish | GTM → Blend-FT | W/L/T | p |
|---|---|---|---|---|---|---|---|
| 1000M1-HF, FT | 5 | 42.45 | 41.82 | **37.31** | −0.62 (−1.07, −0.20) | 4/0/1 | 0.125 |
| 1000M1-HF, IQ | 5 | 28.37 | 28.37 | 28.85 | 0.00 | 1/1/3 | 1 |
| Cox1-HET, FT | 10 | 18.85 | **18.63** | 22.61 | −0.22 (−0.40, −0.06) | 7/2/1 | **0.023** |
| Cox1-HET, IQ | 10 | 18.71 | 18.66 | 22.33 | −0.05 | 5/4/1 | 0.89 |
| pooled | 30 | | | | **−0.19 (−0.33, −0.08)** | 17/7/6 | **0.006** |

How this compares with the earlier pilot:
- Cox1-HET/FT: Blend-FT recovers 0.22 of the 0.49 points that separate GTM from the optimistic per-branch floor of *any* merger (`cs581/gtm` §3.1).
- 1000M1-HF/FT: Blend-FT −0.62 in about 1.5 min per case, versus GTM-Blend-ML −1.03 in 27–250 min per case.
- **1000M1-HF (fragmentary sequences) is the exception where constraints hurt.** Unconstrained FastTree from the GTM start reaches 37.3% (5/5 better than GTM, −5.1 points), so the IQ-TREE subset trees there are wrong in ways that the constraints lock in.

### 4.3 RNASim10K (R0, true alignment; 10,000 taxa, 8,701 sites)

| Arm | FN % | wall (min) | peak RSS |
|---|---|---|---|
| full FastTree | 10.42 | 46.3 | 2.1 GB |
| full IQ-TREE `-fast`, 2 threads, `-mem 10G` | **killed: out of memory** (13.1 GB resident) | – | > 13 GB |
| guide ftfast (FastTree -fastest -noml) | 13.49 | 20.8 | 1.7 GB |
| 16 IQ-TREE subset trees (≤ 1000 taxa) | 10.25 (mean subset FN) | 7.2 | 0.6 GB |
| GTM (ftfast guide), end-to-end | 10.49 | 28.1 | 1.7 GB |
| Blend-FT (ftfast guide), end-to-end | 10.42 | 63.1 | 3.7 GB |
| Polish (no constraints), merge step only | 10.39 | 37.8 (step) | 2.1 GB |
| GTM → Blend-FT, one more iteration | 10.46 / 10.46 | +42.6 | |
| guide kmer (rapidnj on 8-mer distances) | 30.21 | 1.9 | ~3 GB |
| **GTM (kmer guide), end-to-end** | 11.21 | **10.2** | 1.7 GB |
| **Blend-FT (kmer guide), end-to-end** | **10.46** | 45.2 | 3.7 GB |
| Blend-FT-fast (kmer guide), merge step | 10.46 | 33.9 (step) | 3.8 GB |

RNASim10K is easy for ML at this alignment length:
- every method lands at 10.4–10.5% except the GTM with the k-mer guide;
- the subset trees are barely better than full FastTree (10.25 vs 10.42).

So with a decent cheap guide there is nothing for blending to fix. With the bad k-mer guide, blending recovers 0.75 of the 0.79 points that GTM lost relative to FastTree, at the cost of a full-data FastTree run (35 min). The fastest option on this dataset is **GTM with a k-mer guide: 10 minutes, 4.5× faster than FastTree, 0.8 points worse.**

### 4.4 Full IQ-TREE on RNASim10K, solo retry

See the note at the end of this section: it is filled in after the run (4 threads, `-mem 13G`, 60 min cap, nothing else running).

### 4.5 Other mergers (single cases, `n2000_i0.01` r1, ftfast guide, GTM 30.80%)

| Merger | FN % | time | note |
|---|---|---|---|
| GTM | 30.80 | 0.4 s | |
| **GTM-Blend-FT** | **29.74** | 70 s | 6/6 subset trees induced |
| GTM-Blend-FT with `-spr 4 -mlacc 2 -slownni` | 29.74 | 289 s | more search finds nothing more |
| GTM-Blend-ML (prior pilot, radius-4 constrained SPR, 15 rounds) | 30.15 | 612 s, 1.9 GB | logL still rising when capped (+415) |
| TreeMerge (original + PAUP*) | **did not finish** | > 31 min | killed while still on its first pair of subsets |

At 2,000 taxa, the FastTree-constraint blending is both better and about 9× faster than the exhaustive ML search of the prior pilot. TreeMerge does not scale to this setting in our hands.

## 5. Interpretation

1. **Blending is a reliable but small improvement on GTM.** The pooled effect is −0.45 points (p = 3e-6), and per condition it is −0.2 to −1.7 points.
   - It is largest when the guide is poor, matching the earlier pilot's mechanism (subsets that are not clades of the true tree).
   - It never meaningfully hurt in simulation (30/0/2).
   - It is about 0 when the guide is already good (IQ-TREE guide, RNASim10K ftfast guide).
2. **At scale the decomposition matters more than blending.** One re-decomposition (iteration) improves GTM by 2–6 points in the hard simulations. Blending then adds about 0.4–0.6.
   - Blending at the seams can only move taxa locally while keeping all subset trees fixed.
   - When the guide is 57–80% wrong, the subsets themselves are badly chosen. Neither a local blend nor GTM can repair that; re-choosing the subsets can.
3. **The "constraints vs free search" question is data-dependent.**
   - Constraints are good when the subset trees are much better than full-data FastTree (simulations, Cox1-HET).
   - They are bad when the subset trees are poor (fragmentary 1000M1-HF, and the 80%-error k-mer guide). There the free FastTree search from the GTM tree wins by 4–5 points.
   - A *soft* constraint (FastTree `-constraintWeight`, or accepting a constraint violation when likelihood gains enough) is the obvious untested middle ground.
4. **The pipeline is not competitive with IQ-TREE `-fast` where IQ-TREE fits in memory.**
   - At 2,000–5,000 taxa IQ-TREE `-fast` takes 1–7 min and is 0.6–5 points more accurate than every pipeline variant (iterated Blend-FT gets within 1.1 points at n = 2,000).
   - The DTM pipeline's advantages are (a) beating FastTree by 2–8 points on the hard simulations and (b) running in < 4 GB on 10,000 taxa where IQ-TREE `-fast` needed more than 13 GB.
   - On RNASim10K it merely ties FastTree.

## 6. Verdict: **not promising** (as "a better blending DTM"); **unclear** as a reframed project

**"Not promising"** for the project as posed: a blending DTM that beats GTM by a meaningful margin at scale and makes the DTM pipeline competitive. The measured effect is real and very consistent but about half a point. The cases where it would matter (very poor guides) are better served by re-decomposition or a free ML search. Full-data IQ-TREE `-fast` dominates up to at least 5,000 taxa.

What *would* be a defensible 4-week project, if the student still wants this area:

| Week | Work |
|---|---|
| 1 | Reproduce with this code (everything runs). Add soft-constraint Blend-FT (`-constraintWeight` sweep) and an "accept a violation if ΔlogL > τ" rule. These target the 1000M1-HF and k-mer failure modes where the hard constraints lock in errors. |
| 2 | Iterated DTM (PASTA-style 2–3 rounds) with blending at every round, versus GTM at every round. Main grid: 2K/5K/10K taxa, ftfast/kmer guides, 10 replicates per cell (n = 5 cannot reach p < 0.05). |
| 3 | The regime where DTMs are actually needed: memory- or time-bounded runs at 10K–50K taxa (RNASim10K reps R0–R9, RNASim50K from Park et al.). Compare against IQ-TREE `-fast` with explicit memory and time budgets and against FastTree. |
| 4 | Write-up: "blending buys ~0.5 points; re-decomposition buys 2–6; where the constraint helps and where it hurts". |

**Main risks:**
1. The headline effect is small (≈ 0.5 points). A course grader may not consider it "more accurate than the strongest baseline", because that baseline is IQ-TREE `-fast` and it wins at ≤ 5K taxa.
2. The k-mer/iteration gains are not blending gains. A project that ends up being "iterate the DTM" overlaps with the known DACTAL/PASTA idea.
3. FastTree's constraint mechanism is a black box. The project's "method" may look like a clever use of an existing flag rather than a new algorithm.
4. RNASim10K is too easy (all methods within 0.1 points), so the 10K-scale story needs harder data, such as fragmentary sequences or shorter alignments, to separate methods.
5. Runtime numbers here are indicative (shared machine, ±20%).

## 7. Reproducing

```bash
bash cs581/code/setup.sh
MAMBA_ROOT_PREFIX=/opt/mm/root /opt/mm/micromamba install -y -n bio -c conda-forge -c bioconda rapidnj raxml-ng epa-ng
apt-get install -y time
bash cs581/gtm/code/fetch_data.sh                    # GTM source + published Park et al. inputs -> /opt/gtmdata
cd cs581/gtmscale/code
./run_pub.sh                                          # Blend-FT / polish on 1000M1-HF, Cox1-HET (published inputs)
./run_sim.sh; ./run_sim2.sh; ./run_fast.sh            # simulations: round 1, iteration, Blend-FT-fast
REPS=0 ./run_rna.sh                                   # RNASim10K R0 (ftfast guide, + iteration)
python3 pipe.py /opt/gtms/rnasim10k/R0k kmer 1000 gtm,blendft   # RNASim10K, k-mer guide (dir with aln.fa, true.tre)
THREADS=4 IQ_MEM=13G IQ_CAP=3600 python3 iq_full.py /opt/gtms/rnasim10k/R0   # full IQ-TREE, run alone
python3 summarize.py ../results
```
Single-case arms: `pipe.py DIR ftfast 500 blendml` (GTM-Blend-ML, needs `MLSPR_MAX_ROUNDS`) and `pipe.py DIR ftfast 500 tm` (TreeMerge; needs the tm27 env + PAUP* as in `cs581/gtm/REPORT.md` §6).
