# Does the GCM evidence recipe improve ML trees? (pre-registered)

CS581 (UIUC), Fall 2026. Pre-registration: `PREREG.md`, pushed before any tree was computed. Code: `code/`.
Rows: `results*/aln.jsonl` (SP scores), `results*/trees.jsonl` (FastTree), `results/iqtrees.jsonl` (IQ-TREE),
`results/colana.jsonl` (column analysis), `results*/magus.jsonl` (MAGUS runs). All tables come from `code/summarize.py`
(`results/summary.md`).

## Verdict (numbers first)

16 AliSim protein datasets (8 SIMHIGH, 8 SIMMOD; 1,000 taxa), one MAGUS draw each. Every method is a merge-only
variant on MAGUS's own subsets and backbones. Δ = method − MAGUS. RF values are FastTree `-lg -gamma` normalized RF
to the true tree, in %. W/T/L uses a 0.1-point tie band. p = two-sided Wilcoxon signed-rank.

| | n | Δ SP error | ΔSPFN / ΔSPFP | **Δ nRF** | W/T/L | p |
|---|---|---|---|---|---|---|
| **primary: recipe `wsoft0.03:linsi&fftns2#es4`, all** | 16 | −5.69 | −10.79 / −0.58 | **−0.19** | 8/2/6 | **0.63** |
| recipe, SIMHIGH only | 8 | −7.71 | −14.76 / −0.66 | −0.49 | 5/0/3 | 0.46 |
| recipe, SIMMOD only | 8 | −3.66 | −6.82 / −0.51 | +0.11 | 3/2/3 | 0.58 |
| `linsi#es3`, all | 16 | −5.04 | −10.53 / +0.44 | −0.07 | 8/5/3 | 0.26 |
| hard filter `linsi&fftns2-op3`, all | 16 | −2.46 | −4.48 / −0.44 | −0.21 | 8/2/6 | 0.35 |
| hard filter, SIMHIGH | 8 | −3.10 | −5.57 / −0.62 | −0.40 | 5/1/2 | 0.30 |
| *reference: true alignment vs MAGUS*, SIMHIGH | 8 | −24.2 | −31.6 / −16.8 | **−3.51** | 8/0/0 | 0.008 |
| *reference: true alignment vs MAGUS*, SIMMOD | 8 | −10.7 | −14.3 / −7.1 | −0.16 | 4/1/3 | 0.71 |
| IQ-TREE 3 `LG+G4 --fast`, recipe − MAGUS (SIMHIGH R1-R3) | 3 | | | −0.47 | 2/1/0 | 0.25 |

**Does the recipe improve ML tree accuracy? No, not detectably.** The pre-registered primary test fails
(−0.19 nRF points, p = 0.63). None of the secondary tests is significant either. The recipe removes about half of
MAGUS's SPFN on SIMHIGH (−14.8 of 31.6 SPFN points), yet it closes only 0.49 of the 3.51 RF points that separate
MAGUS's tree from the true-alignment tree (14%). Even that 0.49 is indistinguishable from zero. On SIMMOD there is
nothing to gain: FastTree on MAGUS's alignment is already as accurate as on the true alignment (Δ −0.16, p = 0.71).
This replicates protcons's null result (Δ −0.11) with the stronger recipe, twice as many datasets, and new
replicates R5-R8.

**Why not?** Three findings:

1. **Tree changes from these alignment changes are mostly noise.**
   - es3 and the recipe produce alignments within about 1 SP point of each other (SIMHIGH: −6.70 vs −7.71). Yet
     their FastTree trees differ by 0.64 RF points on average (SD 0.69) on SIMHIGH.
   - The recipe's own Δ RF has SD 1.22 on SIMHIGH. That SD sets the scale.
   - Detecting the observed −0.49 on SIMHIGH with 80% power would need about 50 datasets. Detecting the pooled −0.19
     would need about 190.
   - Tree error does follow alignment error *across* datasets: room (RF MAGUS − RF true) vs MAGUS SPFN, Spearman
     ρ = 0.83, p = 8e-5, n = 16. But tree error is a weak, noisy function of SP error *within* a dataset. Over 48
     (dataset × method) pairs, Δ RF tracks Δ SPFN only weakly (ρ = +0.32, p = 0.027) and does not track Δ SPFP
     (ρ = +0.03).
2. **The recovered pairs are in columns that matter, so "they are in irrelevant columns" is not the explanation.**
   - 97-99% of the recipe's recovered true pairs lie in the ~300 dense (≥ 50% of taxa), parsimony-informative
     columns, i.e. the root-derived sites.
   - At the residue level the recipe halves the number of residues split away from their true column (SIMHIGH
     R1-R5: 16-21% → 7-12%). It also slightly lowers misplaced residues (9-12% → 8-11%).
   - The changes are large, and they are in the tree-informative core.
3. **What the recipe gains in correct pairs, it gives back through the pairs it gets wrong.** This comes from an
   oracle diagnostic, `code/refine.py`. It builds split(X), the common refinement of the true alignment and X:
   - split(X) keeps exactly X's true-positive pairs, so it has X's SPFN and SPFP = 0;
   - "FP cost" = RF(X) − RF(split(X)).

   Results on SIMHIGH R1-R5:
   - **Over-splitting is costly.** split(MAGUS) is 2.6 RF points worse than the true alignment (mean 8.81 vs 6.18).
     MAGUS's tree error is mostly SPFN-type, not FP-type.
   - **The recipe's correct pairs would help.** split(recipe) is better than split(MAGUS) by 0.79 RF points on
     average (range −3.9 to +2.6; 4 of 5 better).
   - **The recipe's wrong pairs hurt more.** Their FP cost is +1.85 RF points vs +1.10 for MAGUS's, although the
     recipe has *fewer* FP pairs (ΔSPFP −0.66) and fewer misplaced residues.
   - **The net change is small.** The two effects nearly cancel, giving −0.04 on these five datasets. My reading
     (not proven, n = 5, noisy) is that deleting low-support edges lets MCL merge larger clusters. When such a merge
     is wrong, it puts whole groups of residues from different clades into one column, and ML trees weigh that more
     than the same number of scattered wrong pairs.

**What would make it improve trees?**
- **A larger SPFN gain.** On these data trees respond to alignment error only when a lot of it goes away: the true
  alignment removes ~32 SPFN and ~17 SPFP points for 3.5 RF points. The recipe removes ~15 SPFN points and almost
  no SPFP, which predicts at most ~1 RF point even if the relationship were linear.
- **FP control on the merged clusters.** The oracle split(recipe) gains up to 3.9 RF points on single datasets.
  Something that keeps the recipe's recovered pairs but rejects its wrong merges would be the most direct route,
  e.g. a support or consistency check per MCL cluster, or masking the columns the merge created.
- **Harder data.** Room exists only at SIMHIGH-level divergence.
- **Statistical power.** About 50 SIMHIGH-like datasets, or several MAGUS draws per dataset, to resolve
  effects of ~0.5 RF points.
- **A tree method that uses the alignment more directly,** or pairing with tree-aware masking. Not tested here.

The alignment-level result replicates: the recipe lowers SP error on all 16 datasets (−2.4 to −10.1 points), 16/16
wins. Its overhead is small: 21 s for FFT-NS-2 plus the merge, vs 44 s for MAGUS's own merge, because the
edge-support cut shrinks the graph.

## 1. Data

AliSim (IQ-TREE 3.1.4), LG+G4, root length 300 aa, indels 0.05/0.05 with POW(1.7, 40) lengths, 1,000 taxa.
Yule-Harding tree `iqtree3 -r 1000 -rlen 0.001 MEAN 0.8`, MEAN = 0.06 (SIMMOD) / 0.10 (SIMHIGH). Tree seed
100r+7, sequence seed 100r+13 (`code/gen.sh`).
- R1-R4 are protbench's replicates, regenerated. SIMMOD_R1 has 7,348 columns, identical to protbench, and SIMHIGH_R2's
  true-alignment FastTree tree reproduces protbench's 5.22% RF.
- R5-R8 are new replicates with the same parameters.

| level | MAGUS SP error % | true-alignment FastTree nRF % | MAGUS FastTree nRF % | room |
|---|---|---|---|---|
| SIMHIGH (8) | 20.6-27.6 | 5.2-9.7 | 8.6-12.1 | +2.1 to +5.8 |
| SIMMOD (8) | 9.6-12.5 | 5.5-8.5 | 5.9-8.3 | −0.6 to +1.2 |

## 2. Methods

**Datasets split across machines.** To finish in time, the orchestrating session split the 16 datasets across four
machines (4 cores / 15 GB each), all running this branch's code:
- this session: SIMHIGH R1-R5 (`results/`);
- helper h1: SIMHIGH R6-R8 (`results_h1/`, branch `claude/cs581-gcmtrees-h1`);
- helper h2: SIMMOD R1-R4 (`results_h2/`);
- helper h3: SIMMOD R5-R8 (`results_h3/`).

Wall times are per machine. On this machine MAGUS shared the CPU with FastTree lanes, so its wall times
(1,090-3,082 s) are inflated. Uncontended MAGUS took 544-1,273 s.

**Per dataset** (`code/run.sh`):
- one MAGUS draw with the paper's flags (`gcmx.bbtool_bench`: 25 subsets, 10 backbones × 200, MCL, minclusters,
  4 threads);
- `pc.py rep` builds the merge-only replicate;
- `gcmgen/code/gg.py run` runs the four variants on the same subsets and backbone sequence sets: `linsi` (MAGUS's
  own merge, the baseline), `wsoft0.03:linsi&fftns2#es4`, `linsi#es3`, `linsi&fftns2-op3`.
- The `linsi` control reproduced MAGUS's own SP error on all 16 datasets: identical on 13, within 0.07 points on the
  other 3.
- SP: FastSP. Error = (SPFN+SPFP)/2.

**Trees** (`code/trees.sh` → `protbench/code/trees.py`):
- FastTree 2.1 `-lg -gamma`, single-threaded, on the true alignment and on the four estimates;
- normalized FN/FP/RF against the true tree with DendroPy. All trees are binary, so FN = FP = RF.
- IQ-TREE 3 `-m LG+G4 --fast -T 1 -seed 1` (`code/iqtrees.py`) on SIMHIGH only, as time allowed.

**Why-not diagnostics:**
- `code/colana.py`: where recovered pairs fall, by true-column class, plus residue-level split/misplaced counts.
- `code/refine.py` → `code/extra.sh`: the zero-SPFP refinements, with FastTree on split(MAGUS) and split(recipe).

## 3. Results per dataset (FastTree nRF %, SP error points)

| dataset | RF true | RF MAGUS | room | RF recipe | RF es3 | RF hard | Δ recipe | Δ es3 | Δ hard | SP err MAGUS | Δ SP recipe | ΔSPFN / ΔSPFP recipe |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SIMHIGH_R1 | 6.12 | 11.94 | +5.82 | 12.34 | 11.63 | 12.34 | +0.40 | −0.31 | +0.40 | 25.32 | −5.32 | −9.66 / −0.98 |
| SIMHIGH_R2 | 5.22 | 9.33 | +4.11 | 7.52 | 8.73 | 7.82 | −1.81 | −0.60 | −1.51 | 22.91 | −8.63 | −17.65 / +0.40 |
| SIMHIGH_R3 | 6.42 | 10.03 | +3.61 | 9.23 | 9.83 | 8.93 | −0.80 | −0.20 | −1.10 | 27.55 | −9.47 | −17.97 / −0.97 |
| SIMHIGH_R4 | 6.62 | 9.63 | +3.01 | 10.53 | 11.94 | 8.63 | +0.90 | +2.31 | −1.00 | 22.95 | −5.52 | −10.29 / −0.74 |
| SIMHIGH_R5 | 6.52 | 8.63 | +2.11 | 9.73 | 9.53 | 9.13 | +1.10 | +0.90 | +0.50 | 20.62 | −6.63 | −13.27 / +0.01 |
| SIMHIGH_R6 | 9.73 | 12.14 | +2.41 | 11.43 | 11.74 | 11.84 | −0.71 | −0.40 | −0.30 | 20.63 | −7.83 | −15.56 / −0.10 |
| SIMHIGH_R7 | 6.62 | 10.63 | +4.01 | 8.32 | 8.53 | 10.63 | −2.31 | −2.10 | +0.00 | 27.28 | −8.14 | −14.51 / −1.77 |
| SIMHIGH_R8 | 8.22 | 11.23 | +3.01 | 10.53 | 11.03 | 11.03 | −0.70 | −0.20 | −0.20 | 26.19 | −10.14 | −19.18 / −1.10 |
| SIMMOD_R1 | 6.52 | 5.92 | −0.60 | 6.62 | 5.82 | 6.22 | +0.70 | −0.10 | +0.30 | 12.50 | −2.39 | −4.11 / −0.68 |
| SIMMOD_R2 | 6.52 | 7.02 | +0.50 | 7.02 | 6.52 | 6.92 | +0.00 | −0.50 | −0.10 | 10.15 | −4.52 | −8.74 / −0.30 |
| SIMMOD_R3 | 6.12 | 7.32 | +1.20 | 7.12 | 7.22 | 7.02 | −0.20 | −0.10 | −0.30 | 11.32 | −4.22 | −7.80 / −0.63 |
| SIMMOD_R4 | 5.52 | 6.22 | +0.70 | 6.92 | 6.82 | 6.62 | +0.70 | +0.60 | +0.40 | 10.69 | −3.82 | −7.52 / −0.12 |
| SIMMOD_R5 | 8.53 | 8.02 | −0.51 | 8.43 | 7.92 | 7.52 | +0.41 | −0.10 | −0.50 | 9.88 | −4.41 | −8.65 / −0.18 |
| SIMMOD_R6 | 7.92 | 8.32 | +0.40 | 7.82 | 8.02 | 8.63 | −0.50 | −0.30 | +0.31 | 9.81 | −2.52 | −5.07 / +0.03 |
| SIMMOD_R7 | 6.82 | 6.92 | +0.10 | 6.62 | 6.82 | 6.52 | −0.30 | −0.10 | −0.40 | 11.49 | −4.72 | −8.19 / −1.25 |
| SIMMOD_R8 | 8.53 | 8.02 | −0.51 | 8.12 | 8.02 | 8.22 | +0.10 | +0.00 | +0.20 | 9.55 | −2.71 | −4.49 / −0.93 |

Paired summaries, including FN/FP and SP breakdowns per level: `results/summary.md`. FN and FP deltas equal the RF
deltas because all trees are fully resolved.

## 4. Why-not analyses

### 4.1 Zero-SPFP refinements (oracle; SIMHIGH R1-R5)

| dataset | true | split(MAGUS) | MAGUS | split(recipe) | recipe | FP cost MAGUS | FP cost recipe | split(recipe) − split(MAGUS) |
|---|---|---|---|---|---|---|---|---|
| SIMHIGH_R1 | 6.12 | 12.24 | 11.94 | 8.32 | 12.34 | −0.30 | +4.02 | −3.92 |
| SIMHIGH_R2 | 5.22 | 6.32 | 9.33 | 8.93 | 7.52 | +3.01 | −1.41 | +2.61 |
| SIMHIGH_R3 | 6.42 | 8.93 | 10.03 | 7.32 | 9.23 | +1.10 | +1.91 | −1.61 |
| SIMHIGH_R4 | 6.62 | 8.43 | 9.63 | 8.22 | 10.53 | +1.20 | +2.31 | −0.21 |
| SIMHIGH_R5 | 6.52 | 8.12 | 8.63 | 7.32 | 9.73 | +0.51 | +2.41 | −0.80 |
| mean | 6.18 | 8.81 | 9.91 | 8.02 | 9.87 | +1.10 | +1.85 | −0.79 |

The recipe's FP cost exceeds MAGUS's on 4 of 5 datasets.

The per-dataset values swing by ±3-4 RF points between alignments that differ only in their FP pairs. That is
another view of the noise floor in finding 1.

### 4.2 Where the alignment changes land (SIMHIGH R1-R5; `results/colana.jsonl`)

| dataset | true pairs in dense cols | MAGUS recall dense / gappy | recipe gain dense / gappy (% of true pairs) | split residues MAGUS / recipe / es3 | misplaced residues MAGUS / recipe / es3 |
|---|---|---|---|---|---|
| SIMHIGH_R1 | 97.4% | 68.7% / 73.6% | +9.32 / +0.34 | 18.8% / 12.2% / 12.3% | 12.2% / 10.9% / 12.1% |
| SIMHIGH_R2 | 98.6% | 68.7% / 70.7% | +17.43 / +0.22 | 19.2% / 7.4% / 8.9% | 9.6% / 8.9% / 9.4% |
| SIMHIGH_R3 | 97.7% | 63.3% / 73.7% | +17.64 / +0.33 | 21.4% / 9.7% / 10.3% | 11.2% / 10.2% / 11.5% |
| SIMHIGH_R4 | 98.8% | 72.3% / 75.9% | +10.18 / +0.11 | 16.2% / 9.2% / 9.4% | 10.9% / 10.1% / 10.8% |
| SIMHIGH_R5 | 98.0% | 72.8% / 75.3% | +13.04 / +0.23 | 16.5% / 7.4% / 8.1% | 8.9% / 8.3% / 8.8% |

"Dense" means true columns holding ≥ 50% of the taxa. Nearly all of them are parsimony-informative, and the gain in
parsimony-informative columns equals the dense-column gain to within 0.3 points. Pair counts are dominated by dense
columns by construction (C(n, 2)), which is why the residue-level columns are given too.

### 4.3 Noise

- FastTree is deterministic: re-running it with the sequences shuffled gave an identical tree on SIMHIGH_R1. So
  search noise cannot be measured that way.
- Instead, nearly equivalent alignments serve as the noise yardstick. es3 vs recipe (about 1 SP point apart) differ
  by SD 0.69 RF points on SIMHIGH and 0.37 on SIMMOD.
- Recipe − MAGUS has SD 1.22 (SIMHIGH) and 0.94 (all).
- A second MAGUS draw per dataset was not run (time).

### 4.4 IQ-TREE (SIMHIGH; as far as time allowed)

| dataset | true | MAGUS | recipe | Δ recipe |
|---|---|---|---|---|
| SIMHIGH_R1 | 6.02 | 13.34 | 12.84 | −0.50 |
| SIMHIGH_R2 | 5.12 | 9.83 | 9.03 | −0.80 |
| SIMHIGH_R3 | — | 10.13 | 10.03 | −0.10 |

IQ-TREE agrees with FastTree on the true-alignment trees: 6.02 vs 6.12 and 5.12 vs 5.22. It puts MAGUS somewhat
higher (R1 13.34 vs 11.94). n = 3 cannot support a test (mean −0.47, 2/1/0, p = 0.25). On the same three datasets FastTree gave
+0.40, −1.81 and −0.80 (mean −0.74). The per-dataset signs disagree between the two tree methods on R1, another sign
of tree-estimation noise at this effect size. Each IQ-TREE run took about
30 min and needed ~4.7 GB of RAM. Four concurrent runs were killed by the OOM killer, which is why only R1-R3 were
done.

## 5. Runtime

| step | mean wall |
|---|---|
| MAGUS end to end, uncontended (helpers; 4 threads) | SIMMOD 544-842 s, SIMHIGH 1,052-1,273 s |
| MAGUS merge only (`linsi`) | 44 s |
| recipe: FFT-NS-2 on 10 backbones + merge | 21 s |
| `linsi#es3` merge | 10 s |
| hard filter: FFT-NS-2 `--op 3` + merge | 30 s |
| FastTree `-lg -gamma`, 1 thread, contended | true 990 s, MAGUS 589 s, recipe 559 s |
| IQ-TREE `--fast`, 1 thread | ~1,750 s |

## 6. Deviations from the pre-registration

- **Split across four machines** at the orchestrator's request; same code. All 16 planned datasets were completed.
- **IQ-TREE** only on SIMHIGH R1-R3 (time and memory). No second MAGUS draw, no harder level.
- **Diagnostics added after the first trees were seen,** so they are exploratory:
  - the oracle split(X) refinements;
  - the residue-level split/misplaced counts;
  - the shuffled-input FastTree control (it showed FastTree is order-invariant and was dropped).
- Commits that imported the earlier sessions' code (`bbevidence`, `protcons`, `gcmgen`, `protbench` code
  directories, plus the newer `gcmx/bbtool_bench.py`) were made with the pre-registration commit.
