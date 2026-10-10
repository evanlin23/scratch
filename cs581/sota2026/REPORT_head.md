# Is MAGUS still state of the art, and is "MAGUS with a better base method" a good 4-week project?

*CS581 project-idea pilot, 2026-10-10. Branch `claude/cs581-sota2026`. Code in `code/`, raw results in `results/`.
About 4.5 hours of wall-clock on one 4-core / 15 GB cloud machine. Every number below was measured in this session,
except where it says "published", which means FastSP rescoring of the MAGUS paper's own output alignments
(`../experiments/validation/published_scores.jsonl`).*

## 1. Questions

1. On MAGUS's own benchmarks, is MAGUS (Smirnov & Warnow 2021) still the best accuracy/runtime trade-off against
   aligners published since? Candidates: TWILIGHT (2025), FAMSA 2, MUSCLE5/Super5, Clustal Omega, MAFFT
   (`--auto`, PartTree, L-INS-i) and learnMSA2.
2. On MAGUS-sized subsets (about 40 and 200 sequences), is any of them more accurate than MAFFT L-INS-i? This is the
   instructor's suggested first step for "study MAGUS with different base methods".
3. (Added, because the harness made it cheap.) If a candidate swaps in for L-INS-i inside MAGUS, does MAGUS's
   final alignment improve?

## 2. Prior art (details and DOIs: `prior_art.md`)

- **No independent 2023–2026 benchmark exists.** No paper runs TWILIGHT, FAMSA2, Muscle5, Clustal Omega and the MAFFT
  modes against MAGUS, PASTA or UPP2 on MAGUS's nucleotide suite (ROSE 1000, RNASim, CRW 16S). Every nucleotide
  comparison was run by the method's own developers.
- **TWILIGHT** (Tseng et al., *Bioinformatics* 2025, doi:10.1093/bioinformatics/btaf212) is the only post-2023 paper
  that compares a fast tool with MAGUS on nucleotides, and only on RNASim (10k–1M) and AliSim.
  - At 10k sequences, MAGUS and Muscle5 are about 8% error. TWILIGHT reaches 6.5% at 1M sequences in 32 CPU-minutes.
  - On AliSim long-branch data, MAGUS is best (8.3% vs about 9%).
  - Competitors were given the *true* tree. There is no ROSE, CRW or BAliBASE.
- **FAMSA2** (Nat Biotechnol 2026, doi:10.1038/s41587-026-03095-3) and **learnMSA2** (Bioinformatics 2024,
  doi:10.1093/bioinformatics/btae381) are protein-only.
  - learnMSA2 "supports only protein alignments". On CPU with language-model embeddings it takes about 3 h per
    HomFam family (about 16× MAGUS).
  - Neither paper compares against MAGUS on nucleotides.
- **Muscle5** (Nat Commun 2022, doi:10.1038/s41467-022-34630-w) has no nucleotide benchmark beyond small Bralibase
  RNA sets and no MAGUS comparison.
- **UPP2** (2023, doi:10.1093/bioinformatics/btad007) is the broadest nucleotide comparison: UPP2 ≈ MAGUS are best on
  ROSE and large 16S, and Clustal Omega is worst. It used MUSCLE 3.8, not 5.
- **Base-method swaps exist only for PASTA, not MAGUS.**
  - PASTA with BAli-Phy subsets (Nute & Warnow 2016, doi:10.1186/s12864-016-3101-8) has the best TC everywhere at
    roughly 24 h/subset.
  - PASTA with ProbCons or G-INS-i subsets has been tested on BAliBASE only (Collins & Warnow 2018,
    doi:10.1093/bioinformatics/bty495).
  - Beating L-INS-i on ≤ 200 nucleotide sequences has only been shown for BAli-Phy with a fixed tree (Gupta et al.
    2021, doi:10.1093/bioinformatics/btab555).

## 3. Data and methods

**Data.** All from the MAGUS paper (Illinois Data Bank IDB-2643961):

| set | replicates used | sequences | notes |
|---|---|---|---|
| ROSE 1000 L1–L3, M1–M4, S1–S3 | R0 of each of the 10 conditions | 1000 | simulated DNA, ~1000 nt; random pairs are near saturation (median p-distance 0.69–0.70 on L1/L2/M2) |
| RNASim 1000 | R0 | 1000 | simulated RNA, ~1500 nt |
| RNASim 10K | R0 | 10,000 | large (TWILIGHT's home turf) |
| CRW 16S.T | R0 | 5,548 | large, biological, length-filtered reference |
| BAliBASE RV100 (8 sets) | – | 195–732 | protein, length-filtered references (`../data/balibase_clean`) |

**Tools** (4 threads, wall-clock per run, nothing else heavy running except where noted):
- MAFFT 7.505 (`--auto`; PartTree `--retree 2 --parttree`; L-INS-i, plus G-INS-i and E-INS-i on subsets).
- Clustal Omega 1.2.4; MUSCLE 5.3 (`-align`; Super5); FAMSA 2.4.1.
- TWILIGHT 0.2.3 (`code/twilight_iter.py`).
  - TWILIGHT 0.2.3's CLI needs a guide tree. I re-implemented its official Snakemake iterative mode:
    a PartTree initial tree, TWILIGHT, then FastTree on the alignment with ≥ 95%-gap columns masked, then TWILIGHT
    again; 3 iterations, the workflow default.
  - DIPPER, the workflow's default tree tool, is not installed. "twilight-1" is one iteration on the PartTree tree.
- learnMSA2 was not run: it is protein-only, and on CPU its own paper reports about 3 h per family.
- Scoring: FastSP against the reference; error = (SPFN+SPFP)/2.

**MAGUS reference points.**
- *MAGUS(pub)*: the authors' published MAGUS(Fast) alignment of the same replicate, rescored. Its runtime is the
  paper's own (different hardware).
- *MAGUS(4c)*: earlier reruns of MAGUS(Fast) with the paper's flags on this same 4-core machine type
  (`../experiments/runs/<rep>/prep.json`; 5 of the 8 BAliBASE sets; 16S.T and 1000M1 not available).
- I did **not** rerun MAGUS fresh in this session. At 15–35 min per 1000-sequence dataset (and over 1 h for 16S.T) it
  did not fit the budget next to the other runs. The 4c measurements and two dedicated idle-machine timing sessions
  (1000L3_R0: 1857 s, 1000S2_R0: 969 s, BBA0101: 798 s; branches `cs581-timing-a/b`) give MAGUS's cost on this
  hardware.

**Statistics.**
- Paired by replicate: mean Δerror in points, W/T/L with a ±0.05 tie band, Wilcoxon signed-rank p.
- With one replicate per ROSE condition, "paired by replicate" means paired across conditions.

**Subset test** (`code/subsets.py`):
- *magus40*: 2 of MAGUS's own 25 decomposition subsets per replicate (≈ 40 sequences, clade-like), from the cached
  MAGUS runs. 30 subsets over 15 replicates: 9 ROSE, RNASim, 16S.M and 4 BAliBASE.
- *bb200*: MAGUS's first 200-sequence backbone (8 random sequences per subset).
- *rand40*: 40 uniformly random sequences.
- Baseline: MAGUS's exact subset command (`mafft --localpair --maxiterate 1000 --ep 0.123`), rerun here. It
  reproduces MAGUS's cached subset alignments: 26 of 30 identical, the rest within ±0.3 points (MAFFT 7.450 vs 7.505).
- Subset jobs ran single-threaded.

**Merge pilot** (`code/merge_pilot.py`):
- For a replicate with a cached MAGUS run, realign MAGUS's 25 subsets with tool X.
- Keep MAGUS's own 10 L-INS-i backbones, run only the GCM merge (`gcmx.run_magus -s -b`, paper flags), and score.
- The control (MAGUS's own L-INS-i subsets) reproduces MAGUS(4c) exactly (e.g. BBA0067: 26.28 vs 26.28).

**Timing caveat.** To fit the budget, one or two single-threaded subset or merge jobs ran alongside the 4-thread
full-dataset runs from about 02:05 UTC onwards. Wall-clock for the fast tools is therefore an upper bound, perhaps
up to ~1.5× inflated. That matters little next to the 10–100× gaps to MAGUS.
