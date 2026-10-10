# CS581 ML project scouting: validation, pilots, recommendation

**TL;DR**

*Validation.* The harness reproduces published numbers:
- the Park et al. 2021 RNASim1000 and 1000M1-HF tables, exactly, from their deposited trees;
- their FastTree numbers from our own runs (14.9% on RNASim1000);
- their 24-hour RAxML-NG result with a 35-CPU-minute RAxML-NG v2 run (24.7% vs 24.9%);
- PASTA's own trees, to within about 1 point.

*The headroom is on fragmentary data, not full-length data.*
- On full-length ROSE/RNASim data, all ML methods are within about 1–2 FN points of each other, and even the true alignment helps only about 1 point.
- Masking is a negative result.
- On fragmentary data (1000M1-HF), FastTree gets 49% vs RAxML-NG's 25%.

*Pilot result.* A simple fragment-aware heuristic matches RAxML-NG at about 1/3 of the CPU:

| method | FN | CPU time |
|---|---|---|
| fragment-aware (full-length backbone tree, then constrained RAxML-NG) | 25.2% | 12 min |
| RAxML-NG | 24.7% | 35 min |

That is n=5 replicates, and it beats the published IQ-TREE 2 (30.2%) and GTM (28.4%) trees.

*Recommendation.* Develop that heuristic (§3.1). Alternatives: a starting-tree / time-budget study (§3.2), and GTM + RAxML-NG polish (§3.3).

Branch `claude/cs581-ml`. Everything here can be regenerated from `cs581/ml/code/`:

| file | what it does |
|---|---|
| `fetch.py` | streams the inputs from the MAGUS paper's `Results.zip` (IDB doi:10.13012/B2IDB-2643961_V1) |
| `validate_park.py` | rescores the published Park et al. 2021 trees (IDB doi:10.13012/B2IDB-7008049_V1) |
| `runtrees.py` | runs FastTree, IQ-TREE 3 and RAxML-NG, and scores each tree with `treeerr.py` (dendropy) |
| `colsupport.py`, `make_variants.py` | masking pilot |
| `frag.py` | fragment-aware pilot |
| `dnc.py` | GTM divide-and-conquer pilot |
| `make_tables.py` | regenerates `results/tables.md` (every number below) |

The literature review is in `literature.md`.

Software and setup:
- FastTree 2.1.11 (double precision), IQ-TREE 3.1.4, RAxML-NG 2.0.3.
- Every job used 1 thread, with 4 jobs at a time on a 4-core VM.
- CPU time is the user+sys time of the tool process.
- FN is the missing-branch rate vs the true tree: the fraction of true internal edges absent from the estimate.

## 1. What was validated (our numbers vs. published ones)

**Park, Zaharias & Warnow 2021 (*Algorithms* 14:148).** These are the benchmarks behind the ML lecture's "RF error / runtime of ML heuristics" slides. The authors deposited inputs *and* output trees (IDB-7008049). Rescoring their trees with our scorer reproduces their tables:

| dataset (paper table) | FastTree | IQ-TREE 2 | RAxML-NG | GTM (IQ-TREE start) |
|---|---|---|---|---|
| RNASim1000, 5 reps (Table 3): paper | 14.9 | 15.1 | 15.1 | 14.4 |
| RNASim1000: published trees rescored by us | **14.9** | **15.1** | **15.1** | **14.4** |
| RNASim1000: our own fresh runs | **14.9** (FastTree) | 15.6 (IQ-TREE 3 `--fast`) | – | – |
| 1000M1-HF, 5 reps (Table 5): paper | 50.9 | 30.2 | 24.9 | 28.4 |
| 1000M1-HF: published trees rescored by us | 51.4 | **30.2** | **24.9** | **28.4** |
| 1000M1-HF: our own fresh runs | 48.9 | 37.5 (`--fast`) | **24.7** (RAxML-NG 2.0.3, 1 start, about 35 CPU-min) | – |

Notes on this table:
- 1000M1-HF is ROSE 1000M1 with 500 of the 1000 sequences cut to fragments (about 25% of the median length).
- The input is not deposited as a single file. We rebuilt it from their per-subset alignments, which are full-width, and checked it against the ROSE originals: all 1000 sequences match, and 500 are identical.
- The published RAxML-NG trees are the `lastTree` of a run stopped at the 24 h cap (RAxML-NG 1.0.1, 20 starts, 2 threads).

**PASTA's own tree (`pasta.tre`, MAGUS paper data).** PASTA's final tree is FastTree `-gtr -gamma -fastest` on its final alignment (PASTA's default config). We re-ran that command on the published `pasta_align.txt`:

| condition | published `pasta.tre` FN (mean of R0–R2) | our FastTree `-fastest` FN |
|---|---|---|
| 1000M2 | 11.2 | 10.2 |
| 1000M3 | 7.2 | 7.2 |
| 1000L1 | 12.9 | 11.7 |
| RNASim 1K | 15.2 | 15.3 |

- Per-replicate differences are up to 2.5 points; RF between the two trees is 0.5–7.7%.
- So the setting is reproduced, but not bit-for-bit. Likely reasons: FastTree 2.1.10 vs 2.1.11, and PASTA seeding FastTree with its previous tree.
- Per-replicate table: `results/pasta_tre_check.md`.

**MAGUS-paper data, FastTree / IQ-TREE `--fast` on true / MAGUS / PASTA alignments (R0–R2):** see `results/tables.md`. The MAGUS paper reports alignment error only, so there is no published tree error to match. The relevant published claim is UPP 2015, Table 2: on these ROSE models, FastTree on PASTA alignments is about 1.3 points of FN worse than on the true alignment. We see the same small gap.

| FN, FastTree | 1000M2 | 1000M3 | 1000L1 | RNASim 1K |
|---|---|---|---|---|
| true alignment | 9.4 | 8.0 | 9.3 | 15.1 |
| MAGUS (`gcm.txt`) | 10.8 | 7.7 | 11.2 | 16.1 |
| PASTA | 10.2 | 7.2 | 12.0 | 15.0 |

Two takeaways for project choice:
1. The pipeline and data are validated against published numbers, so an apples-to-apples project is feasible.
2. On full-length ROSE/RNASim data the whole spread is small. The true alignment is only about 1 point better than MAGUS, and RAxML-NG is barely better than FastTree (§2.2). The big, unsolved gap is on **fragmentary data**: on 1000M1-HF, FastTree gets 49–51% FN, while RAxML-NG reaches about 22–25%.

## 2. Pilots

Every number below is regenerated by `code/make_tables.py` into `results/tables.md`. The "CPU min" column includes the cost of any start or backbone tree a method uses.

### 2.1 Pilot A: fragment-aware two-phase ML (the most promising)

**Idea.** On fragmentary data the hard part is placing the fragments. The full-length sequences alone give an easy, accurate tree.
1. Build a "backbone" tree on the full-length sequences only (ungapped length ≥ 50% of the median; about 520 of the 1000) with FastTree (`-gtr -gamma`, about 1 minute).
2. Run RAxML-NG (GTR+G, one parsimony start) on **all** sequences, with the backbone as a non-comprehensive topological constraint (`--tree-constraint`). Only the fragments move. This is `frag_constr`.
3. Optional: an unconstrained RAxML-NG search started from that tree (`frag_polish`).

Code: `code/frag.py`.

**Data.** 1000M1-HF, R0–R4: the exact inputs and the published trees of Park et al. 2021 (IDB-7008049). Comparisons are paired on each replicate against RAxML-NG with one parsimony start (our runs, 1 thread).

| method | n | mean FN | Δ vs RAxML-NG (wins/ties/losses) | mean CPU min | mean ΔlnL vs RAxML-NG |
|---|---|---|---|---|---|
| FastTree | 5 | 48.9% | +24.2 (0/0/5) | 1.9 | – |
| IQ-TREE 3 `--fast` | 5 | 37.5% | +12.8 (0/0/5) | 7.2 | – |
| RAxML-NG, 1 parsimony start | 5 | 24.7% | – | 35.4 | – |
| RAxML-NG from IQ-TREE `--fast` tree | 5 | 24.5% | −0.2 (3/0/2) | 49.0 | +1.8 |
| **frag_constr** (FastTree backbone + constrained RAxML-NG) | 5 | 25.2% | +0.5 (2/0/3) | **12.1** | −36.2 |
| frag_polish (+ unconstrained RAxML-NG from frag_constr) | 5 | 24.5% | −0.2 (2/0/3) | 49.6 | +2.3 |
| *published:* IQ-TREE 2 / GTM (IQ-TREE start) / RAxML-NG 24 h | 5 | 30.2% / 28.4% / 24.9% | | | |

**Reading.**
- `frag_constr` is within 0.5 FN points of a full RAxML-NG search at about **one third of the CPU time** (12 vs 35 min).
- It is far more accurate than the published IQ-TREE 2 trees (30.2%) and GTM trees (28.4%).
- Its likelihood is about 36 log-likelihood units worse, yet its trees are as accurate. This is an instance of "better score ≠ better tree" on fragmentary data, and worth reporting.
- Polishing buys the last 0.7 points but costs all of the time savings.
- Backbone quality drives the error. Replicate 4 has the worst backbone (FN 16.9% vs 11.7–12.9% elsewhere) and is also the worst final tree (31.2%). So a better backbone (IQ-TREE or RAxML-NG on the roughly 520 full-length sequences) is the obvious first lever. It was not tested here.

**Honest verdict.** We see a 3× speed-up at near-parity on n=5, with no accuracy gain yet. The method idea is simple, the data and baselines are published, and there are several untested knobs: backbone method, fragment threshold, budgeted polish, placement instead of constrained search. That makes it a good 4-week project (§3).

**Side finding (validation).** RAxML-NG 2.0.3 with *one* parsimony start and 1 thread gave 24.7% in about 35 CPU-minutes. The published runs (RAxML-NG 1.0.1, 20 starts, 2 threads) were stopped at the 24 h cap and gave 24.9%. Their per-replicate values, 22.9/24.0/23.5/24.4/29.8 vs ours 21.8/23.3/25.7/23.2/29.8, match closely. So the expensive baseline can be reproduced cheaply.

### 2.2 Pilot B: starting trees for RAxML-NG on full-length data

**Idea.** Seed RAxML-NG with the FastTree tree (`--tree ft.nwk`) instead of a parsimony tree.

**Data.** MAGUS alignment, replicate R0 of each condition. This is n=1 per condition, so read it as indicative only.

| condition (R0, MAGUS aln) | FastTree | IQ-TREE `--fast` | RAxML-NG (pars start) | RAxML-NG (FastTree start) | CPU min, pars → FT start (incl. FastTree) | ΔlnL (FT − pars) |
|---|---|---|---|---|---|---|
| 1000M2 | 10.7% | 9.7% | 9.9% | 9.5% | 45 → 49 | −0.3 |
| 1000M3 | 7.7% | 7.7% | 8.0% | 7.2% | 35 → 26 | −0.7 |
| 1000L1 | 11.7% | 11.7% | 12.2% | 11.3% | 39 → 40 | +2.2 |
| RNASim 1K | 14.4% | 15.3% | 13.8% | 13.9% | 63 → 51 | −1.1 |

**Reading.**
- The FastTree start was more accurate on 3 of 4 conditions (by 0.4–0.9 points), tied on RNASim, and was never slower by more than about 10%.
- The likelihood scores differ by only about 1–2 units. Either start lands on essentially equally good optima, with different topologies.
- On full-length data the whole field is compressed. RAxML-NG improves on FastTree by at most 0.8 points here, and is *worse* than FastTree on 1000M3 and 1000L1. This agrees with Liu, Linder & Warnow 2011.

**Verdict.** Cheap and safe, but the effect is small. It makes a good secondary experiment inside the main project, or alternative 1.

### 2.3 Pilot C: alignment-confidence column masking (negative result)

**Idea.**
- Score every column of the MAGUS alignment by the fraction of its homology pairs that the independent PASTA alignment also asserts (`code/colsupport.py`; linear time).
- Drop the columns scoring below 0.5 (`mask50`) or 0.7 (`mask70`).
- Controls:
  - "oracle" masks that score columns against the *true* alignment;
  - a soft ensemble that concatenates the MAGUS and PASTA alignments column-wise (`gcm_pasta`), so each homology is counted once per alignment asserting it.

The support score does predict column correctness: correlation with true precision is 0.58–0.74.

FastTree FN change vs the unmasked MAGUS alignment (R0–R2; negative = better):

| alignment | 1000M2 | 1000M3 | 1000L1 | RNASim | all 12 (wins/ties/losses) |
|---|---|---|---|---|---|
| mask50 | +0.3 | +0.1 | −0.3 | +0.5 | +0.2 (4/0/8) |
| mask70 | +0.8 | +0.1 | +0.2 | +0.4 | +0.4 (3/1/8) |
| oracle50 | −0.0 | +0.4 | +0.1 | −0.0 | +0.1 (4/1/7) |
| oracle70 | −0.5 | +0.3 | −0.1 | +1.4 | +0.3 (5/0/7) |
| MAGUS+PASTA concatenated | −0.9 | +0.3 | +0.3 | −0.5 | −0.2 (7/1/4) |
| PASTA alignment instead | −0.6 | −0.5 | +0.8 | −1.1 | −0.3 (9/0/3) |
| **true alignment instead** | −1.4 | +0.3 | −1.8 | −1.0 | **−1.0 (9/0/3)** |

**Verdict.** Masking does not help, even with perfect knowledge of which columns are wrong. That matches Tan et al. 2015 (Syst Biol 64:778). The ceiling is low anyway: switching all the way to the *true* alignment gains only 1.0 point. **Not recommended** as a project on these data.

### 2.4 Pilot D: GTM divide-and-conquer (+ polish)

DNC_RESULTS

## 3. Recommendation

### 3.1 Recommended: fragment-aware ("backbone-constrained") ML for datasets with fragmentary sequences

**Problem.** Many real alignments mix full-length and fragmentary sequences: amplicons, reads, partial genes. ML heuristics cope badly with this:
- On 1000M1-HF, FastTree gets 49–51% FN and IQ-TREE 30–38%.
- GTM, the leading divide-and-conquer pipeline, gets 28%.
- Only a full RAxML-NG search gets about 25%.

So this is the setting where "current heuristics are either slow or inaccurate" is most visible (ML lecture take-home), and where the gap between methods is about 25 points rather than 1.

**Method (new).** A two-phase heuristic:
1. Choose a backbone of full-length sequences (length threshold τ).
2. Estimate a backbone tree with a fast-but-good method.
3. Run constrained ML (RAxML-NG `--tree-constraint`) so that only the fragments are placed, jointly and under the full ML criterion.
4. Optionally, run an unconstrained ML polish with a budget.

The pilot shows parity with RAxML-NG at about 1/3 of the CPU. The project tries to turn parity into a win by varying:
- (a) the backbone estimator: FastTree / IQ-TREE `--fast` / IQ-TREE / RAxML-NG on the roughly 50% full-length sequences;
- (b) τ;
- (c) the polish budget: RAxML-NG `--spr-radius` or a fixed number of rounds, vs. full;
- (d) placement-then-polish: EPA-ng placement of each fragment onto the backbone, grafted with gappa, then a short unconstrained polish. This is cheaper, but fragments are placed independently.

**Is it new?**
- Smirnov & Warnow 2021 compared "align then ML" with "backbone tree + independent placement (pplacer)". Placement was clearly worse.
- Park et al. 2021 compared GTM, TreeMerge and Constrained-INC.
- We (literature review, `literature.md` §5) found no paper that uses the full-length tree as a *constraint for a joint ML placement of the fragments*, or as the seed for a polish, on these benchmarks.
- It is a small, well-defined variation. That is the right size for 4 weeks.

**Datasets (all public, with baselines):**
1. **1000M1-HF** and RNASim1000 (Park, Zaharias, Warnow 2021, *Algorithms* 14:148, doi:10.3390/a14050148). Inputs and published trees are at Illinois Data Bank doi:10.13012/B2IDB-7008049_V1. Our scripts rebuild the exact input alignments and rescore the published trees to the paper's Table 5 numbers (`code/validate_park.py`).
2. **Smirnov & Warnow 2021** (*Syst Biol* 70:268, doi:10.1093/sysbio/syaa058) low/high-fragmentation versions of ROSE 1000M1–M4 and RNASim 1K, on Dryad, doi:10.5061/dryad.8pk0p2nj8 (DOI as given in the paper; it did not resolve in DataCite from this VM, so verify it before relying on it).
   - Published baselines (FN %, high fragmentation; M1/M2/M3/M4/RNASim): UPP-RAxML 37.0/30.4/23.7/16.7/37.7; UPP-pplacer 48.8/43.7/38.0/32.0/50.7.
   - These are on *estimated* (UPP) alignments. UPP is in bioconda (`sepp`) if we want to run on estimated alignments too.
   - We could not open the Dryad landing page from the API. Fallback: re-create the fragmentation on the ROSE originals in the MAGUS data bank entry (doi:10.13012/B2IDB-2643961_V1) with the paper's protocol (25%/50% of the sequences cut to about 50%/25% of the median length).

**Criteria.**
- FN (missing-branch) rate vs the true tree. This is the primary criterion, matching Park et al. and Smirnov & Warnow.
- RF.
- The GTR+G log-likelihood of the final tree, re-evaluated with `raxml-ng --evaluate` on the same alignment.
- CPU time and peak memory.
- Paired per-replicate comparisons (wins/ties/losses and a Wilcoxon signed-rank test).

**Baselines.**
- FastTree 2 and IQ-TREE 3 (`--fast` and default).
- RAxML-NG: one parsimony start, and the v2 default `--tree auto`.
- Published GTM / IQ-TREE 2 / RAxML-NG trees (IDB-7008049).
- Published UPP-RAxML / UPP-pplacer numbers (Smirnov & Warnow).

**Plan.**

| week | work | done when |
|---|---|---|
| 1 | Set up: data, scripts (all in `cs581/ml/code/`). Reproduce Park Table 5 from published trees (done) and from our own runs. Obtain or re-create the Smirnov–Warnow fragmentary conditions. Run baselines, single-threaded, on 1000M1-HF R0–R4 plus 2 more conditions. | baseline table matches the published numbers within noise |
| 2 | Implement the variants (a)–(d) above as flags of `frag.py`. Install EPA-ng and gappa from bioconda for (d). Tune on R0–R1 of 1000M1-HF only. | pick ≤ 3 variants |
| 3 | Main experiment on the held-out replicates (R2–R9) of 1000M1-HF and 2–4 Smirnov–Warnow conditions (high and low fragmentation). Budget: about 35 CPU-min per RAxML-NG-class run, so 8 reps × 4 conditions × about 6 methods ≈ 100 CPU-h. Feasible on a 4-core laptop over a week, or faster on the campus cluster. | all runs done |
| 4 | Analysis: accuracy vs CPU curves, lnL vs FN, and backbone error vs final error. Optionally run on UPP alignments for the estimated-alignment setting. Write-up. | report and slides |

**Risks and mitigations.**
- *Only a speed-up, no accuracy gain.* This is already the pilot's outcome, and a 3× speed-up at parity with RAxML-NG, beating IQ-TREE and GTM by 3–5 points, is a publishable-style result for a course project.
- *Backbone errors propagate* (replicate 4: backbone 16.9% → final 31%). Variant (a) tests a better backbone, and the unconstrained polish (c) can repair them.
- *The constraint lowers lnL* (−36 units) even when FN does not drop. Report both, as in the pilot. "Better score ≠ better tree" is itself a finding.
- *Compute.* Use single-threaded jobs, run in parallel, and keep 1000-taxon datasets only. RAxML-NG v2 with 1 start already reproduces the published 24 h runs in about 35 CPU-min.
- *Estimated alignments.* Start with true alignments (as Park et al. did). Add UPP alignments only if time allows.

### 3.2 Alternative 1: starting trees and time budgets for RAxML-NG / IQ-TREE on the Warnow-lab benchmarks

Measure FN, lnL and CPU for RAxML-NG and IQ-TREE started from FastTree, IQ-TREE `--fast`, parsimony and random trees. Do this at several search budgets on ROSE 1000M1–M4/L1, RNASim 1K (MAGUS IDB) and 1000M1-HF (IDB-7008049).

- **Pilot (§2.2):** a FastTree start was more accurate than parsimony on 3 of 4 conditions (by 0.4–0.9 points) at similar CPU, on n=1.
- **Novelty:** low to medium. Starting-tree effects are known on empirical data (Morel et al. 2021, MBE 38:1777; adaptive RAxML-NG, Togkousidis et al. 2023, MBE 40:msad227), but not with true trees on these benchmarks.
- **Risk:** very low. This is mostly a comparison study, which the course allows.
- It can also be folded into the main project as a side experiment.

### 3.3 Alternative 2: GTM pipeline with a short RAxML-NG polish ("blending" after a disjoint tree merger)

This is the instructor's divide-and-conquer open problem: "a better DTM that allows blending". The idea:
1. GTM-merge RAxML-NG subset trees.
2. Use the merged tree as the start of a *budgeted* unconstrained RAxML-NG search, so SPR moves can repair the merge.

Data: IDB-7008049 (RNASim 1K/10K/50K, Cox1-HET, 1000M1-HF, with published GTM and RAxML-NG trees). The published gaps are GTM 28.4% vs RAxML-NG 24.9% on 1000M1-HF, and 18.7% vs 18.2% on Cox1-HET.

- The pipeline is implemented (`code/dnc.py`; GTM at github.com/vlasmirnov/GTM). Pilot D (§2.4) is a first data point.
- **Risk:** medium. The interesting regime is 10K–50K taxa, where RAxML-NG takes many hours per run. A student would need the campus cluster.

### 3.4 Not recommended (from the pilots and the literature)

- **Alignment masking/trimming before ML:** negative in our pilot, even with oracle masks, and in Tan et al. 2015. The ceiling is ≤ 1 point.
- **Deep-learning tree estimation:** published models are evaluated on ≤ about 200 taxa. Too risky for 4 weeks at 1000 taxa.

## 4. The two slide questions passed on by the orchestrating session

**"We need better heuristics" (ML lecture take-home).** Our pilots put numbers on the trade-off the slide describes:
- On full-length 1000-taxon data, the slow methods barely help. RAxML-NG improves on FastTree by ≤ 0.8 FN points, and is sometimes worse, at about 20–40× the CPU.
- On fragmentary data they help enormously: 49% → 25%.
- A heuristic that exploits the structure of the input (a full-length backbone, then constrained placement) recovers RAxML-NG accuracy at about 1/3 of the CPU.
- **Verdict:** "better heuristics" is most tractable and most valuable on heterogeneous inputs (fragments, gappy sequences), not on clean full-length data. That is what §3.1 proposes.

**"How can we scale concatenation to large datasets?" (phylogenomics part 2, slide 35)** — assessment only, no pilot.

*The bottleneck.* Concatenation (CA-ML) = ML on a supermatrix of many genes, so it inherits every ML scaling problem plus huge alignment width. NJMerge (Molloy & Warnow 2019, *Algorithms Mol Biol* 14:14, doi:10.1186/s13015-019-0151-x) reports two things on 1000 species × 1000 genes:
- RAxML could not run within 64 GB;
- the last RAxML checkpoint came only after more than 2,250 minutes on the datasets where it did run.

NJMerge + RAxML on disjoint subsets reduced the running time without sacrificing accuracy. GTM (Smirnov & Warnow 2020, *BMC Genomics* 21(Suppl 2):235, doi:10.1186/s12864-020-6605-1) is the faster DTM for the same pipeline.

*What a pilot would be.* GTM + RAxML-NG on the NJMerge datasets (Illinois Data Bank doi:10.13012/B2IDB-1424746_V1 and doi:10.13012/B2IDB-0569467_V2), vs RAxML-NG on the full supermatrix. `code/dnc.py` already implements decompose → subset RAxML-NG → GTM, and would only need a supermatrix as input.

*Why we did not run it.* A single 1000 × 1000-gene supermatrix RAxML-NG baseline would exceed this VM's 15 GB of memory and our CPU budget.

*Verdict.*
- Tractable as a course project only if it reuses the published NJMerge baseline numbers rather than re-running full-supermatrix RAxML.
- The genuinely open part is accuracy. CA-ML is not statistically consistent under ILS (lecture). A scalable CA-ML is useful mainly as a fast, strong baseline.
- This makes it a weaker choice than §3.1 for a 4-week ML project.
