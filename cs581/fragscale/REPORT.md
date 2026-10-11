# fragscale: do two-step pipelines BEAT RAxML-NG on fragmentary data? (CS581 pilot)

**STOPPED EARLY, INCOMPLETE.** The orchestrator asked to stop on 2026-10-11 ~03:50 UTC because the student
decided not to pursue this direction. All jobs were killed. Only the rows in `results/runs.jsonl` finished.
No test reached its planned n, so everything below is descriptive.

Branch `claude/cs581-fragscale`. Pre-registration: `PREREG.md` (pushed before any tree was scored;
3 amendments, listed there; amendment 2 is post hoc and only adds a secondary comparator).
Code: `code/` (`pipe.py` arms, `anytime.py` anytime recorder, `align_witch.py` estimated alignments,
`summarize.py` tables/plot). Results: `results/runs.jsonl` (one row per finished arm), `results/tables.md`
(generated), `results/anytime.tsv`, `results/anytime.png`.

## What finished (1000M1-HF, true alignment unless stated)
| arm | n | mean FN | mean CPU min |
|---|---|---|---|
| RAxML-NG full (1 parsimony start) | 6 (R0–R5) | 24.55% | 39.0 |
| **pipeline `place_ft_0.5_fix_rxfast`** | 4 (R0–R3) | 24.21% | 13.8 (0.37× RAxML-NG) |
| `constr_ft_0.5` | 3 (R0–R2) | 23.33% | 13.0 |
| RAxML-NG fast mode (post hoc) | 4 (R0–R3) | 28.42% | 14.2 |
| IQ-TREE 3 `--fast` | 3 | 37.32% | 7.1 |
| FastTree | 4 | 48.45% | 1.9 |

* **Primary (n = 4 of the planned 6):** pipeline vs RAxML-NG@T_pipe, the tree RAxML-NG would report if killed at the pipeline's CPU:
  24.2% vs 48.3%, 4/0/0 wins, Wilcoxon p = 0.125 (the smallest possible p at n = 4).
  This comparator is weak. By T_pipe, RAxML-NG often still has only its parsimony tree, because its first SPR round alone takes
  11–15 CPU-min. At 1/3× and 1/2× of its own full CPU, RAxML-NG@T is 39.4% FN (n = 6, 0/0/6 vs full, p = 0.031).
* **Fairer comparator (post hoc), RAxML-NG fast mode at about the same CPU** (14.2 vs 13.8 min): pipeline 22.7/22.9/27.1/24.2 vs
  fast mode 24.8/27.4/33.6/27.8. That is 4/0/0, with mean ΔFN −4.2 points. n = 4, so no p is meaningful.
* Pipeline vs full RAxML-NG, same 4 reps: 24.2% vs 23.5%, at 0.37× the CPU. That repeats fragml2's "tie at about 1/3 time".
* Not run or not finished: IQ-TREE anytime curves; RAxML-NG on 10K (only 3 checkpoints by 1,916 CPU-s); stock-EPA-ng arms on 10K;
  all of Q3. Only one WITCH alignment finished (M1HF R0, 17.8 CPU-min, MAFFT backbone); no tree was scored on it.
* RNASim10KHF R0: FastTree 63.3% (13.5 CPU-min); patched graft 35.9% (11.5 CPU-min, 7.9 GB peak); graft + FastTree polish 30.6% (21.4 CPU-min).

**Verdict (provisional, low n): unclear to promising.** At equal CPU the pipeline clearly beats both RAxML-NG
killed at that time and RAxML-NG's own fast mode, in 4 of 4 replicates. But with full time it only ties RAxML-NG.
"Faster AND at least as accurate" holds, as in fragml2. "Beats the best method at equal time" is shown only against
truncated or fast RAxML-NG, with n = 4 and no significance. Neither the scale question nor the estimated-alignment question was tested.

Risks noted: RAxML-NG's anytime curve is very coarse (checkpoint granularity), so the "equal time" framing depends on
the comparator. UPP/SEPP 4.5.6 is broken on bioconda. MAGUS backbones are slow on a 4-core VM. CPU-share
contention (Linux autogroups) inflated wall times.

## 1. Data
| dataset | source | what was done |
|---|---|---|
| **M1HF** (1000M1-HF), R0–R4 | exact Park, Zaharias & Warnow 2021 inputs, Illinois Data Bank IDB-7008049 (`1000M1_HF_Analysis.tar.gz`) | rebuilt from their per-subset alignments (as in fragml2); true tree = ROSE `rose.tt` |
| M1HF R5 (R6–R7 if time) | ROSE 1000M1 R5 (MAGUS Datasets.zip, IDB-2643961) | `make_frag.py`: 500/1000 sequences cut to one fragment, length ~ N(0.25 × median, 60), uniform start (protocol measured on the published 1000M1-HF) |
| **RNASim10KHF** R0 | RNASim 10K R0 (MAGUS Datasets.zip) | same fragmentation (5,000 fragments); true alignment restricted to columns with < 95% gaps among full-length sequences (8,700 → 1,626 columns), as in fragml2/epang. The Smirnov & Warnow 2021 Dryad deposit is not reachable from this VM (403), so their exact HF files were not used |
| RNASimHF R0 | RNASim 1K R0 | same fragmentation; used only in Q3 if time allowed |

## 2. Methods and versions
All tools single-threaded; CPU = user + sys of every step of an arm (shared steps counted in full for each arm);
peak RSS = max over steps. 4-core VM, several jobs at once (CPU time is robust to this, wall is not).

* RAxML-NG 2.0.3 (bioconda), GTR+G, one parsimony start, seed 1. IQ-TREE 3.1.4, GTR+G, seed 1. FastTree 2.1 (`-nt -gtr -gamma`).
* EPA-ng 0.3.8 stock (bioconda) and patched (`code/epa-ng-fix.patch`, from `claude/cs581-epang`); gappa `examine graft --fully-resolve`.
* WITCH 1.0.10 (witch-msa, pip) with a MAFFT 7.505 `--auto` backbone (amendment 3); UPP/SEPP 4.5.6 crashes in its backbone step (amendment 1).
* **Primary pipeline** `place_ft_0.5_fix_rxfast`: backbone = sequences ≥ 0.5 × median length → FastTree → RAxML-NG `--evaluate`
  (branch lengths, GTR+G) → patched EPA-ng places the rest → graft best placement → RAxML-NG fast-mode polish
  (`--opt-topology simplified --stop-rule kh-mult`) from the grafted tree on all sequences.
* `constr_ft_0.5`: RAxML-NG (1 parsimony start) on all sequences with the FastTree backbone as `--tree-constraint`.
* **Anytime curves** (`code/anytime.py`): RAxML-NG rewrites `.raxml.lastTree.TMP` (current best tree) at each checkpoint;
  a poller copies it with the process CPU time. "RAxML-NG@T" = the tree a run killed at CPU T would report (topology).
  IQ-TREE likewise from its checkpoint (`-cptime 10`), capped at the replicate's RAxML-NG CPU.
* `base_rxfastmode` (post hoc, amendment 2): RAxML-NG fast mode from one parsimony start.
* FN = missing-branch rate vs the true tree (dendropy), as in fragml2.

Full generated tables: `results/tables.md`.
