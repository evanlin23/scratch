# fragscale: do two-step pipelines BEAT RAxML-NG on fragmentary data? (CS581 pilot)

Branch `claude/cs581-fragscale`. Pre-registration: `PREREG.md` (pushed before any tree was scored;
3 amendments, listed there; amendment 2 is post hoc and only adds a secondary comparator).
Code: `code/` (`pipe.py` arms, `anytime.py` anytime recorder, `align_witch.py` estimated alignments,
`summarize.py` tables/plot). Results: `results/runs.jsonl` (one row per finished arm), `results/tables.md`
(generated), `results/anytime.tsv`, `results/anytime.png`.

@@VERDICT@@

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

@@RESULTS@@
