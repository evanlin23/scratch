## Methods

**Data (44 alignments, 52-445 taxa).**
* *Simulated, 24*: RAxML Grove (github.com/angtft/RAxMLGrove) empirical DNA trees with their fitted GTR(+I)+G
  parameters, stratified 8 each in 50-99 / 100-199 / 200-499 taxa, original number of sites (495-4,747); simulated with
  AliSim (IQ-TREE 3.1.4) on the RG tree (no indels). True tree = RG tree with branches <= 1e-5 contracted.
  (`code/simulate.py`, seed 581.)
* *Empirical DNA, 12*: random taxon subsamples (100/150/200/300) of the curated CRW 16S rRNA reference alignments
  16S.3, 16S.T, 16S.B.ALL (MAGUS paper data, Illinois Data Bank IDB-2643961), fragments (<1000 non-gap) removed,
  all-gap columns removed (1,770-4,738 sites).
* *Empirical protein, 8*: random 100-200-sequence subsamples of the 8 BAliBASE RV100 reference alignments (107-4,038 columns).
  (`code/prep_empirical.py`.) TreeBASE and the paper's Dryad archive were not reachable from this machine
  (Dryad is behind a bot wall; treebase.org timed out), so the paper's own 300 MSAs could not be used.

**Runs** (`code/run_jobs.sh`; all single-threaded, 4 jobs concurrently on the 4-core VM, so wall times include
memory-bandwidth contention but every method saw the same conditions). Fixed model, no ModelFinder:
GTR+F+I+G4 (DNA), LG+G4 (protein).
* `iqdef1`: IQ-TREE 3.1.4 default search (100 parsimony starts, 20 NNI-optimized candidates, stochastic perturbation,
  stop after 100 unsuccessful iterations), seed 1, run with a 12-line instrumentation patch (`code/iqtree_trace.patch`)
  that writes, for every iteration, the wall time (ms), the iteration's score, the best score, whether a new best
  *topology* was found (the event that resets IQ-TREE's counter), and that tree. Output (trees, lnL, iteration counts)
  is identical to the stock binary on the test case.
* `iqdef2`: same, seed 2 (run-to-run noise of the default itself; datasets with default wall < 700 s only).
* `iqfast`: IQ-TREE `--fast`.
* `rxfast`: RAxML-NG 2.0.3 `--fast` (= simplified SPR search from 1 parsimony tree + KH-mult stopping: the published
  rule as released).
* `rxclassic1`: RAxML-NG 2.0.3 `--search --tree pars{1} --opt-topology classic` (the v1.2-style search from one
  parsimony start; used to reproduce the paper's KH-mult speed-up).
* `iqnstop20`: real IQ-TREE run with `-nstop 20` on 6 datasets, to validate the offline replay.

**Stopping rules replayed offline on the `iqdef1` trace** (`code/replay.py`). Because the search is seeded, a run
stopped at iteration s is identical to the default run truncated at s, so each rule's tree and time can be read off
the trace: time = trace time at s + the default run's post-search time (final model/branch-length optimization).
* `nstopN` (N = 50, 20, 10): IQ-TREE's own rule with a smaller patience (the non-statistical baseline any new
  rule must beat).
* `KHpatN` (N = 100, 50, 20): a new-best topology only counts as a success if a one-sided fast KH test
  (normal approximation on per-site lnL differences, as in Togkousidis et al.; alpha = 0.05) says it is significantly
  better than the last *significant* best tree; stop N iterations after the last significant success. KHpat100 is
  the "drop-in" version with IQ-TREE's patience.
* `KHwinK` (K = 10, 20): RAxML-style: after the 20 initial candidate iterations, every K iterations test the
  current best tree against the best tree K iterations earlier; stop at the first non-significant window.
  `KHwin10mult`: Bonferroni over the number of new-best topologies in the window (the KH-mult analogue).

**Scoring.** For each dataset, all distinct trees (every new-best topology of the trace + all comparator trees,
pruned to IQ-TREE's de-duplicated taxon set) are re-evaluated in a single IQ-TREE `-z ... -n 0 -wsl -zb 1000 -au`
run, so all lnL are under one model fit (re-estimated on the first tree) with branch lengths optimized per tree,
and give the site lnLs used by the KH tests and AU p-values among all candidates. Reported: speed-up vs default
(default wall / method wall), dlnL = lnL(method) - lnL(default), W/T/L on dlnL (tie band +-0.5), Wilcoxon signed-rank
test on dlnL, fraction of AU-plausible trees (p-AU > 0.05), normalized RF to the default tree, and (simulated) FN rate
to the true tree and its paired difference to the default.
