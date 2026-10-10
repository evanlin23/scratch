# EPA-ng large-tree bug: upstream audit (2026-10-10)

Source-code and issue audit of EPA-ng master (b24ea4a, 2026-06-14); nothing was run. Finding under test on
branch `claude/cs581-epang` (`cs581/epang/code/epa-ng-fix.patch`).

**Open problem it would answer.** Wedell, Shen & Warnow, BSCAMPP (IEEE/ACM TCBB 2025; bioRxiv
10.1101/2022.10.26.513936): delta error jumps "from backbone tree size of 2000 to 5000"; "EPA-ng may have
some numeric issues when placing into very large trees ... Further research is needed".

| | pattern 1: `attributes = PLL_ATTRIB_RATE_SCALERS;` (file_io.cpp make_partition) | pattern 2: `scale_buffer[i] += offset` (pll_util.cpp shift_partition_focus) |
|---|---|---|
| present in master / v0.3.8 / v0.3.9-pre | yes | yes |
| introduced | 092521d, 2018-12-11, v0.3.3 ("automatic use of per rate scalers") | line from d9b8277 (2016); wrong since 092521d made per-rate scalers automatic |
| effect | pattern-tip and site-repeat flags are OR-ed in afterwards and survive; only the SIMD arch bits are lost → trees > 2000 tips run on **scalar kernels (speed)** | per-rate scalers hold `rate_cats` entries per site (libpll-2 indexes `scaler[n*rate_cats+i]`; EPA-ng sizes them `sites_alloc * rate_cats`), so the shifted scalers are **misaligned (accuracy)** |
| reached when | `--rate-scalers auto` (default) and `large_tree()` = `tip_nodes > 2000` | the same, plus `call_focused` (Tiny_Tree.cpp): thorough branch-length phase with premasking (default); the focused logl is the placement logl |

- **Prior reports:** none in EPA-ng's 54 issues / 6 PRs; nothing relevant in pll-modules, libpll-2, raxml-ng.
  The v0.3.0 note "pll scalers were copied incorrectly" is a different, older bug.
- **Versions:** bioconda epa-ng 0.3.8 (~86K downloads). BSCAMPP reports v0.3.8 and bundles a 0.3.8 binary,
  called without `--no-pre-mask` / `--rate-scalers`, i.e. on the affected path. TIPP3-fast uses it via BSCAMPP.
- **Reach:** ~690 citations (OpenAlex, Syst Biol 2019). PICRUSt2's default placer is EPA-ng with a ~26,868-tip
  reference tree (affected path). q2-fragment-insertion uses SEPP (not affected).
- **Threshold matches:** BSCAMPP's default subtree of exactly 2,000 tips does not trigger it; 5,000 does.
- **Still to show:** that pattern 2 causes the jump (controls: `--no-pre-mask`, `--rate-scalers off`, fix),
  and what the `|=` fix buys in speed.
