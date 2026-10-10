# Prior art and novelty check (2026-10-10, ~30 min of literature and source-code search)

| question | finding | confidence |
|---|---|---|
| Weighted ASTRID (support/length) for multi-copy / GDL trees? | **Not found.** The wASTRID paper (Liu & Warnow, WABI 2022 / AMB 2023, PMC10355063; bioRxiv 10.1101/2022.05.24.493312) names GDL / MUL-trees and DISCO+wASTRID as future work but does not test them. The wASTRID code (github.com/RuneBlaze/internode) README: "ASTRID-multi (see also DISCO) is still not implemented". | medium-high |
| ASTRAL-Pro3 weighting? | **None.** ASTER v1.25.3.8 `astral-pro3` has no support/length weighting flag. The ASTER paper (Zhang, Nielsen & Mirarab, MBE 2025, 42(8):msaf172) states: "no scalable algorithm has been devised to allow quartet weighting (similar to wASTRAL) for multi-copy genes, and thus, ASTRAL-Pro does not use gene tree branch support or branch length in inferring the topology". Weighted ASTRAL is the separate `wastral` binary (`--mode 1` = hybrid, the default). | high |
| DISCO + wASTRAL / DISCO + wASTRID in a benchmark? | **Not found.** DISCO paper (Willson et al., Syst Biol 2022) used ASTRAL, ASTRID and concatenation. The closest is wQFM-DISCO (PMC11634537). | medium |
| wQFM-GDL weighting? | **No.** Quartet weights are speciation-driven quartet counts with locus-aware normalisation; support and length are used only for optional post-hoc annotation (README; bioRxiv 10.1101/2025.04.04.647228 v4). | high |
| Asteroid | Multi-copy via the minimum distance over copy pairs (no tagging). Options: `--use-gene-bl` (raw branch lengths), `--min-bl X` (contract short branches), per-gene weights. No support weighting. **Closest prior art** for length-based multi-copy distances. | high |
| TREE-QMC | Weighted (`--hybrid`, `-w`, `-c` contraction) and multi-individual (`-a`) but no GDL tagging. | medium |
| wASTRID formulas | wASTRID-s: d_G(u,v) = Σ_{e on path} s(e), leaf edges s = 1, s' = max(0, (s − lb)/(ub − lb)); a fake root edge counts once. wASTRID-pl: Σ ℓ(e) divided by the tree's longest leaf-to-leaf path. No contraction. wASTRID-s with aBayes was best and matched wASTRAL-h at ~100× the speed. | high |

**Conclusion.** A GDL-aware (rooted, tagged, orthologous pairs only), support-weighted or length-weighted
internode-distance method appears to be new.
- The closest published pieces are Asteroid `--use-gene-bl`, untested DISCO+wASTRID, and wASTRAL / weighted TREE-QMC
  on single-copy trees.
- Weighted multi-copy quartet scoring is explicitly called unsolved by the ASTER authors.
