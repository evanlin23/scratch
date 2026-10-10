# New project ideas beyond the slides (literature scout, 2026-10-10)

A scan of 2023-2026 papers for open problems with a cheap CPU pilot, excluding everything already
explored (see `ranking_final.md`). Each novelty check was 1-3 searches, not exhaustive. Pilots for the
top five run overnight (`cs581/code/fanout/SESSIONS.md`).

| rank | idea | source open problem | cheap pilot / go-kill test | course fit |
|---|---|---|---|---|
| 1 | Why EPA-ng gets worse on > 2,000-leaf placement subtrees | Wedell, Shen & Warnow, BSCAMPP, IEEE/ACM TCBB 2025 (doi:10.1109/TCBBIO.2025.3562281): delta error more than doubles above 2,000 leaves even with true query alignments; "EPA-ng may have some numeric issues ... Further research is needed" | RNASim backbone 5-10K leaves: delta error vs subtree size for EPA-ng default / `--no-heur` / `--baseball-heur` / pplacer / APPLES-2; does the jump disappear? | ML placement (SEPP/pplacer/EPA-ng) |
| 2 | WITCH-lite for TIPP3 / TIPP-SD | TIPP3 (PLOS CB 2025): the WITCH step is "the biggest contribution to runtime"; faster-but-accurate replacement "a promising direction"; repeated in BSCAMPP 2025 and TIPP-SD 2026 | hierarchical / BLAST-guided / adaptive-k HMM selection vs WITCH: within ~1 SPFN point at < 20% runtime? | UPP/WITCH HMM ensembles |
| 3 | Predicted-3Di (ProstT5) alignment with FoldMason; as MAGUS evidence | FoldMason (Science 2026) gains need real structures; Unicore (2025) runs ProstT5 → FoldMason with no MSA benchmark | BAliBASE SP/TC vs L-INS-i, Clustal, MAGUS; then 3Di backbones in GCM | MSA, protein benchmarks |
| 4 | KH-test early stopping in IQ-TREE 3 | Togkousidis, Stamatakis & Gascuel, Syst Biol 2025 (doi:10.1093/sysbio/syaf043): 3.9-5x faster RAxML-NG; "can seamlessly be integrated into other ... tools" | replay the rule on IQ-TREE runs: >= 2x faster at no lnL/RF loss? | ML tree search heuristics |
| 5 | Identifiability / choosing k in distance-mixture deconvolution | Arasti et al., RECOMB 2026 (PMC12871782): "characterizing unidentifiable cases" and better criteria for k left open | exhaustive small-tree enumeration; BIC/CV for k | additivity, four-point condition |
| 6 | CMAPLE inside divide-and-conquer for mixed-divergence data | CMAPLE (MBE 2024): "only works on low divergence"; plan to combine with classical methods | cluster to low-divergence pieces + graft | DCM / divide and conquer (overlaps GTM work) |
| 7 | Read-technology-aware TIPP-SD thresholds | TIPP-SD 2026 | incremental | profiling |
| 8 | WASTER for deep phylogenies / branch lengths under ILS | Zhang & Nielsen 2025 | heavy simulation | ILS |
| 9 | Ancestral reconstruction accuracy with modern aligners | Vialle et al. MBE 2018 | AliSim → aligners → IQ-TREE `--ancestral` | low novelty |

Not verified: TIPP3 reference-package size, whether FoldMason accepts ProstT5-only databases, hidden
IQ-TREE stopping rules, the coverage of the Santus et al. 2025 nf-core MSA benchmark.
