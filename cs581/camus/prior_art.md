# CAMUS pilot: prior art (checked 2026-10-09, ~30 min)

| work | ref | relation to the ideas |
|---|---|---|
| CAMUS (Willson & Warnow, ISMB 2026) | Bioinformatics 42(S1) btag245, doi:10.1093/bioinformatics/btag245; preprint bioRxiv doi:10.64898/2026.02.01.703143; data doi:10.13012/B2IDB-6892704_V1 / _V2 | Fixed ASTRAL tree. t = 0.5 was picked from {tested values} on 26-taxon data with 20 replicates. Future work proposes "CAMUS inside a heuristic search that explores different constraint trees" and a stopping rule for k. |
| CAMUS code | github.com/jsdoublel/camus (v1.0.1, 2026-03; HEAD fb6a287 used here) | Filter: sort counts f1≥f2≥f3 and keep f2 iff `floor(t·(f2+f3)) < f2−f3` (with `-q 2`, f3 is always dropped). This is a ratio rule, not a test, and does not use sample size. |
| NetCS (Dai & Molloy, WABI 2026) | LIPIcs WABI 2026 vol. 390 paper 1, doi:10.4230/LIPIcs.WABI.2026.1; code in molloy-lab/TREE-QMC | Level-1 blob reconstruction from quartets by majority vote and merge sort. Near-perfect given the true tree of blobs; end-to-end error comes from tree-of-blobs estimation. It does not benchmark CAMUS. **A direct competitor at 200 taxa in minutes.** |
| SNaQ / SNaQ.jl | Solís-Lemus & Ané 2016, doi:10.1371/journal.pgen.1005896; SNaQ.jl v1, bioRxiv doi:10.1101/2025.11.17.688917 | Pseudolikelihood. Search over the tree is free. Most accurate at 16 taxa in the CAMUS paper. |
| PhyloNet-MPL | Yu & Nakhleh 2015, BMC Genomics 16(S10):S10, doi:10.1186/1471-2164-16-S10-S10 (DOI from memory, not verified) | The fixed-tree (FT) and free-tree versions differ, and the free-tree version is better. That is direct evidence that the base tree matters for MPL. |
| NANUQ | Allman, Baños, Rhodes 2019, AMB doi:10.1186/s13015-019-0159-2 | **Per-quartet hypothesis tests** (α tree test, β star test) to classify 4-taxon sets. This is prior art for idea (b). |
| TINNiK / NANUQ+ (MSCquartets) | Allman, Baños, Mitchell, Rhodes 2024, AMB 19, doi:10.1186/s13015-024-00266-2 | Tree of blobs from per-quartet tests, then a cycle per blob. The vignette recommends varying α and β. |
| Identifiability of level-1 networks from quartet CFs | Allman, Baños, Garrote-Lopez, Rhodes 2024, Bull Math Biol 86:110, doi:10.1007/s11538-024-01339-4 | Theory. |
| Squirrel | Holtgrefe et al. 2025, MBE 42(4) msaf067, doi:10.1093/molbev/msaf067 | Level-1 network from quarnets. Also "few quarnets", arXiv:2409.06034. |
| PhyNEST | Kong, Swofford, Kubatko, Syst Biol doi:10.1093/sysbio/syae054 | Sequence-based composite likelihood. |
| InPhyNet | Kolbow, Kong, Solís-Lemus, bioRxiv doi:10.1101/2025.05.05.652278 | Divide-and-conquer network merging. Relevant to scalability. |

**Novelty check:**
- **(a) Base-tree search.** No published work does this for CAMUS, but it is announced as CAMUS future work by the same group.
- **(b) Sample-size-aware threshold.** This is new for CAMUS, but the idea is anticipated by NANUQ/TINNiK per-quartet tests, which must be cited.
- **(c) Scalability.** This is the stated future work. NetCS already claims 200 taxa in minutes for the level-1 blob problem.
