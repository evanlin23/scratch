## Data

| dataset | source | used here |
|---|---|---|
| FastMulRFS simulations (= ASTRAL-Pro "S100") | Molloy & Warnow 2020, Illinois Data Bank **doi:10.13012/B2IDB-5721322_V1** (true-trees-psize-*, estimated-gene-trees-psize-*) | 100 species; DL rate 1e-10 / 2e-10 / 5e-10 × Ne 1e7 / 5e7; reps 01–10; first 100 genes of `genes-with-gt3species.txt`; SimPhy true gene + locus trees; RAxML trees from 100 bp and 25 bp |
| DISCO simulations | Willson, Roddur & Warnow 2022, Illinois Data Bank **doi:10.13012/B2IDB-4050038_V1** (`trees.tar.gz`) | 100 species, 1000 gene families: default, gdl_{1e-10,5e-10,1e-9}_{0,05,1}, ils_1e4, ils_2e8, missing_1000; gtrees_10000_l1 (50–10,000 genes); species_1000 (1000 species). SimPhy true gene + locus trees; gene trees estimated from 100 bp |
| previous pilot | branch `claude/cs581-astridpro2`, `results/disco_runs.jsonl` | DISCO estimated-tree runs (fast methods reps 01–07, ASTRAL-Pro3 reps 01–02), same software versions; reused, see `results/prior/` |

## Methods (`code/bench.py`; all 1 thread unless stated)

- **ASTRID-Pro** (`code/apro.cpp`, from the astridpro2 branch; now with optional OpenMP over genes, `-t`), MinDup
  rooting + species-overlap tags, orthologous pairs, speciation nodes counted; FastME 2.
- **ASTRID-multi** (same code, all pairs), **ASTRID-DISCO** (DISCO v1.4.1 at GitHub HEAD a55df08 + ASTRID),
  **DISCO-ASTRAL** (DISCO + ASTRAL-IV), **ASTRAL-Pro3** (ASTER v1.25.3.8, bioconda), **Asteroid** (GitHub HEAD,
  missing-data correction on), **wQFM-GDL** v1.0.3 (GitHub abdur-rafi/wQFM-GDL, tree mode `-t`, 4 GB heap).
- **True tags** (`code/truetag.py`). A locus-tree node is a speciation iff its time equals a species-tree node time
  (SimPhy places speciations exactly there; checked: ≥ 99.8% of these nodes have disjoint child species sets). A
  gene-tree node is a duplication iff the locus-tree LCA of its leaves' loci is a duplication. The SimPhy root is
  kept. `*-tt` methods use these: ASTRID-Pro with `-T`, and DISCO's own `decompose()` applied at the true
  duplications (`code/disco_tt.py`) followed by ASTRID or ASTRAL-IV.
- **Tag accuracy** (`code/tagacc.py`): on true gene trees, every cross-species leaf pair is an ortholog iff its LCA is
  tagged S. Compared: truth vs ASTRAL-Pro3 `-T` and vs DISCO's MinDup root + overlap tags (the scheme ASTRID-Pro
  and ASTRID-DISCO use).
- **Hybrid** (Q4): `astral-pro3 -g <ASTRID-Pro tree>` with default search (`apro3-guide`) or `-r 1 -s 0`
  (`apro3-guide-fast`); control `-r 1 -s 0` without guide (`apro3-fast`). `code/qscore.py` scores the true tree,
  the ASTRAL-Pro3 tree and the ASTRID-Pro tree under ASTRAL-Pro3's objective (`-C -c`).
- Peak memory: `ru_maxrss` of each method's processes from `os.wait4`.
- Statistics: FN rate (= RF for binary trees), paired by (condition, replicate, level), two-sided Wilcoxon on non-zero
  differences, W/T/L, Holm within each pre-registered family.
