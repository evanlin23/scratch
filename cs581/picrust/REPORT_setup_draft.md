## 1. Setup facts (PICRUSt2 2.6.3, bioconda, installed 2026-10-10)

| item | value | where |
|---|---|---|
| PICRUSt2 | 2.6.3 (`picrust2=2.6.3 pyhdfd78af_2`, bioconda); reference files from the v2.6.2 tarball, as the bioconda post-link script does | `bin/.picrust2-post-link.sh` |
| EPA-ng shipped | **0.3.8** (`epa-ng 0.3.8 h077b44d_7`, bioconda) | install log |
| default placer | `-t/--placement_tool ... default="epa-ng"` | `bin/picrust2_pipeline.py:66-67`, `picrust2/place_seqs.py` |
| alignment | `hmmalign --trim --dna --mapali bac_ref.fna ... bac_ref.hmm study.fna` | `picrust2/place_seqs.py:45` |
| EPA-ng call | `epa-ng --tree T --ref-msa R --query Q --chunk-size 5000 -T p -m bac_ref.model -w out --filter-acc-lwr 0.99 --filter-max 100` (no `--no-pre-mask`, no `--rate-scalers`, so premasking on and rate scalers `auto`) | `picrust2/place_seqs.py:178-187` |
| bacterial reference (default, "PICRUSt2-SC", GTDB r214) | **26,868 tips**, 1,578 alignment columns, model `GTR+FU+G4m` | `default_files/bacteria/bac_ref/` |
| archaeal reference | 1,002 tips (below the 2,000-tip trigger; not affected) | `default_files/archaea/arc_ref/` |
| previous default (v2.0–2.5, "oldIMG", still available) | one tree, 19,493 bacteria + 406 archaea ≈ 19.9K tips (also > 2,000) | PICRUSt2 wiki, PICRUSt2-SC page |
| query alignment start | hmmalign output is 1,588–1,611 columns; **no ASV starts at column 0**. V4 ASVs (ocean) start at column 549; mammal (V4–V5, 400 nt) at 1,005–1,006 | `results/*_placement.json` |

So the default PICRUSt2 path satisfies all three conditions of bug 2: > 2,000 tips (26,868), `--rate-scalers auto`
(EPA-ng log: "Automatic switching of use of per rate scalers"), premasking on, and queries that start mid-alignment.
