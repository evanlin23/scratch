# Does the EPA-ng rate-scaler bug change PICRUSt2's outputs? (CS581 pilot)

## Verdict (short)

**PICRUSt2 is on the affected path, and the bug changes its placements a lot. It barely changes the
predictions people use, and fixing it does not measurably improve accuracy where ground truth exists.**

- **Same path.** PICRUSt2 2.6.3 (current bioconda) ships EPA-ng 0.3.8 as its default placer. It uses a
  26,868-tip bacterial reference tree and EPA-ng's defaults: premasking on, `--rate-scalers auto`.
  hmmalign puts every ASV mid-alignment (first column 372–1,006, never 0). All conditions of bug 2 hold.
- **Placements.** With the patch, 7–43% of real ASVs move to a different edge. The share grows as reads
  get shorter: V4 7%, 400 nt 29%, 130 nt 43%. Stock likelihoods are off by up to ~200,000 log units.
  Stock makes 1.3–2.4× as many "certain" placements (LWR ≥ 0.99) as the patched build.
- **NSTI and per-ASV predictions.** NSTI changes for 31–75% of ASVs, mostly upward: stock places reads
  falsely close to a reference. Per-ASV KO vectors change for 8–39% of ASVs.
- **Sample-level outputs** (what PICRUSt2 users analyse) barely change. Stock and patched predictions
  correlate at Spearman ≥ 0.996 per sample for V4 data, ≥ 0.97 for 400 nt, and median 0.956 (min 0.87)
  for 130-nt reads. KO relative abundance moves 0.2–3.5% (median ½·L1).
- **Accuracy vs paired shotgun metagenomes** (the PICRUSt2 paper's validation data):
  - No significant change on ocean (n = 6) or mammal (n = 11).
  - On indian (130 nt, 108 samples), KO Spearman improves by +0.006 (0.7515 → 0.7574; 97/108 samples
    better, p = 3e-15), but pathway Spearman worsens by −0.003 (p = 4e-4).
- **Accuracy vs known genome content** (leave-out: 300 GTDB genomes removed from the tree, their own 16S
  amplicon placed back):
  - Placement error and KO F1 are statistically identical for stock and patched.
  - 130-nt V3 reads: 15% of placements move. Mean error is 10.59 vs 10.56 edges (p = 0.76); KO F1 is
    0.866 vs 0.866.
  - The fix mainly removes false confidence (LWR ≥ 0.99: 44 → 24) and raises NSTI.

So the bug is real on PICRUSt2's default path and visible in its intermediate outputs (placements, LWR,
NSTI, the closest reference genome, per-ASV predictions). For gene-family and pathway abundances the
effect is small, and its accuracy impact is within noise, with one small, mixed-sign exception. Reasons:

- Hidden-state prediction averages over the neighbourhood.
- Moved ASVs mostly land a few edges away, among genomes with similar gene content.
- Sample profiles are dominated by gene families common to all bacteria.

**Be honest in the project:** PICRUSt2 is a good *reach* case study ("the bug affects a tool with
thousands of citations, here is how much"), not evidence that PICRUSt2 results are wrong.

## 1. Setup facts (PICRUSt2 2.6.3, bioconda, installed 2026-10-10)

| item | value | where |
|---|---|---|
| PICRUSt2 | 2.6.3 (`picrust2 2.6.3 pyhdfd78af_2`, bioconda). Reference files come from the v2.6.2 tarball, as the bioconda post-link script does (the post-link download was blocked here; the same tarball was taken from GitHub tag v2.6.2) | `bin/.picrust2-post-link.sh` |
| EPA-ng shipped | **0.3.8** (`epa-ng 0.3.8 h077b44d_7`, bioconda) | install log |
| default placer | `-t/--placement_tool ... default="epa-ng"` (SEPP is the alternative) | `bin/picrust2_pipeline.py:66-67` |
| alignment | `hmmalign --trim --dna --mapali bac_ref.fna ... bac_ref.hmm study.fna` | `picrust2/place_seqs.py:45` |
| EPA-ng call | `epa-ng --tree T --ref-msa R --query Q --chunk-size 5000 -T p -m bac_ref.model -w out --filter-acc-lwr 0.99 --filter-max 100`. There is no `--no-pre-mask` and no `--rate-scalers`, so premasking is on and rate scalers are `auto` (EPA-ng log: "Automatic switching of use of per rate scalers") | `picrust2/place_seqs.py:178-187` |
| bacterial reference (default "PICRUSt2-SC", GTDB r214) | **26,868 tips**, 1,578 columns, model `GTR+FU+G4m` | `default_files/bacteria/bac_ref/` |
| archaeal reference | 1,002 tips: below the 2,000-tip trigger, so not affected (the pipeline places every ASV in both trees and keeps the domain with lower NSTI) | `default_files/archaea/arc_ref/` |
| previous default (v2.0–2.5, "oldIMG", still shipped as `picrust2_pipeline_oldIMG.py`) | one tree with 19,493 bacteria + 406 archaea ≈ 19.9K tips, also > 2,000. EPA-ng ≥ 0.3.3 (Dec 2018) has the bug, so essentially every PICRUSt2 EPA-ng run since the 2020 paper used the affected path | PICRUSt2 wiki, PICRUSt2-SC page |
| query start column (hmmalign output, 1,575–1,660 columns) | no ASV starts at column 0. Medians: ocean (V4) 549; mammal (V4–V5, 400 nt) 1,006; indian (130 nt after 341F) 372 | `results/*_placement.json` |

## 2. Methods

**Builds.** EPA-ng v0.3.8 (tag, commit cff6a47) was built from source twice: stock, and with
`cs581/epang/code/epa-ng-fix.patch` (bug 2 scaler shift plus bug 1 `|=`). Wrappers on `PATH` make
PICRUSt2's `place_seqs.py` call either build (`code/epa_wrap.sh`). **Check:** on ocean, the stock source
build and the bioconda binary give the same best edge, LWR and likelihood for all 1,148/1,148 ASVs.

**Pipeline.** `picrust2_pipeline.py -s seqs -i table -o out -p 4` with all defaults: EPA-ng, both
domains, mp HSP, KO + EC + MetaCyc pathways, unstratified, max NSTI 2. Two changes, neither affecting the
outputs:

1. `util.check_empty_traits` is disabled after the first run (`code/picrust2_pipeline_nocheck.py`). It
   loads every trait table to check that it is non-empty (~8 min, ~14 GB).
2. 8–16 GB of swap was added. EPA-ng needs ~19 GB of virtual memory for this tree; the machine has 15 GB.

**Data.** The PICRUSt2 manuscript repository (github.com/gavinmdouglas/picrust2_manuscript) provides
each validation dataset's 16S ASVs, ASV table and paired HUMAnN2 KO/pathway tables from shotgun
metagenomes:

| dataset | ASVs | samples (paired MGS) | read type |
|---|---|---|---|
| ocean | 1,148 | 6 | V4, ~250 nt |
| mammal | 323 | 11 | ~400 nt |
| indian | 2,384 (in reference alignment) | 108 | 130 nt after 341F |
| hmp | 1,865 | 137 | V3–V5 ~245 nt; placement only (see §3.3) |

The tutorial server (kronos.pharmacology.dal.ca) was unreachable, so the tutorial data was not used.
Blueberry, cameroon and primate were not run for lack of time (~25 min per pipeline run, dominated by KO
HSP on the 26,868-tip tree).

**Accuracy vs metagenomes** (a simplified version of the manuscript protocol, `code/compare_outputs.py`):

- Per sample, Spearman's ρ between predicted and HUMAnN2 KO abundances over KOs that both PICRUSt2's
  current KO table and HUMAnN2 can report. Missing values are set to 0.
- Pathway accuracy is computed the same way over `possible_path/picrust2_path.txt`.
- Presence/absence precision and recall.
- Stock vs patched compared with a paired Wilcoxon test.

**Leave-out with known gene content** (`code/loo_prep.py`, `loo_prune.R`, `run_loo.sh`, `score_loo.py`):

- 300 random GTDB bacterial genomes (seed 1) are removed from the reference tree (castor), the MSA and
  the KO table. The HMM is rebuilt on the pruned MSA (`hmmbuild --dna`; same length, 1,529), because
  `hmmalign --mapali` requires matching files.
- Each held-out genome's own 16S is cut into two query sets:
  - (a) V4 between the 515F/806R primers (~253 nt);
  - (b) the first 130 nt after 341F, mimicking the indian reads (295/300 genomes).
- Queries are placed with PICRUSt2's `place_seqs.py`, and KOs are predicted with `hsp.py` (mp, NSTI).
- Placement error = edges between the best-LWR edge and the edge where the genome sat in the full tree.
  This is computed with xor-hashed bipartitions.
- KO accuracy = presence precision, recall and F1 against the genome's true KO row.

Caveat: the reference is a GTDB *genome* tree, so 16S placement disagrees with it even without any bug.
The median error is ~7 edges for 130-nt reads.

## 3. Results

### 3.1 Placement (stock vs patched EPA-ng, inside PICRUSt2)

PLACEMENT_TABLE

"logL changed" means the best placement's likelihood differs by > 0.001. Stock likelihoods are inflated
by up to 197K (ocean), 316K (mammal) and 216K (indian) log units. Differences go both ways (mammal: patched up to 160K higher for some ASVs), but on
indian the patched likelihood is never higher. "LWR≥0.99" counts placements
EPA-ng reports as essentially certain. Moved ASVs are falsely confident under the stock build (mean LWR
0.70–0.95 vs 0.35–0.66 patched).

### 3.2 Downstream PICRUSt2 outputs

NSTI_TABLE

NSTI goes **up** with the patch for most changed ASVs: stock places reads falsely close to a reference
tip. In ocean, 9 ASVs switch domain (bacteria ↔ archaea) because the bacterial NSTI changes.

### 3.3 Accuracy vs paired metagenomes

ACC_TABLE

HMP_TEXT

### 3.4 Accuracy vs known genome content (leave-out, 300 genomes)

LOO_TABLE

### 3.5 Runtime

EPA-ng wall clock inside PICRUSt2 (4 threads, 26,868-tip tree):

| dataset | stock | patched |
|---|---|---|
| ocean | 101 s | 93 s |
| mammal | 110 s | 98 s |
| indian | 146 s | 127 s |
| leave-out V4 | 88 s | 83 s |
| leave-out V3 | 118 s | 91 s |

A direct back-to-back rerun on ocean gave 110 s vs 111 s. The difference is within run-to-run noise:
loading the tree and the reference CLVs dominates, so bug 1's SIMD loss is not visible at these query
counts. EPA-ng is only ~10% of a ~20–25 min PICRUSt2 run; KO hidden-state prediction is ~60%.

## 4. How this fits a 4-week CS581 project on the EPA-ng bug

PICRUSt2 is the most-cited downstream user of EPA-ng on the affected path, so it makes a good fourth-week
"impact" section. It cannot carry the project alone.

- **Week 1–2 (core):** the bug, the fix, and the RNASim/BSCAMPP accuracy results (cs581/epang).
- **Week 3 (impact on a real pipeline, this pilot extended):**
  - Show that PICRUSt2's default path triggers the bug.
  - Report placement and NSTI changes on the manuscript's validation datasets: all seven if run on a
    bigger machine (hmp needs > 30 GB for the patched build at chunk size 5,000).
  - Report the leave-out accuracy.
  - Optionally add SEPP (PICRUSt2's alternative placer, unaffected) as a third column.
- **Week 4:** write-up. The honest message: "a likelihood bug that doubles placement error on
  RNASim fragments changes 7–43% of PICRUSt2 placements and inflates their confidence. Because
  hidden-state prediction is robust to small placement errors, functional predictions are almost
  unchanged." That is a useful lesson about where placement accuracy matters (phylogenetic placement
  benchmarks, BSCAMPP subtree sizes) and where it does not (gene-content averaging).
- Before reporting upstream, compare with SEPP or with `--no-pre-mask` on one PICRUSt2 dataset. This
  checks that the patched placements match an unaffected code path (here they were checked only
  indirectly, through the RNASim results on cs581/epang).

## 5. Files

- `code/`
  - `epa_wrap.sh`: PATH wrapper.
  - `run_pipeline*.sh`, `picrust2_pipeline_nocheck.py`: pipeline runs.
  - `cmp_jplace.py`, `summ_jplace.py`: placement comparison.
  - `compare_outputs.py`: downstream outputs and MGS accuracy.
  - `loo_*.py`, `loo_prune.R`, `run_loo.sh`, `score_loo.py`: leave-out experiment.
  - `make_tables.py`: tables.
  - `queue*.sh`, `after_loo*.sh`: job order.
- `results/`
  - `*_placement.json`, `*_compare.json`, `*_jplace_stock_vs_fix.tsv.gz`: per-ASV best edge, LWR and
    logL for stock vs patched.
  - `*_accuracy_per_sample.tsv`.
  - `loo_*_summary.json`, `loo_*_scores.tsv`.
