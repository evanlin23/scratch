# uDance baseline on M1HF fragmentary data (CS581 pilot)

Date: 2026-10-10. Time-boxed (60 min). R0 ran successfully; R1 was NOT run (it was out of time: R0 took about 22 min wall time single-threaded).

## Version
- github.com/balabanmetin/uDance, master commit `d0747600944f332202ba650584420184f759614f` (2024-04-02), README says v1.6.5.
- Clone: `/opt/udance_work/uDance` (the code is patched locally, see below; nothing was committed).
- Env `udance` (/opt/mm/root/envs/udance): python 3.9, snakemake 7.32.4, raxml 8.2.12, iqtree 2.1.2, fasttree 2.1.10,
  treeshrink 1.3.9, julia 1.7.1 (ArgParse), apples 2.0.11, kmeans1d 0.3.1, ASTRAL 5.17.2 (built from the bundled source), openjdk.

## Install commands
```
export MAMBA_ROOT_PREFIX=/opt/mm/root
git clone https://github.com/balabanmetin/uDance.git /opt/udance_work/uDance
# = install.sh's conda line, + openjdk (for ASTRAL)
/opt/mm/micromamba create -y -n udance -c bioconda -c conda-forge -c smirarab python=3.9 pip newick_utils=1.6 \
  setuptools seqkit=2.1.0 scipy dendropy=4.5.2 pandas=1.3.0 snakemake raxml=8.2.12 iqtree=2.1.2 treeshrink=1.3.9 \
  fasttree=2.1.10 julia=1.7.1 gappa=0.7.1 trimal=1.4.1 raxml-ng tqdist openjdk
/opt/mm/micromamba install -y -n udance -c conda-forge time rsync   # rsync is needed by process_a_marker.sh, not in install.sh
export PATH=/opt/mm/root/envs/udance/bin:$PATH
pip install apples==2.0.11 kmeans1d==0.3.1
cd /opt/udance_work/uDance && julia uDance/correction_multi.jl datasmall/alignments/p0309.fasta >/dev/null  # installs ArgParse
cd uDance/tools/ASTRAL && ./make.sh && java -Djava.library.path=lib/ -jar native_library_tester.jar
```

## Local patches needed for a single-gene, 1-thread run (`udance_patches.diff`)
1. `uDance/create_concat_alignment.sh`: `seqkit concat` fails with "at least 2 files needed" when there is only one gene. The patch copies the single alignment instead.
2. `uDance/PoolAstralWorker.py` and the `udance.smk` blinference rule: ASTRAL-MP aborts with "at least two threads are needed" when `-T 1`. The patch uses the single-threaded `astral.5.17.2.jar` instead and drops `-C -T`.
3. `uDance/PoolAstralWorker.py`: with one gene, TreeShrink can remove a backbone (anchor) taxon from the only gene tree, so ASTRAL aborts on the constraint tree ("SIBE was not seen in main input trees"). The patch prunes the constraint trees to the taxa present in the gene trees. This was triggered once, in partition 5 (SIBE).

## Data / setup (`prep_udance.py`, `run_udance_rep.sh`, `config_template.yaml`)
- `WD/alignments/gene.fasta` = all 1000 seqs of `true_align.fasta` (single gene, DNA).
- Backbone = seqs with ungapped length >= 0.5 x median (R0: median 707.5, giving 522 backbone seqs). `WD/backbone.nwk` = `trees/cache/true_align.bb_ft_0.5.tre`
  (the existing FastTree tree; its leaf set matches these 522 exactly). Queries = the other 478 (shortest is 79 nt ungapped).
- Config = repo `config.yaml` with these changes: `backbone: "tree"`, `resources.cores: 1`, `large_memory: 8000`,
  `refine_config.occupancy: 1` (the default of 2 would drop every non-anchor taxon when there is only one gene). Everything else is at the defaults:
  `infer_config.method: raxml-8`, `numstart: 2`, `numthread: 1`; `trim percent_nongap: 0.05` (trimAl + TAPER correction);
  `backbone_filtering: False`; APPLES-2 FM, `filter 0.2`, `base 25`, `overlap 0.05`; `cluster_size: auto` (which gave 200, so 6 partitions);
  `sublength 100`, `fraglength 75` (the fragment filter in PoolAlignmentWorker is a hard-coded 75 bp); `contract 0.33`; `infer_branchlen: True`.
- Command: `nice -n 10 snakemake --cores 1 --configfile WD/config.yaml --snakefile udance.smk --rerun-incomplete all`, run from the repo dir with
  `OMP/MKL/OPENBLAS_NUM_THREADS=1` and wrapped in GNU `time -v`. Workdir: `/opt/udance_work/R0`.

## Results R0 (true tree /opt/data/mlcache/M1HF/R0/true_tree.tre, scored with treeerr.py on the common leaf set)
| tree | leaves | queries kept | FN | FP |
|---|---|---|---|---|
| udance.updates.nwk (= maxqs) | 922 / 1000 | 401 / 478 | 0.2516 | 0.2524 |
| udance.incremental.nwk | 922 / 1000 | 401 / 478 | 0.2571 | 0.2579 |
| (input backbone FastTree, 522 leaves) | 522 | n/a | 0.1236 | 0.1252 |

Fragments were dropped: 78 taxa are missing (77 queries plus backbone SIBE). Where they were lost:
- The 75 bp fragment filter removed 0.
- APPLES `--exclude` left 3 queries unplaced (475/478 placed).
- TreeShrink removed 7 (1-2 per partition).
- The partition ASTRAL trees cover 990 taxa, but the stitched output has only 922. That means about 68 were lost at stitching, and I did not have time to find the cause. Use this baseline with care.

Note: the FN/FP figures are measured only on the 922 retained leaves, so they are not directly comparable with methods that place all 1000 sequences.

Compute, single thread, `nice 10`, with 4 other heavy jobs on the machine: the run took 4 resumed snakemake invocations because of the failures above. Completed gene trees and placement were reused between them; the early runs that failed at placement_prep and at the rsync step were restarted from scratch and are not counted.
- Summed wall time = 1315 s (21.9 min). Summed CPU (user+sys) = 1315 s. Each partition gene tree (FastTree + TreeShrink + 2 x [IQ-TREE -fast, RAxML-8, IQ-TREE aBayes]) took about 3.5 min.
- Peak RSS (largest single process, GNU time) = 491 MB.
- Per-part `time -v` files: `/opt/udance_work/R0/time_v.part{1,2a,3,4}.txt`. Snakemake logs: `snakemake.part*.log`. uDance log: `udance.log`.

## R1
Not run (out of time). To run it: `TAG=.r1 /home/user/scratch/cs581/fragml2/udance/run_udance_rep.sh R1`. The workdir `/opt/udance_work/R1` is already prepared (528 backbone, 472 queries). Expect about 22 min.

## R1–R4 (run by the main session with `run_udance_rep.sh`, same config, single thread, nice 10)
- R1, R2, R4 finished; scores in `results/udance.md` (`code/udance_score.py`).
- R3 failed in the final `stitch` rule (snakemake "Error in rule stitch"; all 6 partition ASTRAL trees
  exist). Not debugged (time box); R3 is reported as a failure.
