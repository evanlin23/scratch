Disclosure: AI-assisted (Claude). The code and these instructions were written with Claude for the CS581 project.

# Reproducing the GCM evidence-filter experiment

This package reproduces the core experiment: filtering low-support evidence in MAGUS's merge step (GCM). Only the
alignment graph changes. MAGUS's subsets, backbones, MCL clustering and min-clusters trace stay the same.

| file | what it does |
|---|---|
| `filter_gcm.py` | One merge with one variant (`magus`, `esK`, `fracF`, `vote-hard`, `vote-hard-bb`, `vote-hardT[-bb]`) on a MAGUS work dir. With `--from-scratch` it runs MAGUS end to end with the filter switched on. With `--ref` it scores the result with FastSP. Needs only MAGUS, numpy and scipy (plus java and FastSP for `--ref`). |
| `fetch_bank.sh` | Collects the rep bank tarballs and the logged result rows from the experiment branches into `data/`, checks sha256, and fixes subset order and paths. |
| `run_all.sh` | Recomputes everything from the bank in four restartable stages: `aln`, `bb`, `trees`, `tables`. |
| `tables.py` | Writes `out/TABLES.md` (rebuilt tables) and `out/MISMATCHES.md` (each rebuilt number next to the logged one). |

`data/` and `out/` are regenerable and git-ignored.

## 0. Machine

- Tested on Ubuntu 24.04.5 LTS, 4 cores, 15 GB RAM, about 35 GB free disk (bank 2.1 GB, outputs about 6 GB).
- The largest merge (16S.T, 5,548 sequences) needs about 4 GB RAM.

## 1. Install (fresh machine, about 5 minutes)

```bash
git clone https://github.com/evanlin23/scratch.git && cd scratch
git checkout claude/cs581-repro

sudo apt-get update
sudo apt-get install -y git python3-pip openjdk-17-jre-headless fasttree mafft mcl

# FastSP (alignment SP error); the runs used commit 146c409 (HEAD of the repository, Aug 2021)
sudo mkdir -p /opt/tools && sudo chown "$(id -u)" /opt/tools
git clone https://github.com/smirarab/FastSP.git /opt/tools/FastSP
git -C /opt/tools/FastSP checkout 146c4096d2ec96c1227188b3aed18d09440bf4c9

# MAGUS, pinned to the commit every run used
git clone https://github.com/vlasmirnov/MAGUS.git cs581/code/MAGUS
git -C cs581/code/MAGUS checkout 39041fc8da5dcb44c95e90c212c667cf225ec129
pip install -e cs581/code/MAGUS numpy==2.5.3 scipy==1.18.1   # add --break-system-packages on a system Python
```

`cs581/code/setup.sh` does the same, plus the extra tools the other sub-projects use.

### Versions used for the verification below

| tool | version | note |
|---|---|---|
| MAGUS | git 39041fc (2024-03-25) | |
| MAFFT used by MAGUS | 7.450 (bundled in `MAGUS/magus/tools/mafft`) | the apt MAFFT 7.505 is not used by the merge |
| MCL used by MAGUS | 14-137 (bundled in `MAGUS/magus/tools/mcl`) | |
| FastTree | 2.1.11, double precision (Ubuntu package `fasttree 2.1.11-2`) | called as `FastTree -lg -gamma -quiet` |
| FastSP | git 146c409 | `java -Xmx4g -jar FastSP.jar` |
| Java | OpenJDK 17.0.20 or 21.0.12 | both give identical FastSP output |
| Python / numpy / scipy | 3.13.16 / 2.5.3 / 1.18.1 | the vote mixture is fitted with scipy optimisers |

Set `FASTSP_JAR=/path/FastSP.jar` if FastSP is somewhere else.

## 2. Get the rep bank (about 30 s)

```bash
bash cs581/repro/fetch_bank.sh
```

What this collects:
- Every `cs581/gcmvote/bank/*.tar.gz` (branches `claude/cs581-gcmvote-*`) and `cs581/gcmtrees/bank/*.tar.gz`
  (branches `claude/cs581-gcmtrees*`).
- Each tarball is checked against its `MANIFEST.tsv` sha256 and extracted to `data/bank/<branch>/<REP>/`.
- Every branch's logged result rows go to `data/logged/<branch>/`.

What a replicate holds:
- `inputs/subalignments/` (MAGUS's 25 L-INS-i subset alignments)
- `inputs/backbones/` (MAGUS's 10 L-INS-i backbones, 200 sequences each)
- `true.fasta`, `unaligned.fasta`, and `true_tree.nwk` for the simulated sets
- for the backbone-count reps: `bb5/`, `bb10/`, `bb20/`

**Subset order.** MAGUS numbers subsets in `os.listdir` order, and the merged alignment depends on that order in
about the 6th decimal. For example, BBA0154 `magus` gives 0.2122879 with the logged order and 0.2122887 with
sorted order. `fetch_bank.sh` therefore writes `subsets.json` (relative paths) for every replicate:
- from the run's own `subsets.json` (16 replicates);
- otherwise from the order of the files inside the tarball (tar reads directories in the same order). Where both
  exist, they agree on 10 of 11 replicates.

`subsets.source` says which one was used.

## 3. Quick check: three replicates, logged numbers exactly (about 15 minutes)

```bash
cd cs581/repro
B=data/bank
python3 filter_gcm.py $B/gcmvote-h3/BBA0154 magus        --out out/check/BBA0154.magus.fasta --ref $B/gcmvote-h3/BBA0154/true.fasta
python3 filter_gcm.py $B/gcmvote-h3/BBA0154 es4          --out out/check/BBA0154.es4.fasta   --ref $B/gcmvote-h3/BBA0154/true.fasta
python3 filter_gcm.py $B/gcmvote-h3/BBA0154 vote-hard-bb --out out/check/BBA0154.hbb.fasta   --ref $B/gcmvote-h3/BBA0154/true.fasta
python3 filter_gcm.py $B/gcmvote-h4/1000L3_R0 magus        --out out/check/1000L3.magus.fasta --ref $B/gcmvote-h4/1000L3_R0/true.fasta
python3 filter_gcm.py $B/gcmvote-h4/1000L3_R0 es4          --out out/check/1000L3.es4.fasta   --ref $B/gcmvote-h4/1000L3_R0/true.fasta
python3 filter_gcm.py $B/gcmvote-h4/1000L3_R0 vote-hard-bb --out out/check/1000L3.hbb.fasta   --ref $B/gcmvote-h4/1000L3_R0/true.fasta
python3 filter_gcm.py $B/gcmvote-h9b/16S.M_R0_B10 magus        --out out/check/16SM.magus.fasta --ref $B/gcmvote-h9b/16S.M_R0_B10/true.fasta
python3 filter_gcm.py $B/gcmvote-h9b/16S.M_R0_B10 vote-hard-bb --out out/check/16SM.hbb.fasta   --ref $B/gcmvote-h9b/16S.M_R0_B10/true.fasta
```

Each command prints one JSON row. `avgErr` must equal the logged value to every printed digit:

| replicate (data type) | variant | expected avgErr | logged in |
|---|---|---|---|
| BBA0154 (protein, BAliBASE) | magus | 0.21228786278830686 | `claude/cs581-gcmvote-h3`: `cs581/gcmvote/results_h3/BBA0154/vote.jsonl` |
| | es4 | 0.21001288748547298 | same |
| | vote-hard-bb | 0.18909024828356435 | same |
| 1000L3_R0 (DNA, ROSE) | magus | 0.1186143423904078 | `claude/cs581-gcmvote-h4`: `cs581/gcmvote/results_h4/1000L3_R0/vote.jsonl` |
| | es4 | 0.11549729184680912 | same |
| | vote-hard-bb | 0.12379252207291541 | same |
| 16S.M_R0, B = 10 (real rRNA) | magus | 0.13034605651232317 | `claude/cs581-gcmvote-h9b`: `cs581/gcmvote/results_h9b/sens.jsonl` |
| | vote-hard-bb | 0.12365449279700502 | same |

## 4. Everything from the bank

```bash
bash cs581/repro/run_all.sh                # = run_all.sh aln bb trees tables; JOBS=4 merges at a time
```

The run is restartable: rerun the same command after an interruption. It writes:
- `out/rows/<branch>/<REP>/<variant>@B<B>.json`: one scored merge;
- `out/aln/...`: merged alignments;
- `out/trees/...`: FastTree trees;
- `out/TABLES.md`, `out/MISMATCHES.md`.

Measured on the 4-core machine:

| stage | jobs | wall time |
|---|---|---|
| aln (every replicate × magus / es4 / vote-hard-bb, B = 10) | RUNTIME_ALN | |
| bb (h6, h6b, h9, h9b × 6 variants × B = 5 / 10 / 20) | RUNTIME_BB | |
| trees (FastTree on true / magus / es4 / vote-hard-bb, SIMHIGH and 1000M1 reps) | RUNTIME_TREES | |
| tables | | < 1 min |

Single merges take about 10–60 s on 1,000-sequence replicates and about 8–10 min on 16S.3 / 16S.T. FastTree
`-lg -gamma` takes about 10–15 min per 1,000-taxon protein alignment.

What `out/TABLES.md` contains:
- training and held-out SP error for magus / es4 / vote-hard-bb, by data type (proteins, simulated DNA/RNA,
  real rRNA), with W/T/L and Wilcoxon;
- the remaining banked replicates;
- the check-in-12 pooled held-out sets;
- the backbone-count table (B = 5 / 10 / 20);
- FastTree nRF on the SIMHIGH and 1000M1 reps.

`out/MISMATCHES.md` lists every rebuilt merge and tree next to the logged row of the same MAGUS draw. It also
recomputes the gcmvote REPORT's summary cells from the logged raw rows.

Known limits:
- The gcmvote main branch's own MAGUS draws were never banked. These are the training sets and 12 of the 14
  held-out sets in REPORT.md §4. Their numbers can be re-aggregated from logged rows but not re-merged; the bank
  has other draws of those datasets. See `REPRO_STATUS.md`.
- Fixed-fraction / es rows that were logged by `gcmgen/code/gg.py` (h6, h6b, h9, h9b, the gcmtrees branches) ran
  MAGUS with its native dict graph. Its graph file has a different line order, so MCL can land on a slightly
  different clustering than `filter_gcm.py` / `vote.py`. These rows are marked `gg.py` in MISMATCHES.md.

## 5. Without the bank: from raw sequences

MAGUS is unseeded: every end-to-end run draws a new decomposition and new backbone sequence sets. A rerun from raw
sequences therefore reproduces the *experiment*, not the logged digits; the bank exists for the digits.

### 5.1 Data sources

| data | source | DOI / URL | used as |
|---|---|---|---|
| ROSE 1000-taxon DNA (1000L1–L3, M1–M4, S1–S3), RNASim 1000, 16S.M / 16S.3 / 16S.T (Gutell CRW) | MAGUS paper data, Illinois Data Bank IDB-2643961, `Datasets.zip` (~390 MB) | doi:10.13012/B2IDB-2643961_V1; file `https://databank.illinois.edu/datafiles/u373n/download` | ROSE: `ROSE/1000XX/R<i>/rose.aln.true.fasta` with true tree `rose.tt`. RNASim: `RNASim/1000/R<i>/true_align.txt`. 16S: `Gutell/16S.X/R0/true_align_clean.txt` |
| BAliBASE RV100 (BBA0039, 0067, 0081, 0101, 0117, 0134, 0154, 0190), length-filtered as in the MAGUS paper (sequences more than 20 % from the median length removed) | MAGUS paper data (as above); BAliBASE: Thompson et al. 1999 doi:10.1093/bioinformatics/15.1.87, BAliBASE 3: doi:10.1002/prot.20527 | | committed in `cs581/data/balibase_clean/RV100_<NAME>.fasta` |
| HomFam (blmb, aat, Acetyltransf, PDZ) | SALMA/EMMA data release, Illinois Data Bank IDB-2567453, `salma_paper_datasets.zip`, folder `homfam/<family>/`; HomFam: Sievers et al. 2011 doi:10.1038/msb.2011.75 | doi:10.13012/B2IDB-2567453_V1 | `python3 cs581/protbench/code/make_homfam.py FAMILY_DIR OUT_DIR 2000 1` gives `ref.fasta` (Homstrad seeds = reference) and `unaln.fasta` (seeds + random family members, 2,000 total); scored on the seeds only |
| AliSim proteins SIMMOD, SIMHIGH (1,000 taxa) | IQ-TREE 3.1.4 AliSim (Ly-Trong et al. 2022, doi:10.1093/molbev/msac092) | `micromamba create -n bio -c conda-forge -c bioconda iqtree=3.1.4` | `cs581/gcmtrees/code/gen.sh` (branch `claude/cs581-gcmtrees`), below |

AliSim settings, replicate r (tree seed 100·r + 7, sequence seed 100·r + 13; MEAN = 0.06 for SIMMOD, 0.10 for SIMHIGH):

```bash
iqtree3 -r 1000 tree.nwk -rlen 0.001 $MEAN 0.8 -seed $((100*r+7)) -redo
iqtree3 --alisim sim -t tree.nwk -m LG+G4 --length 300 --indel 0.05,0.05 \
        --indel-size POW{1.7/40},POW{1.7/40} -seed $((100*r+13)) -af fasta -redo
# sim.fa = true alignment (true.fasta), tree.nwk = true tree (true_tree.nwk)
```

### 5.2 A fresh MAGUS draw with the filter switched on (paper settings: 25 subsets, 10 L-INS-i backbones × 200)

```bash
python3 cs581/repro/filter_gcm.py --from-scratch unaligned.fasta vote-hard-bb --subsets 25 \
        --out run/vote-hard-bb.fasta --ref true.fasta --keep-work
```

- Use `--subsets 100` for 16S.3 / 16S.T, as the logged runs did.
- `unaligned.fasta` is the reference with gaps removed. For HomFam it is `unaln.fasta`, and `--ref ref.fasta`
  scores on the seeds.
- The filter is applied inside MAGUS's own graph construction (`filter_gcm.install_filter` replaces
  `graph_builder.buildMatrix`). Everything else is MAGUS's own code.

A paired experiment needs one MAGUS draw shared by all variants:
1. Run `--from-scratch ... magus --keep-work` once.
2. Turn its work dir into a bank-style replicate. The subsets are in `run/magus.fasta.work/subalignments/` and
   the backbones are `run/magus.fasta.work/graph/backbone_*_mafft.txt`:

   ```bash
   W=run/magus.fasta.work; R=myrep
   mkdir -p $R/inputs/subalignments $R/inputs/backbones
   cp $W/subalignments/*.txt $R/inputs/subalignments/
   cp $W/graph/backbone_*_mafft.txt $R/inputs/backbones/
   cp run/magus.fasta.subsets.json $R/subsets.json     # MAGUS's subset order (only the file names are used)
   cp true.fasta $R/
   ```

3. Run every variant merge-only on it (section 3 commands). A merge-only run on such a replicate reproduces the
   end-to-end run exactly. Checked on BBA0154: `--from-scratch ... vote-hard-bb` gave avgErr 0.19289341027677787
   (a new draw, 26 min at `-np 2`), and the merge-only rerun on the replicate built this way gave the same value
   in 11 s.

The original pipeline is still available: `gcmx.bbtool_bench` → `protcons/code/pc.py rep` →
`gcmvote/code/run.py`, as described in `gcmgen/code/fresh.sh`.
