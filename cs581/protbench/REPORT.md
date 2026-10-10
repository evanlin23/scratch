# Does the Clustal-backbone effect generalize beyond the 8 BAliBASE RV100 sets?

Question: MAGUS builds its GCM evidence from 10 backbone alignments of 200 sequences aligned by MAFFT L-INS-i.
On 5 BAliBASE RV100 sets, aligning only these backbones with Clustal Omega lowered MAGUS's error by 1.6-1.9 SP
points. Does that carry over to other protein benchmarks (10AA, HomFam, simulated proteins), when does it help vs
hurt, and does it change tree accuracy?

<!-- RESULTS -->

## Data (sources and licenses)

All protein. One MAGUS draw (draw 0) per dataset: more datasets rather than more draws.

| set | datasets | size | reference | source |
|---|---|---|---|---|
| 10AA | 1GADBL_100, coli_epi_100 (new); the 8 RV100 sets in their **unfiltered** form (with fragments) | 303-807 seqs | structural, all sequences | SALMA/EMMA data release, Illinois Data Bank [IDB-2567453](https://doi.org/10.13012/B2IDB-2567453_V1) (`salma_paper_datasets.zip`, `10aa/`), CC0 |
| HomFam | aat, p450, sdr, adh, blmb, rrm, PDZ, Acetyltransf, rvp, zf-CCHH | 2,000 seqs each (subsample) | Homstrad seeds only (5-20 seqs) | same release, `homfam/` (the 10 largest HomFam families, 15k-94k seqs, Sievers et al. 2011), CC0 |
| SIMMOD / SIMHIGH | 4 replicates each | 1,000 seqs, root length 300 aa | true alignment and true tree | simulated here with AliSim (IQ-TREE 3.1.4) |

Notes on the data:
- **10AA is mostly the same BAliBASE sets.** The "10AA" collection used in UPP/PASTA-era and later Warnow-lab
  papers is the 8 BAliBASE RV100 sets plus 1GADBL_100 and coli_epi_100. Only those two are new data here. The
  8 RV100 sets in this release are the full RV100 files (303-807 sequences, including fragmentary sequences),
  whereas the bbtool sessions use the length-filtered copies in `cs581/data/balibase_clean/` (195-732
  sequences). I queued the unfiltered RV100 sets last, as a "same families, with fragments" variant.
- **HomFam**: clustal.org's full HomFam tarball (`homfam-20110613-25.tar.gz`) returned HTTP 403 from this
  machine, so I used the 10 largest families as distributed with the SALMA/EMMA data. They are far larger than
  5,000 sequences, so each was subsampled to 2,000 sequences: all Homstrad seed sequences plus a uniform random
  sample of the rest (`code/make_homfam.py`, seed 1). Scores are on the estimated alignment restricted to the
  seed sequences (all-gap columns removed), the standard HomFam protocol; with 5-20 seeds per family each
  HomFam score is noisy. HomFam backbones rarely contain seed sequences, so the backbone SPFN/SPFP column is
  computed only on backbones holding >= 2 seeds (`bb_scored` = how many).
- **Simulation** (`/opt/data/sim/SIM{MOD,HIGH}/R{1..4}`): random Yule-Harding tree on 1,000 taxa
  (`iqtree3 -r 1000 -rlen 0.001 MEAN 0.8`, MEAN = 0.06 moderate / 0.10 high), then
  `iqtree3 --alisim -m LG+G4 --length 300 --indel 0.05,0.05 --indel-size POW{1.7/40},POW{1.7/40}`, tree seed
  100r+7, sequence seed 100r+13. Calibration (one MAGUS run each, not in the tables): MEAN 0.03 gave MAGUS
  1.3% error (too easy), MEAN 0.06 gave 11.8%. MEAN 0.10 was not calibrated in advance.

## Methods

Driver: `cs581/code/gcmx/bbtool_bench.py` (`code/run_all.sh`, job list `jobs.txt`), one job at a time, 4 threads,
4-core / 15 GB cloud VM, nothing else heavy running. Per dataset:

- `magus`: MAGUS end to end with the paper's flags (25 subsets, 10 backbones x 200, MCL, minclusters), timed.
- `e2e-clustalo`: MAGUS end to end with its backbones aligned by Clustal Omega 1.2.4 (`--threads 1`, the 10
  backbones run in parallel as MAGUS schedules them), timed. Its own decomposition and backbone draw.
- `merge-mafft`: control, GCM merge only on `magus`'s own subsets and L-INS-i backbones.
- `merge-clustalo`: paired, the same subsets and the same backbone sequence sets realigned by Clustal Omega,
  then the GCM merge only.

To fit the budget I ran only `--tools clustalo --union ''` (no mafft-auto backbones, no L-INS-i+Clustal union).

Change to the driver (this branch): optional 4th job column UNALIGNED; when given, MAGUS aligns it and every
score uses the estimate restricted to the reference sequences (`acc_ref`); TC is now recorded.

Scoring: FastSP; error = (SPFN+SPFP)/2 in %. Differences: Clustal backbones minus L-INS-i backbones, in error
points; W/T/L with |d| < 0.05 points a tie; two-sided Wilcoxon signed-rank (`code/summarize.py`).

Baselines (`code/baselines.py`): MAFFT L-INS-i and Clustal Omega alone, 4 threads, 30-min cap.

Trees (`code/trees.py`, simulated data): FastTree 2 `-lg -gamma` on the true alignment and on each estimate;
normalized FN / FP / RF against the true tree (DendroPy 4.4).
