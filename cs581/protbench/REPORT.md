# Does the Clustal-backbone effect generalize beyond the 8 BAliBASE RV100 sets?

Question: MAGUS builds its GCM evidence from 10 backbone alignments of 200 sequences aligned by MAFFT L-INS-i.
On 5 BAliBASE RV100 sets, aligning only these backbones with Clustal Omega lowered MAGUS's error by 1.6-1.9 SP
points. Does that carry over to other protein benchmarks (10AA, HomFam, simulated proteins), when does it help vs
hurt, and does it change tree accuracy?

## Verdict (numbers first)

22 protein datasets, one MAGUS draw each (2 new 10AA sets, 2 unfiltered RV100 sets, 10 HomFam families at
2,000 sequences, 8 AliSim simulations). Difference = MAGUS with Clustal Omega backbones minus MAGUS, in SP-error
points ((SPFN+SPFP)/2 x 100); negative = Clustal backbones better. W/T/L = Clustal-backbone wins/ties/losses.

| set | n | MAGUS error | d merge-only (paired) | W/T/L | d end-to-end | W/T/L | MAGUS wall | e2e wall ratio |
|---|---|---|---|---|---|---|---|---|
| 10AA new (1GADBL, coli_epi) | 2 | 3.7 | +0.23 | 0/1/1 | +0.22 | 0/0/2 | 80 s | 0.38 |
| RV100 unfiltered (BBA0039, BBA0067) | 2 | 16.2 | -0.35 | 1/0/1 | -0.09 | 1/0/1 | 972 s | 0.29 |
| HomFam (2,000 seqs, seed-scored) | 10 | 18.2 | +0.76 | 2/3/5 | +1.11 | 4/1/5 | 201 s | 0.49 |
| SIMMOD (AliSim, MAGUS ~11.6%) | 4 | 11.6 | **+3.01** | 0/0/4 | **+2.93** | 0/0/4 | 805 s | 0.35 |
| SIMHIGH (AliSim, MAGUS ~23.6%) | 4 | 23.6 | **+4.08** | 0/0/4 | **+3.55** | 0/0/4 | 1,176 s | 0.38 |
| all biological (10AA + RV100 + HomFam) | 14 | | +0.52 (p = 0.39) | 3/4/7 | +0.81 (p = 0.60) | 5/1/8 | | |
| **all** | 22 | | **+1.62 (p = 0.004)** | 3/4/15 | **+1.70 (p = 0.016)** | 5/1/16 | | 0.38 |

(p = two-sided Wilcoxon signed-rank; with n = 4 the smallest possible p is 0.125, so 4/4 losses is the strongest
statement the simulated sets allow per level.)

**Does the effect generalize? No.** Outside BAliBASE, Clustal Omega backbones never gave MAGUS a meaningful,
consistent gain:

- **Simulated proteins (LG+G4 with indels): a clear, consistent loss**, 8/8 paired merges worse, by 1.6-5.1
  points, at both difficulty levels; both SPFN and SPFP get worse (SIMMOD +3.0 / +3.0, SIMHIGH +4.4 / +3.7
  SPFN / SPFP points). Clustal Omega's backbones themselves are terrible here: backbone SPFP 25-29% vs 9-13%
  for L-INS-i on SIMMOD, 43-47% vs 20-25% on SIMHIGH. Curiously TC goes *up* with Clustal backbones on the
  simulations (+7 to +14 points), so the final alignment gets more columns exactly right while aligning many
  more residue pairs wrongly; I did not chase this.
- **HomFam: neutral on average, with one large loss.** Paired merge-only: mean +0.76, 2 wins, 3 ties, 5 losses;
  without PDZ (+7.0 points) the mean is +0.06. HomFam is scored on only 5-20 seed sequences per family, so
  single-family numbers are noisy; the end-to-end differences swing by up to +-6 points (rvp +6.0 end-to-end
  vs 0.00 paired; sdr -3.0 vs +0.06), which is draw noise from MAGUS's random decomposition, not the backbones.
- **10AA new sets** (1GADBL, coli_epi): easy for MAGUS (3-4% error); Clustal backbones tie or lose slightly.
- **Unfiltered RV100 (same families as the pilot, fragments kept)**: BBA0067 improves (-0.83 paired, all from
  SPFP: -3.0 SPFP, +1.3 SPFN), BBA0039 is a near-tie (+0.12). Only 2 of the 8 finished (~30 min each, out of
  time), so this neither confirms nor refutes the BAliBASE result.

**When does it help vs hurt?** The paired effect tracks how much worse Clustal's backbones are than
L-INS-i's: Spearman rho = 0.79 (p = 1e-5, n = 22) between the backbone error gap (Clustal minus L-INS-i
backbone (SPFN+SPFP)/2, on the reference sequences the backbones contain) and the final paired difference.
When Clustal's backbones are about as accurate as L-INS-i's (10AA, RV100, most HomFam families: gap within
about +-7 points; zf-CCHH, +12.7, is the exception at +0.86 paired), the swap is roughly neutral and can help a little through SPFP, which is the pilot's mechanism
(fewer aligned pairs, fewer wrong ones: Clustal backbones give GCM 0-22% fewer residue pairs). When Clustal's
backbones are much worse (simulated indel-rich data: gap +21 to +29 points; PDZ: +7.7), MAGUS gets clearly
worse. GCM does not "filter" bad evidence: it follows the backbones.

Caveat on the HomFam backbone gap: only backbones holding >= 2 seed sequences can be scored (1-5 of 10 per
family), so that column is rough for HF; the correlation is driven mainly by the simulated sets.

**Runtime.** End-to-end MAGUS with Clustal backbones took 0.38x MAGUS's wall-clock on aggregate (0.21-0.90x per
dataset; biggest savings on long sequences and many-sequence backbones, smallest on the very short zf-CCHH
and rvp). Realigning 10 backbones with Clustal took 7-36 s vs minutes for L-INS-i. The speedup is real and
robust; the accuracy gain is not.

**Trees.** On the simulated data the Clustal-backbone swap barely moves FastTree error: paired merge-only +0.05 RF
points (SIMMOD) and +0.9 (SIMHIGH, 4/4 worse); end-to-end differences are dominated by MAGUS's random draw.

**Bottom line for the project.** The BAliBASE gain looks benchmark-specific (BAliBASE references score only
core blocks, where Clustal's lower-SPFP backbones are not penalized), and on simulated proteins with a known
true alignment the swap costs 3-4 points. If the project keeps this direction, frame it as a speed/accuracy
trade-off or as "which backbone aligner suits which data", and test the union of L-INS-i and Clustal backbones
(not run here) rather than a replacement.

### Trees (simulated data only)

FastTree 2.1 `-lg -gamma` on each alignment; normalized RF (%) against the true 1,000-taxon tree (all trees
binary, so FN = FP = RF). `merge-mafft` has exactly MAGUS's SP scores and length, so MAGUS's tree is the
paired-merge baseline. Raw rows: `results/trees.jsonl`.

| dataset | true aln | MAGUS | e2e Clustal-bb | merge Clustal-bb | d RF e2e | d RF merge |
|---|---|---|---|---|---|---|
| SIMMOD_R1 | 6.52 | 5.92 | 6.02 | 6.52 | +0.10 | +0.60 |
| SIMMOD_R2 | 6.52 | 7.22 | 6.52 | 7.22 | -0.70 | +0.00 |
| SIMMOD_R3 | 6.12 | 7.52 | 7.02 | 7.12 | -0.50 | -0.40 |
| SIMMOD_R4 | 5.52 | 6.62 | 6.82 | 6.62 | +0.20 | +0.00 |
| SIMHIGH_R1 | 6.12 | 13.54 | 12.04 | 14.74 | -1.50 | +1.20 |
| SIMHIGH_R2 | 5.22 | 6.42 | 9.93 | 8.22 | +3.51 | +1.80 |
| SIMHIGH_R3 | 6.42 | 11.03 | 8.93 | 11.13 | -2.10 | +0.10 |
| SIMHIGH_R4 | 6.62 | 9.03 | 10.93 | 9.63 | +1.90 | +0.60 |

| level | comparison | n | mean d RF | W/T/L | p |
|---|---|---|---|---|---|
| SIMMOD | end-to-end | 4 | -0.23 | 2/0/2 | 0.63 |
| SIMMOD | merge-only (paired) | 4 | +0.05 | 1/2/1 | 1 |
| SIMHIGH | end-to-end | 4 | +0.45 | 2/0/2 | 0.88 |
| SIMHIGH | merge-only (paired) | 4 | +0.93 | 0/0/4 | 0.125 |

**Trees follow the alignment error only weakly.** At moderate divergence the 3-point SP loss from Clustal
backbones does not change tree error at all (within +-0.7 RF points; MAGUS's tree is already as good as the
true-alignment tree, 5.9-7.5% vs 5.5-6.5%). At high divergence the paired merge with Clustal backbones gives a
worse tree on 4/4 replicates, but only by +0.1 to +1.8 RF points (mean +0.9), while end to end the
decomposition draw dominates (-2.1 to +3.5). Note that MAGUS on SIMHIGH loses 0.2-7.4 RF points relative to
the true alignment, so alignment error does matter for trees there; it is the Clustal-vs-L-INS-i backbone
choice that moves trees little. (The pilot's SPFP gains on BAliBASE have no tree counterpart to test; there
is no true tree.)


### Per-dataset results

Full tables (including TC, wall-clock and backbone SPFP): `results/tables.md`; raw rows: `results/bbtool.jsonl`.

| dataset | seqs | ref seqs | MAGUS err | e2e Clustal-bb err | d e2e | d merge (paired) | MAGUS wall s | e2e wall s | bb SPFP L-INS-i / Clustal |
|---|---|---|---|---|---|---|---|---|---|
| 10AA_1GADBL | 561 | 561 | 3.14 | 3.28 | +0.13 | -0.02 | 129 | 48 | 3.1 / 3.6 |
| SIMMOD_R1 | 1000 | 1000 | 13.89 | 18.17 | +4.28 | +4.55 | 882 | 310 | 12.5 / 27.7 |
| HF_aat | 2000 | 10 | 13.67 | 14.52 | +0.84 | -0.03 | 500 | 210 | 22.2 / 25.0 |
| SIMHIGH_R1 | 1000 | 1000 | 24.55 | 29.60 | +5.05 | +4.76 | 912 | 421 | 25.1 / 46.1 |
| 10AA_coliepi | 320 | 320 | 4.32 | 4.63 | +0.31 | +0.48 | 31 | 14 | 3.2 / 4.0 |
| HF_p450 | 2000 | 12 | 21.47 | 21.23 | -0.24 | -0.34 | 717 | 326 | 11.9 / 13.7 |
| SIMMOD_R2 | 1000 | 1000 | 10.74 | 12.21 | +1.48 | +2.43 | 840 | 288 | 9.1 / 25.8 |
| HF_sdr | 2000 | 13 | 25.60 | 22.56 | -3.04 | +0.06 | 147 | 79 | 18.5 / 22.3 |
| SIMHIGH_R2 | 1000 | 1000 | 21.48 | 23.11 | +1.63 | +1.87 | 1252 | 426 | 20.3 / 42.5 |
| HF_adh | 2000 | 5 | 1.03 | 1.03 | +0.00 | +0.00 | 105 | 55 | 0.5 / 0.1 |
| SIMMOD_R3 | 1000 | 1000 | 10.85 | 11.33 | +0.47 | +1.58 | 738 | 241 | 10.1 / 25.5 |
| HF_blmb | 2000 | 6 | 27.47 | 25.61 | -1.86 | +0.70 | 301 | 160 | 2.4 / 3.3 |
| SIMHIGH_R3 | 1000 | 1000 | 24.93 | 28.03 | +3.10 | +4.56 | 1298 | 464 | 25.0 / 45.6 |
| HF_rrm | 2000 | 20 | 20.87 | 21.26 | +0.39 | +0.06 | 41 | 26 | 14.1 / 19.7 |
| HF_PDZ | 2000 | 6 | 11.16 | 20.21 | +9.06 | +7.03 | 61 | 35 | 6.2 / 12.7 |
| SIMMOD_R4 | 1000 | 1000 | 10.93 | 16.44 | +5.51 | +3.47 | 760 | 286 | 10.3 / 28.8 |
| HF_Acetyltransf | 2000 | 6 | 32.35 | 33.49 | +1.14 | -0.76 | 96 | 52 | 29.6 / 27.0 |
| SIMHIGH_R4 | 1000 | 1000 | 23.53 | 27.95 | +4.42 | +5.11 | 1244 | 500 | 24.8 / 47.2 |
| HF_rvp | 2000 | 6 | 14.33 | 20.33 | +5.99 | +0.00 | 30 | 23 | 29.3 / 23.5 |
| HF_zf-CCHH | 2000 | 15 | 13.99 | 12.82 | -1.16 | +0.86 | 12 | 10 | 6.9 / 20.0 |
| 10AAfull_BBA0039 | 807 | 807 | 6.69 | 7.02 | +0.34 | +0.12 | 624 | 129 | 7.8 / 7.5 |
| 10AAfull_BBA0067 | 410 | 410 | 25.68 | 25.17 | -0.52 | -0.83 | 1320 | 435 | 28.7 / 28.0 |

`merge-mafft` (the control: GCM merge on MAGUS's own subsets and backbones) reproduced MAGUS's score on all 22
datasets to within 0.01 error points (identical on 16; largest difference BBA0067, 25.68 vs 25.69), so the paired difference isolates the
backbone aligner.

### What was not done

- Only 2 of the 8 unfiltered RV100 sets (out of time; ~30 min per set end to end).
- No mafft-auto backbones and no L-INS-i + Clustal union (`--tools clustalo --union ''`), to fit the budget.
- No stand-alone MAFFT L-INS-i / Clustal Omega baselines and no PASTA runs (out of time; `code/baselines.py` is
  ready). FastTree on the merge-mafft alignments was skipped (same SP scores as MAGUS).
- No additional BAliBASE 3 reference sets.
- One MAGUS draw per dataset: the end-to-end differences contain decomposition noise (see HomFam); the
  paired merge-only differences do not.


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
  sequences). I queued the unfiltered RV100 sets last, as a "same families, with fragments" variant; only BBA0039 and
  BBA0067 finished.
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
  1.3% error (too easy), MEAN 0.06 gave 11.8%. MEAN 0.10 was not calibrated in advance and gave 21.5-24.9%.

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
