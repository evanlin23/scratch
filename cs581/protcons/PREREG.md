# Pre-registration: does consistency-filtered GCM evidence generalize?

Written and committed before any new data was run (CS581, UIUC, Fall 2026). Everything below is frozen;
deviations will be reported as such in REPORT.md.

## Methods under test (implemented exactly as in `cs581/bbevidence/code/bbe.py`)

All are merge-only: MAGUS's own 25 L-INS-i subset alignments and its own 10 backbone sequence sets (200
sequences each) are kept; only the evidence files given to GCM change. Merge flags as in bbe.py
(`--graphclustermethod mcl --graphtracemethod minclusters --graphtraceoptimize false -f 4`, `-np 4`).

- **Primary:** `linsi|cons0.7`. MAGUS's own 10 L-INS-i backbones; every backbone column whose cross-backbone
  consistency is < 0.7 is turned into an insertion column (lowercase letters, `.` for gaps).
- **Secondary 1:** `linsi&fftns2-op3`. Per backbone, only the residue pairs that both L-INS-i and MAFFT FFT-NS-2
  `--op 3` (same sequences) align.
- **Secondary 2:** `(linsi+clustalo)|cons0.7`, written `linsi+clustalo|cons0.7` in bbe.py's grammar. MAGUS's 10
  L-INS-i backbones plus the same 10 sequence sets aligned by Clustal Omega 1.2.4 (`--threads 1`), 20 backbones,
  consistency-masked at 0.7.
- **Control:** `linsi`, the GCM merge on MAGUS's own backbones (should reproduce MAGUS's own output).

## Gate (reference-free regime test, bbevidence REPORT §5 / §4.3)

Statistic: `support_a_only` of `cs581/bbevidence/code/pairdiff.py` with A = `linsi`, B = `clustalo`: among the
cross-subset residue pairs that an L-INS-i backbone aligns but the Clustal Omega alignment of the same sequence
set does not, the fraction of (other L-INS-i backbone containing both residues) cases in which that other
backbone also aligns them. It uses no reference (I recompute it without the reference file; equivalence to
pairdiff.py is checked on a dataset with a full reference).

The source gives a range (helped ≤ 0.59, hurt ≥ 0.64), so the threshold is the midpoint, **τ = 0.615**:
**filter (use the method) if support < 0.615; otherwise keep MAGUS's own merge.**

## Held-out data (priority order)

1. Simulated proteins with true trees: AliSim LG+G4 with indels as in `cs581/protbench` (SIMMOD, SIMHIGH; R1–R4
   each; same seeds), 1,000 sequences, MAGUS with the paper's flags (25 subsets, 10 × 200 backbones).
2. HomFam (10 families, 2,000-sequence subsamples via `protbench/code/make_homfam.py`, seed 1; scored on the
   Homstrad seeds) and 10AA 1GADBL_100 / coli_epi_100.
3. Fresh MAGUS draws on the 8 BAliBASE RV100 sets (`cs581/data/balibase_clean`), one each.
4. Nucleotide controls (cached ROSE / RNASim / 16S MAGUS inputs), where the gate should say "do not filter".

One MAGUS draw per dataset. Every comparison is paired on that draw's subsets and backbone sets.

## Endpoints

- **Primary:** paired Δ error = error(method) − error(`linsi` control), error = (SPFN + SPFP)/2 in % (FastSP;
  HomFam on the seed sequences only), on held-out **protein** data (families 1–3 pooled), two-sided Wilcoxon
  signed-rank test; W/T/L with |Δ| < 0.05 points a tie. Success for the primary method = mean Δ < 0 and
  p < 0.05 on the pooled protein data.
- **Secondary:**
  - the same for the two secondary methods, and per dataset family;
  - Δ SPFN, Δ SPFP, Δ TC;
  - simulated data: FastTree 2 `-lg -gamma` trees on the control (MAGUS) alignment, the primary method's alignment
    and the true alignment; normalized RF to the true tree, paired Δ;
  - gate accuracy: predicted "help" (support < 0.615) vs observed help (Δ of the primary method < 0), over all
    held-out datasets including the nucleotide controls; and the gated policy (method if gate says filter, else
    control) evaluated as its own Δ;
  - runtime of each recipe (extra seconds over MAGUS, and end to end where measured).

## Analysis notes fixed in advance

- If a run fails for a dataset, it is reported as missing, not replaced.
- No threshold, method parameter or dataset is changed after seeing held-out results; any extra exploratory
  analysis is labelled exploratory.
