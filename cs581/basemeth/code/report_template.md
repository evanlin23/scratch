# MAGUS (and PASTA) with different subset aligners: pilot report

Projects addressed (instructor list, `cs581/notes/CS581-project-suggestions.txt`): (a) "study MAGUS with different
base methods", (b) "Prank or ProbCons within PASTA or MAGUS", (c) "Regressive vs PASTA and MAGUS". Pre-registration:
`PREREG.md` (committed before any result). Raw rows: `results/*.jsonl`; all tables: `results/tables.md`
(`code/summarize.py`). One ~6 h session on a 4-core VM.

## Verdict (numbers first)

__VERDICT__

## Main table: MAGUS final SP error (%) by subset aligner

Merge-only: the same 25 MAGUS subsets and the same 10 L-INS-i backbones in every column; only the subset aligner
changes. Bold = best in the row. `linsi` = MAGUS's own subset command re-run (the control).

__MAINTABLE__

## Paired tests vs L-INS-i subsets (SP-error points; negative = better than MAGUS's default)

__PAIRED__

## Pre-registered primary test

__PRIMARY__

## Why: subset accuracy carries through to the final alignment

__SUBSET__

GCM keeps each subset alignment intact (the final alignment restricted to a subset is that subset's alignment), so
a better subset aligner should give a better MAGUS alignment, and it does in direction almost everywhere: on
BAliBASE MUSCLE5, ProbCons and Clustal Omega align MAGUS's subsets 0.2-0.4 points better than L-INS-i and the
final alignment improves by 0.3-0.5; on the simulated sets the same tools align the subsets 1-10 points worse and
the final alignment is 1-9.5 points worse. G-INS-i is the only method that is never clearly worse on subsets.

## Runtime

__COST__

All jobs shared the 4-core VM with other jobs of this pilot (MAGUS, PASTA), so wall times are inflated; CPU
seconds are the fair comparison. In end-to-end MAGUS the subsets are a minor part of the cost (MAGUS's own
backbones are 10 L-INS-i runs on 200 sequences each; the merge takes 1-54 s). Swapping in MUSCLE5 or ProbCons
multiplies the subset stage by 4-10x; FAMSA, Kalign and Clustal Omega make it nearly free, but they are the
least accurate.

## Trees (simulated proteins)

FastTree 2.1.11 `-lg -gamma`; normalized RF (%) vs the true 1,000-taxon tree (`code/trees.py`).

__TREES__

## PASTA 1.8.3 with its built-in aligners, and whole-dataset baselines

__BASE__

## Nucleotide side (secondary)

__DNA__

__METHODS__

## Data

| set | datasets | size | reference | MAGUS inputs |
|---|---|---|---|---|
| BAliBASE RV100 (length-filtered, `cs581/data/balibase_clean`) | BBA0039, 0067, 0101, 0154, 0190 | 274-732 seqs | structural, all sequences | cached MAGUS runs (`cs581/experiments/runs/*_R0`) |
| HomFam (SALMA/EMMA release IDB-2567453, 2,000-seq subsamples, seed 1, `code/make_homfam.py`) | aat, sdr, rrm, PDZ, Acetyltransf, zf-CCHH | 2,000 seqs | Homstrad seeds only (6-20 seqs) | new MAGUS runs |
| AliSim proteins (as in the protbench branch: Yule tree, LG+G4, 300 aa, indel 0.05/0.05 POW 1.7/40; seeds 100r+7 / 100r+13) | SIMMOD R1-R2 (mean branch 0.06), SIMHIGH R1-R2 (0.10) | 1,000 seqs | true alignment and tree | new MAGUS runs |
| ROSE 1000M2, 1000L1 R0; RNASim 1000 R0; 16S.M R0 (MAGUS paper data) | 4 | 740-1,000 | true / CRW | cached MAGUS runs |

AliSim R3/R4 (in PREREG scope) were dropped for time: each MAGUS run took ~20-40 min under load.
I could not verify that the regenerated AliSim files are byte-identical to protbench's (no checksums were kept);
same commands and seeds, IQ-TREE 3.1.4, and MAGUS's error is in the same range (SIMMOD_R1 12.7% here vs 13.9% there;
MAGUS draws differ).

## Plan if pursued (4 weeks)

__PLAN__

## Risks / caveats

__RISKS__
