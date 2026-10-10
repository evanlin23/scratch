# Pilot: structure-aware protein MSA from sequence alone (ProstT5-predicted 3Di)

Overnight pilot, 2026-10-10, 4-core CPU machine, no GPU. Code in `code/`, tables in `results/`.
Everything here is reproducible with `code/predict3di.sh`, `code/run_methods.sh`, `code/summarize.py`, `code/rv100.py`.

## Question and source

> FoldMason (Gilchrist et al., Science 2026; bioRxiv 10.1101/2024.08.01.606130) shows structure-based MSA gains beyond the
> twilight zone, but with real structures; Unicore (2025, PMC12203212) runs ProstT5 -> FoldMason without structures but
> reports no MSA-accuracy benchmark ("database creation is Unicore's biggest bottleneck, which requires 3Di conversion by ProstT5").

The FoldMason preprint itself names the follow-up: "we intend to use ProstT5 to predict 3Di directly from amino acid
sequences ... eliminating the need for pre-computed structures".

**Pilot question.** Does aligning ProstT5-predicted 3Di strings (with FoldMason, or with MAFFT L-INS-i using Foldseek's 3Di
substitution matrix) beat MAFFT L-INS-i on low-identity BAliBASE sets, and does predicted-3Di evidence help MAGUS?

PRIOR_ART

## Methods

**Data.** BAliBASE 3.0 (downloaded from lbgi.fr): RV11 (<20% identity; 38 full-length "BB" + 38 truncated "BBS" sets) and
RV12 (20-40%; 44 BB sets; the 44 BBS sets were still running at the deadline, see the end). 120 sets, 4-30 sequences each.
Large sets: RV100 BBA0067, BBA0101, BBA0117, BBA0081 from `cs581/data/balibase_clean` (the MAGUS-paper references).

**3Di prediction.** `foldseek createdb in.fa db --prostt5-model weights --threads 4` (Foldseek 10.941cd33, ProstT5 f16 gguf,
CPU). The database (AA + 3Di, no coordinates) is kept and reused by every 3Di method so all methods see identical 3Di.

**Methods compared (all single-threaded, same inputs; paired per set).**

| name | what |
|---|---|
| `linsi` | MAFFT 7.505 L-INS-i on amino acids (**baseline**) |
| `clustalo`, `mafft` | Clustal Omega 1.2.4 (`--threads 1`), MAFFT default (FFT-NS-2) |
| `fm` | FoldMason 4.dd3c235 `structuremsa` on the ProstT5 database (works without coordinates; refinement needs lDDT, so off) |
| `fm_aa` | control: same FoldMason run with the 3Di term switched off (`--bitfactor-3di 0.01`; 0 is rejected) |
| `linsi3di` | MAFFT L-INS-i on the **predicted 3Di strings** with Foldseek's `mat3di.out` via `--aamatrix`, mapped back to AA |
| `mc2` | M-Coffee (T-Coffee 12.00.7fb08c2 consistency library) of `linsi` + `linsi3di` |
| `mc3di` | M-Coffee of `linsi` + `linsi3di` + `fm` |
| `mcaa` | control: M-Coffee of three AA-only alignments `linsi` + `clustalo` + `mafft` (separates "ensemble" from "3Di") |
| `*_true` | oracle: same pipelines with 3Di from the **experimental PDB structure** (`code/true3di.py` maps each BAliBASE sequence to its best PDB chain; used only on the 73 RV11 sets where every sequence mapped, min coverage >= 0.9) |

**Scoring.** `code/bbscore.py` reimplements BAliBASE `bali_score` (SP and TC on reference core blocks). Checked against the
official `bali_score` binary on 10 sets: SP identical to 3 decimals, TC identical before bali_score's 2-decimal truncation.
**Published-number check:** our L-INS-i on all 76 RV11 sets has SP 0.678; the DIALIGN-TX paper (Subramanian et al. 2008,
AMB 3:6, Table 6) reports 0.671 for MAFFT 6.240 L-INS-i on RV11. Our Clustal Omega RV11 SP 0.608 vs 0.590 in a published
comparison (BMC Bioinf. tables). RV100 sets are scored with FastSP (SPFN+SPFP)/2 like the MAGUS paper.
Significance: two-sided Wilcoxon signed-rank on paired per-set differences.

## Results: BAliBASE RV11/RV12 (SP and TC on core blocks, higher is better)

RESULTS_TABLE

## Results: large RV100 sets and MAGUS evidence

RV100_TABLE

## Runtime

RUNTIME

## Verdict

VERDICT
