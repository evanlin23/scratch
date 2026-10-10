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

## Prior art and novelty check (searched 2026-10-10)

Searches (web, bioRxiv/PMC/GitHub): "ProstT5 predicted 3Di multiple sequence alignment benchmark BAliBASE FoldMason";
"3Di ProstT5 MSA accuracy HOMSTRAD"; "FoldMason Science 2026 ProstT5 sequence-only"; "Santus nf-core multiplesequencealign";
"learnMSA2 / vcMSA"; "Muscle-3D ProstT5"; "Unicore ProstT5 FoldMason"; "3Di substitution matrix MAFFT --aamatrix";
"structure-aware MSA 3Di language model 2026".

| work | what it does with predicted 3Di | MSA accuracy vs a reference? |
|---|---|---|
| ProstT5 (Heinzinger et al., NAR Genom Bioinf 2024, lqae150) | predicts 3Di from sequence; benchmarked for remote-homology search (SCOPe40) | no MSA benchmark |
| FoldMason (Gilchrist, Mirdita, Steinegger, Science 391:485, 2026; bioRxiv 10.1101/2024.08.01.606130) | preprint: "we intend to use ProstT5 to predict 3Di directly from amino acid sequences"; current release (4.dd3c235) exposes `--prostt5-model` in `easy-msa`/`createdb` | benchmarks (HOMSTRAD, BAliBASE-type references, lDDT) use real/AlphaFold structures; I found no ProstT5-input benchmark. Could not read the final Science text (paywall / bioRxiv 429): **student should check its supplement before the proposal.** |
| Unicore (Kim, Park, Steinegger, GBE 2025, evaf109; PMC12203212) | ProstT5 -> Foldseek clustering -> FoldMason per core gene -> AA alignment -> IQ-TREE | tree congruence only, no MSA accuracy |
| Garg & Hochberg, "A general substitution matrix for structural phylogenetics" (MBE 2025, msaf124; PMC12198762) | ProstT5 3Di aligned with MAFFT G-INS-i `--aamatrix` (Foldseek 3Di matrix), for 3Di substitution models | no SP/TC against references; AA and 3Di aligned separately |
| Puente-Lelievre et al. 2024 (cited by Garg & Hochberg for partitioned AA+3Di models; not read) | 3Di for tree inference | not checked |
| Muscle-3D (Edgar & Tolstoy, bioRxiv 10.1101/2024.10.26.620413) | structure alphabet (Reseek "Mega") from real structures | benchmarked on structures, not sequence-only |
| Santus et al., nf-core/multiplesequencealign (NAR Genom Bioinf 2025, lqaf104; PMC12311786) | framework includes FoldMason, mTM-align, 3D-Coffee with structures | no ProstT5/3Di mention |
| learnMSA2 (Becker & Stanke, Bioinformatics 2024) | ProtT5 embeddings inside a pHMM (not 3Di); ~+6 TC points on HomFam, needs many sequences | different route (embeddings), large families only |
| PROMALS / MSACompro (older) | predicted secondary structure / contacts as alignment evidence | the classic analogue; predicted 3Di is the modern version |

Conclusion: the tool chain (ProstT5 -> FoldMason, ProstT5 -> MAFFT with a 3Di matrix) exists and is used for phylogenomics,
but I found **no published measurement of its MSA accuracy on BAliBASE/HOMSTRAD**, no paired comparison with MAFFT L-INS-i,
and nothing on predicted-3Di evidence inside divide-and-conquer aligners (MAGUS/GCM). The idea looks open; the risk is that
the FoldMason Science supplement or a 2026 preprint I could not see contains such a table.

## Methods

**Data.** BAliBASE 3.0 (downloaded from lbgi.fr): RV11 (<20% identity; 38 full-length "BB" + 38 truncated "BBS" sets) and
RV12 (20-40%; 44 BB + 44 BBS). 164 sets, 4-30 sequences each.
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

Mean SP / TC per group; Δ = paired difference to L-INS-i; W/T/L = sets where the method is better / tied / worse than
L-INS-i on SP; p = Wilcoxon signed-rank on SP. Full output (also TC tests, per-subgroup BB vs BBS): `results/summary_vs_linsi.txt`;
per-set scores: `results/scores.csv`.

**RV11 (<20% identity, n = 76)**

| method | SP | TC | ΔSP | W/T/L | p | ΔTC |
|---|---|---|---|---|---|---|
| L-INS-i (baseline) | 0.678 | 0.463 | | | | |
| Clustal Omega | 0.608 | 0.378 | −0.070 | 19/1/56 | 8e-8 | −0.085 |
| MAFFT default | 0.577 | 0.312 | −0.101 | 9/1/66 | 4e-9 | −0.151 |
| FoldMason on predicted 3Di (`fm`) | 0.703 | 0.473 | +0.024 | 41/1/34 | 0.21 | +0.010 |
| FoldMason, 3Di off (`fm_aa`) | 0.417 | 0.121 | −0.261 | 9/1/66 | 2e-12 | −0.342 |
| **L-INS-i on predicted 3Di (`linsi3di`)** | **0.727** | 0.504 | **+0.048** | 48/3/25 | **0.002** | +0.041 |
| M-Coffee linsi+linsi3di (`mc2`) | 0.708 | 0.499 | +0.030 | **62/9/5** | 2e-9 | +0.036 |
| M-Coffee linsi+linsi3di+fm (`mc3di`) | 0.725 | **0.513** | +0.047 | 53/4/19 | 9e-6 | **+0.050** |
| M-Coffee of 3 AA aligners (`mcaa`, control) | 0.650 | 0.410 | −0.028 | 22/3/51 | 3e-4 | −0.053 |

By subgroup, `linsi3di` vs L-INS-i: RV11 full-length +0.057 (p = 0.02), truncated +0.040 (p = 0.05).

**RV12 (20-40% identity, n = 88)**

| method | SP | TC | ΔSP | W/T/L | p | ΔTC |
|---|---|---|---|---|---|---|
| L-INS-i (baseline) | 0.939 | 0.849 | | | | |
| Clustal Omega | 0.905 | 0.789 | −0.034 | 20/1/67 | 3e-8 | −0.060 |
| FoldMason on predicted 3Di | 0.900 | 0.776 | −0.039 | 11/0/77 | 2e-13 | −0.073 |
| L-INS-i on predicted 3Di | 0.912 | 0.791 | −0.027 | 24/0/64 | 4e-8 | −0.058 |
| M-Coffee linsi+linsi3di (`mc2`) | 0.939 | 0.850 | −0.000 | 52/1/35 | 0.28 | +0.001 |
| M-Coffee linsi+linsi3di+fm | 0.930 | 0.827 | −0.009 | 22/4/62 | 7e-6 | −0.022 |

**Predicted vs experimental 3Di (oracle), RV11 sets where every sequence maps to a PDB chain (n = 73), paired:**
L-INS-i on true 3Di is +0.012 SP over L-INS-i on predicted 3Di (30/2/41, p = 0.23); FoldMason on true 3Di is +0.008 over
FoldMason on predicted 3Di (p = 0.33). ProstT5's 3Di is about as useful as the real thing here; the ceiling is the
alphabet/aligner, not the prediction. (The RV12 oracle is not usable: most RV12 names do not map to downloadable PDB chains,
`results/true3di_coverage.tsv`.)

**A reference-free gate.** Mean pairwise identity of the L-INS-i alignment (computed without the reference) separates the
two regimes (`results/identity_gain.png`, `results/identity_switch.txt`): `linsi3di` − L-INS-i is +0.056 SP for identity
15-20% (n = 29), +0.035 for 20-25% (n = 45), −0.004 for 25-30% (n = 17), −0.025 for 30-40% (n = 70). Rule "use the 3Di
alignment if identity < 0.25, else L-INS-i" over all 164 sets: SP 0.838 vs 0.818, +0.020, W/T/L 47/93/24, p = 0.002; with
`mc3di` as the low-identity choice: +0.022, 52/94/18, p = 6e-6. **Caveat: the threshold was picked on the same data**
(the plateau 0.25-0.30 is broad, and on BAliBASE the identity estimate mostly recovers the RV11/RV12 label), so this needs a
held-out test (HOMSTRAD, BAliFam, QuanTest) before it is a claim.

Other observations:
- FoldMason's engine is the weak part, not 3Di: with 3Di switched off it collapses (0.417 on RV11), and with 3Di it is below
  L-INS-i-on-3Di everywhere (paired, RV11: −0.024, p = 0.001). MAFFT's iterative consistency on the 3Di string does better
  than FoldMason's progressive AA+3Di profile alignment without its lDDT refinement.
- Combination is not just "ensembling": M-Coffee of three AA aligners is *worse* than L-INS-i (−0.028 on RV11), while
  M-Coffee with one 3Di alignment is better on 62 of 76 RV11 sets and loses on only 5.
- 3Di gap-open sensitivity (`results/gap_sensitivity.txt`): `--op` 1.0 / 1.53 (default) / 2.5 changes RV11-BB SP by ≤ 0.003.
- Failure cases exist: on 3 RV11 sets (BB11006, BB11017, BB11035, both versions) `linsi3di` loses 0.12-0.74 SP; `mc2` stays
  within 0.04 of L-INS-i on all six (better on two), which is why the consensus is the safer default.

## Results: large RV100 sets and MAGUS evidence

RV100_TABLE

## Runtime

RUNTIME

## Verdict

VERDICT
