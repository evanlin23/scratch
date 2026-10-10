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
| Puente-Lelievre et al. 2024 (structural phylogenetics, partitioned AA+3Di) | 3Di from (predicted) structures for tree inference | no MSA reference benchmark |
| Muscle-3D (Edgar & Tolstoy, bioRxiv 10.1101/2024.10.26.620413) | structure alphabet (Reseek "Mega") from real structures | benchmarked on structures, not sequence-only |
| Santus et al., nf-core/multiplesequencealign (NAR Genom Bioinf 2025, lqaf104; PMC12311786) | framework includes FoldMason, mTM-align, 3D-Coffee with structures | no ProstT5/3Di mention |
| learnMSA2 (Becker & Stanke, Bioinformatics 2024) | ProtT5 embeddings inside a pHMM (not 3Di); ~+6 TC points on HomFam, needs many sequences | different route (embeddings), large families only |
| PROMALS / MSACompro (older) | predicted secondary structure / contacts as alignment evidence | the classic analogue; predicted 3Di is the modern version |

Conclusion: the tool chain (ProstT5 -> FoldMason, ProstT5 -> MAFFT with a 3Di matrix) exists and is used for phylogenomics,
but I found **no published measurement of its MSA accuracy on BAliBASE/HOMSTRAD**, no paired comparison with MAFFT L-INS-i,
and nothing on predicted-3Di evidence inside divide-and-conquer aligners (MAGUS/GCM). The idea looks open; the risk is that
the FoldMason Science supplement or a 2026 preprint I could not see contains such a table.
