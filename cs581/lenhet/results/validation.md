# Harness validation: WITCH-NG Table 2, ROSE 1000M3-HF

Target: Liu & Warnow 2023, WITCH-NG, doi:10.1093/bioadv/vbad024, Table 2, 1000M3-HF.
- Scored over all sequences, with the same backbone for every method.
- The backbone is the full-length sequences. It is most likely MAGUS-estimated (the
  method's default); the paper does not say.
- Mean of 20 replicates.

Our runs (`code/validate.sh`):
- Data: ROSE-HF (doi:10.13012/B2IDB-6128941_V1) 1000M3, replicates R0–R1.
- Backbone: the 500 full-length sequences; queries are the 500 fragments.
- Scored: FastSP against `true_align_fragged.txt`.

| method | backbone | n | SPFN | SPFP | published SPFN | published SPFP |
|---|---|---|---|---|---|---|
| UPP | MAGUS + FastTree | 2 | 0.055 | 0.044 | 0.054 | 0.038 |
| WITCH | MAGUS + FastTree | 2 | 0.052 | 0.044 | 0.048 | 0.039 |
| UPP | true | 2 | 0.012 | 0.002 | – | – |
| WITCH | true | 2 | 0.010 | 0.002 | – | – |

Numbers are unmasked (upper-cased). Masking lower-case insertion letters changes SPFP by
at most 0.001 here.

- **SPFN** reproduces within 0.004.
- **SPFP** is about 0.005 above the published value. Plausible causes: only 2 of 20
  replicates, MAGUS randomness, or a different MAGUS version or settings for the
  backbone.
- **WITCH vs UPP:** the published ordering (WITCH ≤ UPP) holds.
- **Fragments only, with the true backbone:** UPP SPFN 0.057 and WITCH 0.046.
- **Most of the all-sequence error comes from the backbone.** Restricted to backbone
  sequences, the MAGUS backbone has an error of about 0.05.
