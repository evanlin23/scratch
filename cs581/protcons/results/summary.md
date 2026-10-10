## Paired Δ vs MAGUS's own merge (error points; negative = better)

### Pooled

| group | method | n | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | Δ TC |
|---|---|---|---|---|---|---|---|---|
| **all held-out protein** | primary: L-INS-i \| cons0.7 | 27 | -0.60 | 14/5/8 | 0.158 | +0.92 | -2.12 | -0.19 |
| **all held-out protein** | L-INS-i ∩ FFT-NS-2 --op 3 | 27 | -1.23 | 22/1/4 | 0.000479 | -0.75 | -1.72 | +1.06 |
| **all held-out protein** | (L-INS-i + Clustal) \| cons0.7 | 27 | +0.94 | 11/3/13 | 0.454 | +5.18 | -3.30 | +2.07 |
| simulated (SIMMOD+SIMHIGH) | primary: L-INS-i \| cons0.7 | 8 | -0.85 | 5/0/3 | 0.109 | -1.31 | -0.38 | +1.35 |
| simulated (SIMMOD+SIMHIGH) | L-INS-i ∩ FFT-NS-2 --op 3 | 8 | -2.28 | 8/0/0 | 0.00781 | -4.13 | -0.42 | +4.17 |
| simulated (SIMMOD+SIMHIGH) | (L-INS-i + Clustal) \| cons0.7 | 8 | +4.52 | 2/0/6 | 0.0781 | +9.68 | -0.65 | +10.55 |
| protein excl. BAliBASE | primary: L-INS-i \| cons0.7 | 20 | -0.01 | 8/4/8 | 0.913 | +0.67 | -0.68 | -0.26 |
| protein excl. BAliBASE | L-INS-i ∩ FFT-NS-2 --op 3 | 20 | -0.97 | 15/1/4 | 0.00701 | -1.21 | -0.74 | +1.46 |
| protein excl. BAliBASE | (L-INS-i + Clustal) \| cons0.7 | 20 | +2.48 | 4/3/13 | 0.0126 | +6.12 | -1.17 | +2.77 |

### Per dataset family

| group | method | n | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | Δ TC |
|---|---|---|---|---|---|---|---|---|
| SIMMOD | primary: L-INS-i \| cons0.7 | 4 | -0.79 | 3/0/1 | 0.25 | -1.56 | -0.03 | +1.92 |
| SIMMOD | L-INS-i ∩ FFT-NS-2 --op 3 | 4 | -2.07 | 4/0/0 | 0.125 | -3.91 | -0.24 | +5.50 |
| SIMMOD | (L-INS-i + Clustal) \| cons0.7 | 4 | +0.10 | 2/0/2 | 1 | -0.31 | +0.51 | +14.71 |
| SIMHIGH | primary: L-INS-i \| cons0.7 | 4 | -0.90 | 2/0/2 | 0.625 | -1.07 | -0.73 | +0.78 |
| SIMHIGH | L-INS-i ∩ FFT-NS-2 --op 3 | 4 | -2.48 | 4/0/0 | 0.125 | -4.35 | -0.61 | +2.85 |
| SIMHIGH | (L-INS-i + Clustal) \| cons0.7 | 4 | +8.93 | 0/0/4 | 0.125 | +19.66 | -1.81 | +6.38 |
| 10AA | primary: L-INS-i \| cons0.7 | 2 | +0.02 | 0/2/0 | 0.5 | +0.03 | +0.01 | +0.00 |
| 10AA | L-INS-i ∩ FFT-NS-2 --op 3 | 2 | +0.10 | 1/0/1 | 1 | +0.35 | -0.16 | +0.00 |
| 10AA | (L-INS-i + Clustal) \| cons0.7 | 2 | +0.02 | 1/0/1 | 1 | +0.42 | -0.39 | +0.00 |
| HomFam | primary: L-INS-i \| cons0.7 | 10 | +0.66 | 3/2/5 | 0.25 | +2.38 | -1.06 | -1.60 |
| HomFam | L-INS-i ∩ FFT-NS-2 --op 3 | 10 | -0.15 | 6/1/3 | 0.426 | +0.82 | -1.11 | -0.42 |
| HomFam | (L-INS-i + Clustal) \| cons0.7 | 10 | +1.34 | 1/3/6 | 0.0547 | +4.42 | -1.74 | -2.89 |
| BAliBASE (fresh draws) | primary: L-INS-i \| cons0.7 | 7 | -2.29 | 6/1/0 | 0.0312 | +1.66 | -6.25 | +0.01 |
| BAliBASE (fresh draws) | L-INS-i ∩ FFT-NS-2 --op 3 | 7 | -1.98 | 7/0/0 | 0.0156 | +0.56 | -4.53 | -0.07 |
| BAliBASE (fresh draws) | (L-INS-i + Clustal) \| cons0.7 | 7 | -3.46 | 7/0/0 | 0.0156 | +2.49 | -9.41 | +0.08 |
| nucleotide (held out) | primary: L-INS-i \| cons0.7 | 5 | +4.66 | 0/1/4 | 0.125 | +10.73 | -1.40 | -0.00 |
| nucleotide (held out) | L-INS-i ∩ FFT-NS-2 --op 3 | 5 | +4.76 | 0/0/5 | 0.0625 | +9.18 | +0.34 | -0.67 |
| nucleotide (held out) | (L-INS-i + Clustal) \| cons0.7 | 5 | +24.17 | 0/0/5 | 0.0625 | +51.11 | -2.78 | +0.00 |
| 16S.M (in-sample) | primary: L-INS-i \| cons0.7 | 1 | -0.37 | 1/0/0 | nan | +0.18 | -0.92 | +0.16 |
| 16S.M (in-sample) | L-INS-i ∩ FFT-NS-2 --op 3 | 1 | -0.40 | 1/0/0 | nan | +0.92 | -1.73 | +0.16 |
| 16S.M (in-sample) | (L-INS-i + Clustal) \| cons0.7 | 1 | +3.67 | 0/0/1 | nan | +10.28 | -2.95 | +0.69 |

### Per dataset

| dataset | family | MAGUS err | merge ctrl err | Δ L\|cons | Δ L∩F | Δ (L+C)\|cons | ΔSPFN / ΔSPFP (L\|cons) | support | gate | frac L-only |
|---|---|---|---|---|---|---|---|---|---|---|
| SIMMOD_R1 | SIMMOD | 13.13 | 13.13 | +0.14 | -0.97 | +1.51 | +0.10 / +0.17 | 0.8247 | keep | 0.391 |
| SIMMOD_R2 | SIMMOD | 11.27 | 11.27 | -1.63 | -3.79 | -2.75 | -3.02 / -0.24 | 0.8743 | keep | 0.3425 |
| SIMMOD_R3 | SIMMOD | 11.16 | 11.16 | -1.32 | -2.51 | -0.52 | -2.59 / -0.06 | 0.8422 | keep | 0.3495 |
| SIMMOD_R4 | SIMMOD | 10.26 | 10.26 | -0.36 | -1.01 | +2.16 | -0.71 / -0.00 | 0.8338 | keep | 0.3875 |
| SIMHIGH_R1 | SIMHIGH | 25.92 | 26.01 | +0.27 | -0.89 | +11.71 | +1.04 / -0.50 | 0.6576 | keep | 0.5811 |
| SIMHIGH_R2 | SIMHIGH | 22.83 | 22.83 | -1.99 | -4.17 | +3.57 | -3.68 / -0.31 | 0.7565 | keep | 0.5112 |
| SIMHIGH_R3 | SIMHIGH | 25.73 | 25.73 | -2.05 | -3.25 | +8.34 | -3.05 / -1.05 | 0.6833 | keep | 0.5372 |
| SIMHIGH_R4 | SIMHIGH | 22.70 | 22.70 | +0.16 | -1.60 | +12.08 | +1.41 / -1.08 | 0.6683 | keep | 0.5539 |
| 10AA_1GADBL | 10AA | 3.21 | 3.21 | +0.01 | -0.14 | -0.08 | +0.02 / +0.00 | 0.8753 | keep | 0.042 |
| 10AA_coliepi | 10AA | 4.09 | 4.09 | +0.03 | +0.33 | +0.11 | +0.04 / +0.02 | 0.6962 | keep | 0.0433 |
| HF_Acetyltransf | HomFam | 29.04 | 29.04 | -2.31 | -1.60 | -1.86 | -0.05 / -4.56 | 0.5641 | filter | 0.2626 |
| HF_PDZ | HomFam | 14.64 | 14.64 | +3.54 | -0.65 | +4.11 | +6.76 / +0.32 | 0.5368 | filter | 0.299 |
| HF_aat | HomFam | 12.60 | 12.60 | +0.79 | -1.00 | +1.16 | +1.69 / -0.10 | 0.5953 | filter | 0.2869 |
| HF_adh | HomFam | 1.03 | 1.03 | +0.00 | +0.00 | +0.00 | +0.00 / +0.00 | 0.5882 | filter | 0.2702 |
| HF_blmb | HomFam | 20.31 | 20.31 | +3.42 | +2.09 | +7.69 | +8.90 / -2.05 | 0.4716 | filter | 0.4579 |
| HF_p450 | HomFam | 20.31 | 20.31 | +1.10 | +0.43 | +1.87 | +6.14 / -3.94 | 0.5059 | filter | 0.3168 |
| HF_rrm | HomFam | 19.76 | 19.76 | -0.09 | -0.33 | +0.01 | +0.05 / -0.22 | 0.7266 | keep | 0.1375 |
| HF_rvp | HomFam | 17.80 | 17.80 | +0.00 | -0.09 | +0.15 | +0.00 / +0.00 | 0.6393 | keep | 0.0131 |
| HF_sdr | HomFam | 23.17 | 23.17 | +0.46 | +0.13 | +0.05 | +0.88 / +0.04 | 0.7289 | keep | 0.1813 |
| HF_zf-CCHH | HomFam | 12.83 | 12.83 | -0.33 | -0.44 | +0.23 | -0.61 / -0.05 | 0.8 | keep | 0.0369 |
| BB_BBA0039 | BAliBASE (fresh draws) | 4.40 | 4.40 | +0.01 | -0.10 | -0.31 | +0.02 / -0.01 | 0.6095 | filter | 0.0447 |
| BB_BBA0067 | BAliBASE (fresh draws) | 25.97 | 25.97 | -1.87 | -0.38 | -2.17 | +2.22 / -5.95 | 0.4744 | filter | 0.2938 |
| BB_BBA0081 | BAliBASE (fresh draws) | 60.35 | 60.64 | -9.78 | -9.60 | -14.87 | +3.62 / -23.18 | 0.2059 | filter | 0.7224 |
| BB_BBA0101 | BAliBASE (fresh draws) | 28.71 | 28.71 | -2.49 | -1.24 | -3.50 | +1.91 / -6.89 | 0.431 | filter | 0.3006 |
| BB_BBA0117 | BAliBASE (fresh draws) | 13.15 | 13.15 | -0.34 | -0.44 | -0.68 | -0.71 / +0.02 | 0.6081 | filter | 0.0892 |
| BB_BBA0134 | BAliBASE (fresh draws) | 18.61 | 18.58 | -0.60 | -0.81 | -0.44 | +3.21 / -4.40 | 0.3643 | filter | 0.2818 |
| BB_BBA0154 | BAliBASE (fresh draws) | 21.40 | 21.40 | -0.99 | -1.29 | -2.26 | +1.34 / -3.31 | 0.4347 | filter | 0.1634 |
| 1000L1_R0 | nucleotide (held out) | 7.47 | 7.47 | +9.16 | +7.22 | +41.76 | +21.42 / -3.10 | 0.6607 | keep | 0.893 |
| 1000M2_R1 | nucleotide (held out) | 11.01 | 11.01 | +10.25 | +11.76 | +38.53 | +23.24 / -2.74 | 0.6302 | keep | 0.8349 |
| 1000M3_R0 | nucleotide (held out) | 4.52 | 4.52 | +0.16 | +0.42 | +0.79 | +0.27 / +0.04 | 0.7753 | keep | 0.1898 |
| 1000S1_R0 | nucleotide (held out) | 9.78 | 9.78 | +3.77 | +3.89 | +39.61 | +8.70 / -1.15 | 0.7061 | keep | 0.8067 |
| RNASim_R1 | nucleotide (held out) | 8.93 | 8.93 | -0.02 | +0.52 | +0.15 | +0.03 / -0.06 | 0.7751 | keep | 0.1884 |
| 16S.M_R0 | 16S.M (in-sample) | 12.86 | 12.86 | -0.37 | -0.40 | +3.67 | +0.18 / -0.92 | 0.8281 | keep | 0.2358 |

## Gate (support < 0.615 → filter) vs observed effect of the primary method

Observed help = Δ(L\|cons0.7) < 0 (strict sign).

- all: n = 33; accuracy 18/33 = 55% (filter & helped 7, filter & hurt 6, keep & would-hurt 11, keep & would-help 9)
- protein: n = 27; accuracy 14/27 = 52% (filter & helped 7, filter & hurt 6, keep & would-hurt 7, keep & would-help 7)
- nucleotide (held out + 16S): n = 6; accuracy 4/6 = 67% (filter & helped 0, filter & hurt 0, keep & would-hurt 4, keep & would-help 2)

Gated policy (primary method if the gate says filter, else MAGUS's own merge), Δ vs MAGUS:

| group | n | mean Δ | W/T/L | p | n filtered | mean Δ if always filtering |
|---|---|---|---|---|---|---|
| **all held-out protein** | 27 | -0.35 | 7/16/4 | 0.583 | 13 | -0.60 |
| simulated (SIMMOD+SIMHIGH) | 8 | +0.00 | 0/8/0 | nan | 0 | -0.85 |
| protein excl. BAliBASE | 20 | +0.33 | 1/15/4 | 0.225 | 6 | -0.01 |
| SIMMOD | 4 | +0.00 | 0/4/0 | nan | 0 | -0.79 |
| SIMHIGH | 4 | +0.00 | 0/4/0 | nan | 0 | -0.90 |
| 10AA | 2 | +0.00 | 0/2/0 | nan | 0 | +0.02 |
| HomFam | 10 | +0.65 | 1/5/4 | 0.312 | 6 | +0.66 |
| BAliBASE (fresh draws) | 7 | -2.29 | 6/1/0 | 0.0312 | 7 | -2.29 |
| nucleotide (held out) | 5 | +0.00 | 0/5/0 | nan | 0 | +4.66 |
| 16S.M (in-sample) | 1 | +0.00 | 0/1/0 | nan | 0 | -0.37 |

## POST-HOC (not pre-registered): global edge-support threshold

### Edge support (GCM edges kept iff >= k backbones contribute)

| group | method | n | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | Δ TC |
|---|---|---|---|---|---|---|---|---|
| **all held-out protein** | edge support ≥ 2 | 22 | -1.34 | 11/5/6 | 0.0496 | -2.53 | -0.14 | +1.39 |
| **all held-out protein** | edge support ≥ 3 | 22 | -1.83 | 11/4/7 | 0.063 | -3.44 | -0.23 | +3.21 |
| **all held-out protein** | edge support ≥ 5 | 22 | -2.16 | 11/4/7 | 0.0386 | -3.71 | -0.61 | +4.96 |
| simulated (SIMMOD+SIMHIGH) | edge support ≥ 2 | 8 | -3.73 | 8/0/0 | 0.00781 | -7.67 | +0.22 | +3.89 |
| simulated (SIMMOD+SIMHIGH) | edge support ≥ 3 | 8 | -5.18 | 8/0/0 | 0.00781 | -10.62 | +0.25 | +8.92 |
| simulated (SIMMOD+SIMHIGH) | edge support ≥ 5 | 8 | -5.81 | 8/0/0 | 0.00781 | -11.60 | -0.03 | +13.76 |
| protein excl. BAliBASE | edge support ≥ 2 | 20 | -1.48 | 11/4/5 | 0.0329 | -2.77 | -0.19 | +1.53 |
| protein excl. BAliBASE | edge support ≥ 3 | 20 | -2.04 | 11/3/6 | 0.0364 | -3.77 | -0.30 | +3.53 |
| protein excl. BAliBASE | edge support ≥ 5 | 20 | -2.38 | 11/3/6 | 0.0269 | -4.06 | -0.70 | +5.46 |
| SIMMOD | edge support ≥ 2 | 4 | -2.74 | 4/0/0 | 0.125 | -5.61 | +0.12 | +5.16 |
| SIMMOD | edge support ≥ 3 | 4 | -3.87 | 4/0/0 | 0.125 | -7.90 | +0.15 | +11.32 |
| SIMMOD | edge support ≥ 5 | 4 | -4.24 | 4/0/0 | 0.125 | -8.50 | +0.03 | +18.00 |
| SIMHIGH | edge support ≥ 2 | 4 | -4.71 | 4/0/0 | 0.125 | -9.74 | +0.31 | +2.62 |
| SIMHIGH | edge support ≥ 3 | 4 | -6.50 | 4/0/0 | 0.125 | -13.34 | +0.35 | +6.51 |
| SIMHIGH | edge support ≥ 5 | 4 | -7.39 | 4/0/0 | 0.125 | -14.70 | -0.08 | +9.53 |
| 10AA | edge support ≥ 2 | 2 | -0.00 | 0/2/0 | 0.5 | -0.01 | +0.00 | +0.00 |
| 10AA | edge support ≥ 3 | 2 | -0.00 | 0/2/0 | 1 | -0.01 | +0.00 | +0.00 |
| 10AA | edge support ≥ 5 | 2 | +0.01 | 0/2/0 | 0.5 | -0.01 | +0.02 | +0.24 |
| HomFam | edge support ≥ 2 | 10 | +0.02 | 3/2/5 | 0.426 | +0.60 | -0.56 | -0.06 |
| HomFam | edge support ≥ 3 | 10 | +0.08 | 3/1/6 | 0.426 | +0.97 | -0.81 | -0.08 |
| HomFam | edge support ≥ 5 | 10 | -0.11 | 3/1/6 | 0.82 | +1.17 | -1.38 | -0.15 |
| BAliBASE (fresh draws) | edge support ≥ 2 | 2 | +0.13 | 0/1/1 | 1 | -0.11 | +0.37 | +0.00 |
| BAliBASE (fresh draws) | edge support ≥ 3 | 2 | +0.19 | 0/1/1 | 0.5 | -0.12 | +0.49 | +0.00 |
| BAliBASE (fresh draws) | edge support ≥ 5 | 2 | +0.04 | 0/1/1 | 1 | -0.21 | +0.28 | +0.00 |
| nucleotide (held out) | edge support ≥ 2 | 0 | | | |   |  |  |
| nucleotide (held out) | edge support ≥ 3 | 0 | | | |   |  |  |
| nucleotide (held out) | edge support ≥ 5 | 0 | | | |   |  |  |
| 16S.M (in-sample) | edge support ≥ 2 | 0 | | | |   |  |  |
| 16S.M (in-sample) | edge support ≥ 3 | 0 | | | |   |  |  |
| 16S.M (in-sample) | edge support ≥ 5 | 0 | | | |   |  |  |

| dataset | edge support ≥ 2 | edge support ≥ 3 | edge support ≥ 5 |
|---|---|---|---|
| SIMMOD_R1 | -2.32 | -3.36 | -3.76 |
| SIMMOD_R2 | -3.73 | -5.26 | -5.55 |
| SIMMOD_R3 | -2.78 | -3.88 | -4.25 |
| SIMMOD_R4 | -2.13 | -3.00 | -3.40 |
| SIMHIGH_R1 | -4.92 | -6.41 | -7.23 |
| SIMHIGH_R2 | -6.15 | -8.63 | -9.31 |
| SIMHIGH_R3 | -3.76 | -5.84 | -7.27 |
| SIMHIGH_R4 | -4.02 | -5.10 | -5.75 |
| 10AA_1GADBL | -0.01 | -0.01 | +0.01 |
| 10AA_coliepi | -0.00 | +0.00 | +0.01 |
| HF_Acetyltransf | -3.86 | -2.86 | -3.58 |
| HF_PDZ | -0.38 | -0.65 | -0.72 |
| HF_aat | +0.02 | -0.11 | -0.80 |
| HF_adh | +0.00 | +0.00 | +0.00 |
| HF_blmb | +1.98 | +1.05 | +0.27 |
| HF_p450 | +0.47 | +0.88 | +1.42 |
| HF_rrm | +0.17 | +0.37 | +0.28 |
| HF_rvp | +0.40 | +0.40 | +0.40 |
| HF_sdr | +1.53 | +1.51 | +1.48 |
| HF_zf-CCHH | -0.16 | +0.20 | +0.20 |
| BB_BBA0039 | -0.02 | +0.01 | -0.01 |
| BB_BBA0067 | +0.28 | +0.36 | +0.08 |
| BB_BBA0081 |  |  |  |
| BB_BBA0101 |  |  |  |
| BB_BBA0117 |  |  |  |
| BB_BBA0134 |  |  |  |
| BB_BBA0154 |  |  |  |
| 1000L1_R0 |  |  |  |
| 1000M2_R1 |  |  |  |
| 1000M3_R0 |  |  |  |
| 1000S1_R0 |  |  |  |
| RNASim_R1 |  |  |  |
| 16S.M_R0 |  |  |  |

## Tree error (FastTree -lg -gamma, nRF to the true tree, %)

| dataset | true | magus | cons | fftint | lccons |
|---|---|---|---|---|---|
| SIMMOD_R1 | 6.52 | 5.72 | 6.52 | 5.82 | 6.12 |
| SIMMOD_R2 | 6.52 | 6.92 | 6.92 | 6.92 | 7.22 |
| SIMMOD_R3 | 6.12 | 7.42 | 7.62 | 6.82 | 7.42 |
| SIMMOD_R4 | 5.52 | 6.82 | 6.52 | 6.82 | 7.12 |
| SIMHIGH_R1 | 6.12 | 11.23 | 13.54 | 10.63 | 12.34 |
| SIMHIGH_R2 | 5.22 | 10.33 | 10.53 | 9.43 | 8.93 |
| SIMHIGH_R3 | 6.42 | 8.53 | 8.43 | 9.13 | 8.53 |
| SIMHIGH_R4 | 6.62 | 9.03 | 8.83 | 9.53 | 9.83 |

Δ nRF cons − MAGUS: n = 8, mean +0.36, W/T/L 3/1/4, p = 0.469

Δ nRF fftint − MAGUS: n = 8, mean -0.11, W/T/L 3/2/3, p = 0.625

Δ nRF lccons − MAGUS: n = 8, mean +0.19, W/T/L 1/2/5, p = 0.406

## Runtime (seconds, 4 threads; one job at a time for proteins)

Extra = prep (new backbone alignments + masking/intersection) + (merge − control merge). Gate = Clustal backbones + support statistic.

| dataset | MAGUS e2e | ctrl merge | extra L\|cons | extra L∩F | extra (L+C)\|cons | gate |
|---|---|---|---|---|---|---|
| SIMMOD_R1 | 948 | 37 | -10 | -22 | 21 | 40 |
| SIMMOD_R2 | 828 | 30 | -1 | -15 | 18 | 36 |
| SIMMOD_R3 | 764 | 27 | -5 | -14 | 23 | 39 |
| SIMMOD_R4 | 837 | 27 | -1 | -8 | 35 | 42 |
| SIMHIGH_R1 | 999 | 56 | -14 | -26 | 3 | 48 |
| SIMHIGH_R2 | 1261 | 40 | 20 | 1 | 25 | 49 |
| SIMHIGH_R3 | 1147 | 61 | -31 | -30 | 16 | 52 |
| SIMHIGH_R4 | 1464 | 45 | -11 | -12 | 40 | 66 |
| 10AA_1GADBL | 229 | 2 | 7 | 4 | 36 | 48 |
| 10AA_coliepi | 34 | 1 | 3 | 1 | 10 | 12 |
| HF_Acetyltransf | 341 | 6 | 2 | 2 | 6 | 9 |
| HF_PDZ | 266 | 4 | 2 | 2 | 7 | 8 |
| HF_aat | 1097 | 11 | 4 | 7 | 42 | 55 |
| HF_adh | 332 | 2 | 2 | 2 | 12 | 13 |
| HF_blmb | 734 | 21 | -5 | 8 | 18 | 33 |
| HF_p450 | 1734 | 14 | 6 | 13 | 51 | 59 |
| HF_rrm | 241 | 1 | 1 | 2 | 6 | 6 |
| HF_rvp | 264 | 1 | 2 | 2 | 8 | 9 |
| HF_sdr | 450 | 2 | 3 | 3 | 17 | 19 |
| HF_zf-CCHH | 203 | 1 | 0 | 2 | 1 | 1 |
| BB_BBA0039 | 578 | 3 | 8 | 6 | 47 | 59 |
| BB_BBA0067 | 2473 | 42 | -26 | 11 | 79 | 118 |
| BB_BBA0081 | 3770 | 272 | -237 | -75 | -107 | 131 |
| BB_BBA0101 | 1890 | 42 | -24 | 12 | 55 | 95 |
| BB_BBA0117 | 58 | 1 | 1 | 2 | 5 | 5 |
| BB_BBA0134 | 1609 | 115 | -101 | -42 | -41 | 68 |
| BB_BBA0154 | 965 | 39 | -20 | -9 | 43 | 74 |
| 1000L1_R0 | 1870 | 30 | 3 | 67 | 675 | 744 |
| 1000M2_R1 | 1562 | 26 | 8 | 84 | 675 | 738 |
| 1000M3_R0 | 1145 | 8 | 22 | 53 | 171 | 210 |
| 1000S1_R0 | 1563 | 24 | 12 | 64 | 529 | 585 |
| RNASim_R1 | 2286 | 18 | 32 | 74 | 519 | 594 |
| 16S.M_R0 | 1741 | 13 | 16 | 30 | 181 | 231 |

Nucleotide MAGUS times are from the cached worker runs (contended); not comparable.
