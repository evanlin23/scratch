## Protein (BAliBASE)

#### Single-tool backbones (10, MAGUS's sequence sets)

| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |
|---|---|---|---|---|---|---|---|---|---|
| clustalo | 8 | 22.90 | -1.13 | 6/0/2 | 0.109 | +0.34 | -2.60 | 62 | 75 |
| fftns2 | 6 | 20.29 | -0.30 | 3/0/3 | 0.562 | +2.10 | -2.70 | 27 | 73 |
| famsa | 6 | 22.37 | +1.77 | 0/1/5 | 0.0625 | +2.41 | +1.14 | 15 | 54 |
| clustalo-iter | 6 | 19.97 | -0.62 | 5/0/1 | 0.219 | -0.15 | -1.10 | 852 | 59 |

#### Amount and diversity of evidence

| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |
|---|---|---|---|---|---|---|---|---|---|
| linsi~5 | 8 | 23.95 | -0.08 | 3/2/3 | 0.742 | +0.08 | -0.23 | – | 69 |
| clustalo~5 | 8 | 22.82 | -1.21 | 6/0/2 | 0.109 | +0.56 | -2.97 | 31 | 92 |
| clustalo@s1 | 6 | 19.62 | -0.98 | 5/0/1 | 0.0938 | +0.20 | -2.16 | 133 | 69 |
| clustalo@n20 | 6 | 19.76 | -0.83 | 5/0/1 | 0.219 | +0.07 | -1.73 | 172 | 59 |
| clustalo@n30 | 6 | 19.71 | -0.88 | 5/0/1 | 0.219 | +0.07 | -1.84 | 307 | 65 |
| linsi+clustalo | 8 | 23.36 | -0.67 | 7/0/1 | 0.0234 | -0.39 | -0.94 | 62 | 84 |
| linsi~5+clustalo^5 | 6 | 20.11 | -0.49 | 5/0/1 | 0.156 | -0.06 | -0.92 | 19 | 51 |

#### Consensus: pair intersections (same sequence sets)

| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |
|---|---|---|---|---|---|---|---|---|---|
| linsi&clustalo | 8 | 22.01 | -2.02 | 7/0/1 | 0.0391 | +0.66 | -4.70 | 62 | 50 |
| linsi&fftns2 | 8 | 22.41 | -1.63 | 6/1/1 | 0.0234 | +1.25 | -4.50 | 24 | 59 |
| linsi&fftns2-op3 | 8 | 22.20 | -1.83 | 7/1/0 | 0.0156 | +0.51 | -4.17 | 17 | 55 |
| linsi&famsa | 6 | 19.69 | -0.91 | 4/1/1 | 0.219 | +1.66 | -3.47 | 15 | 37 |
| clustalo&fftns2 | 8 | 21.76 | -2.28 | 7/0/1 | 0.0391 | +1.93 | -6.48 | 77 | 53 |
| linsi&clustalo&fftns2 | 8 | 21.63 | -2.40 | 7/0/1 | 0.0391 | +2.02 | -6.82 | 77 | 18 |
| linsi&clustalo+clustalo@s1&fftns2@s1 | 6 | 18.81 | -1.79 | 6/0/0 | 0.0312 | +0.50 | -4.08 | 190 | 45 |
| clustalo&fftns2+clustalo@s1&fftns2@s1+clustalo@s2&fftns2@s2 | 6 | 18.34 | -2.25 | 6/0/0 | 0.0312 | +1.40 | -5.90 | 356 | 46 |

#### Consensus: column masks on MAGUS's own L-INS-i backbones (MAFFT only, no new alignments)

| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |
|---|---|---|---|---|---|---|---|---|---|
| linsi\|cons0.3 | 6 | 20.10 | -0.49 | 4/2/0 | 0.0625 | -0.05 | -0.94 | – | 26 |
| linsi\|cons0.5 | 8 | 23.19 | -0.84 | 7/0/1 | 0.0234 | +0.14 | -1.82 | – | 21 |
| linsi\|cons0.7 | 8 | 22.30 | -1.73 | 7/1/0 | 0.00781 | +1.01 | -4.47 | – | 9 |
| linsi\|cons0.8 | 8 | 22.15 | -1.88 | 7/0/1 | 0.0234 | +2.81 | -6.57 | – | 6 |
| linsi\|agree-clustalo-0.5 | 8 | 22.90 | -1.13 | 7/0/1 | 0.0156 | +1.10 | -3.35 | – | 16 |

#### Union plus consistency mask

| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |
|---|---|---|---|---|---|---|---|---|---|
| linsi+clustalo\|cons0.5 | 8 | 22.18 | -1.86 | 8/0/0 | 0.00781 | +0.21 | -3.92 | 62 | 33 |
| linsi+clustalo\|cons0.7 | 8 | 21.51 | -2.53 | 7/0/1 | 0.0156 | +2.01 | -7.06 | 62 | 15 |
| linsi+linsi&clustalo | 8 | 23.80 | -0.24 | 6/0/2 | 0.109 | +0.04 | -0.51 | 62 | 86 |
| linsi+linsi&fftns2 | 6 | 20.35 | -0.24 | 4/1/1 | 0.0938 | +0.00 | -0.49 | 27 | 64 |
| linsi+linsi&clustalo&fftns2 | 6 | 20.26 | -0.34 | 5/0/1 | 0.0938 | -0.04 | -0.63 | 52 | 67 |

#### MAFFT settings and gappy-column masks

| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |
|---|---|---|---|---|---|---|---|---|---|
| linsi-op3@h | 2 | 24.37 | +0.60 | 0/0/2 | nan | +0.91 | +0.30 | 1201 | 108 |
| fftns2@h | 2 | 23.81 | +0.05 | 1/0/1 | nan | +3.51 | -3.42 | 8 | 122 |
| linsi@h&fftns2@h | 6 | 19.20 | -1.40 | 5/0/1 | 0.0625 | +1.45 | -4.24 | – | 28 |
| ginsi@h | 2 | 22.84 | -0.92 | 2/0/0 | nan | -0.49 | -1.36 | 1224 | 105 |
| ginsi-ul8@h | 2 | 22.45 | -1.31 | 1/1/0 | nan | +2.07 | -4.70 | 3303 | 124 |
| ginsi-ul4@h | 1 | 27.96 | -1.22 | 1/0/0 | nan | -0.28 | -2.17 | 1396 | 56 |
| linsi\|gap0.5 | 6 | 20.72 | +0.13 | 3/0/3 | 0.688 | -0.30 | +0.55 | – | 25 |
| linsi\|gap0.7 | 6 | 20.72 | +0.12 | 2/0/4 | 0.312 | -0.04 | +0.28 | – | 36 |
| linsi\|gap0.9 | 6 | 20.91 | +0.32 | 1/1/4 | 0.156 | +0.42 | +0.21 | – | 42 |
| clustalo\|gap0.7 | 6 | 19.87 | -0.72 | 5/0/1 | 0.438 | +0.46 | -1.91 | 39 | 43 |

#### Per replicate: control error (%), then Δ vs control (points)

| variant | BBA0039 | BBA0067 | BBA0081 | BBA0101 | BBA0117 | BBA0134 | BBA0154 | BBA0190 |
|---|---|---|---|---|---|---|---|---|
| linsi | 4.69 | 26.28 | 56.28 | 29.19 | 12.41 | 18.34 | 21.66 | 23.40 |
| clustalo | -0.25 | -1.56 | -3.20 | -2.83 | +0.76 | +1.50 | -2.27 | -1.21 |
| linsi+clustalo | -0.33 | -0.39 | -0.69 | -1.57 | +0.36 | -0.49 | -1.34 | -0.89 |
| linsi&clustalo | -0.28 | -1.40 | -6.19 | -2.98 | +1.00 | -0.12 | -2.83 | -3.36 |
| linsi&fftns2 | -0.02 | -1.04 | -6.05 | -1.91 | +0.59 | -0.97 | -1.69 | -1.91 |
| linsi&fftns2-op3 | -0.16 | -0.96 | -7.01 | -1.59 | +0.03 | -1.26 | -1.76 | -1.96 |
| clustalo&fftns2 | -0.33 | -1.81 | -7.19 | -3.58 | +1.64 | -0.73 | -2.93 | -3.29 |
| linsi\|cons0.5 | -0.05 | -1.76 | -0.88 | -2.32 | +0.06 | -1.21 | -0.46 | -0.07 |
| linsi\|cons0.7 | -0.04 | -2.47 | -5.52 | -2.94 | -0.10 | -1.50 | -1.05 | -0.23 |
| linsi\|cons0.8 | -0.13 | -2.60 | -6.51 | -3.09 | +0.57 | -0.61 | -1.81 | -0.87 |
| linsi+clustalo\|cons0.5 | -0.28 | -2.68 | -4.46 | -3.83 | -0.13 | -0.82 | -1.92 | -0.74 |
| linsi+clustalo\|cons0.7 | -0.39 | -3.20 | -7.63 | -4.69 | +0.35 | -0.58 | -2.52 | -1.54 |

## Nucleotide (RNASim, 16S.M, ROSE)

#### Single-tool backbones (10, MAGUS's sequence sets)

| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |
|---|---|---|---|---|---|---|---|---|---|
| clustalo | 5 | 22.76 | +13.11 | 1/0/4 | 0.125 | +20.21 | +6.02 | 837 | 231 |

#### Amount and diversity of evidence

| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |
|---|---|---|---|---|---|---|---|---|---|
| linsi+clustalo | 5 | 9.64 | -0.01 | 2/1/2 | 1 | +0.10 | -0.12 | 837 | 193 |

#### Consensus: pair intersections (same sequence sets)

| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |
|---|---|---|---|---|---|---|---|---|---|
| linsi&clustalo | 3 | 15.05 | +4.71 | 0/0/3 | 0.25 | +11.04 | -1.62 | 755 | 24 |
| linsi&fftns2 | 3 | 18.82 | +8.49 | 1/0/2 | 0.75 | +18.93 | -1.96 | 92 | 37 |
| linsi&fftns2-op3 | 5 | 15.04 | +5.39 | 1/0/4 | 0.125 | +11.14 | -0.36 | 71 | 38 |
| clustalo&fftns2 | 3 | 25.34 | +15.00 | 0/0/3 | 0.25 | +31.81 | -1.81 | 847 | 60 |
| linsi&clustalo+clustalo@s1&fftns2@s1 | 1 | 10.57 | +0.65 | 0/0/1 | nan | +4.19 | -2.90 | 1910 | 31 |
| clustalo&fftns2+clustalo@s1&fftns2@s1+clustalo@s2&fftns2@s2 | 1 | 11.30 | +1.38 | 0/0/1 | nan | +5.98 | -3.22 | 3242 | 49 |

#### Consensus: column masks on MAGUS's own L-INS-i backbones (MAFFT only, no new alignments)

| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |
|---|---|---|---|---|---|---|---|---|---|
| linsi\|cons0.3 | 1 | 9.88 | -0.04 | 0/1/0 | nan | -0.12 | +0.04 | – | 20 |
| linsi\|cons0.5 | 5 | 10.60 | +0.95 | 1/2/2 | 0.438 | +2.12 | -0.22 | – | 17 |
| linsi\|cons0.7 | 3 | 13.35 | +3.02 | 2/0/1 | 1 | +6.94 | -0.91 | – | 15 |
| linsi\|cons0.8 | 3 | 15.62 | +5.29 | 1/0/2 | 0.75 | +12.43 | -1.86 | – | 14 |

#### Union plus consistency mask

| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |
|---|---|---|---|---|---|---|---|---|---|
| linsi+clustalo\|cons0.5 | 3 | 22.22 | +11.89 | 1/0/2 | 0.5 | +24.95 | -1.17 | 755 | 16 |
| linsi+clustalo\|cons0.7 | 3 | 25.41 | +15.07 | 0/0/3 | 0.25 | +32.89 | -2.75 | 755 | 14 |

#### MAFFT settings and gappy-column masks

| variant | n | mean err | mean Δ err | W/T/L | p | Δ SPFN | Δ SPFP | bb s (wall, 10 bb) | merge s |
|---|---|---|---|---|---|---|---|---|---|

#### Per replicate: control error (%), then Δ vs control (points)

| variant | 1000L1 | 1000M2 | 1000S1 | 16S.M | RNASim |
|---|---|---|---|---|---|
| linsi | 7.47 | 8.23 | 9.78 | 12.86 | 9.92 |
| clustalo | +26.81 | +19.77 | +13.01 | +6.04 | -0.06 |
| linsi+clustalo | -0.03 | +0.21 | -0.08 | +0.55 | -0.70 |
| linsi&clustalo |  | +12.13 |  | +1.82 | +0.18 |
| linsi&fftns2 |  | +25.58 |  | -0.66 | +0.54 |
| linsi&fftns2-op3 | +7.22 | +15.65 | +3.89 | -0.40 | +0.58 |
| clustalo&fftns2 |  | +42.30 |  | +1.09 | +1.62 |
| linsi\|cons0.5 | +1.24 | +3.56 | +0.04 | -0.02 | -0.08 |
| linsi\|cons0.7 |  | +9.49 |  | -0.37 | -0.08 |
| linsi\|cons0.8 |  | +16.04 |  | -0.65 | +0.47 |
| linsi+clustalo\|cons0.5 |  | +35.36 |  | +1.08 | -0.78 |
| linsi+clustalo\|cons0.7 |  | +41.47 |  | +3.67 | +0.08 |

