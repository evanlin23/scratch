# Does consistency-filtered GCM evidence generalize? A pre-registered held-out test

CS581 (UIUC), Fall 2026. Pre-registration: `PREREG.md`, committed before any new data was run. Code: `code/`.
Raw rows: `results/*.results.jsonl`, `*.gate.json`, `*.trees.jsonl`, `magus.jsonl`. All tables below come from
`code/summarize.py` (`results/summary.md`).

## TL;DR

| method (merge-only, paired on MAGUS's own subsets and backbone sets) | held-out protein, n = 28: mean Δ err | W/T/L | p | protein excluding BAliBASE (n = 20) | held-out nucleotide (n = 5) |
|---|---|---|---|---|---|
| **primary: `linsi\|cons0.7`** | −0.59 | 15/5/8 | 0.13 | −0.01 (8/4/8, p = 0.91) | +4.66 (0/1/4) |
| `L ∩ FFT-NS-2 --op 3` | **−1.28** | **23/1/4** | **0.0003** | **−0.97** (15/1/4, p = 0.007) | +4.76 (0/0/5) |
| `(L + Clustal)\|cons0.7` | +0.86 | 12/3/13 | 0.55 | +2.48 (4/3/13, p = 0.013) | +24.17 (0/0/5) |
| *post hoc:* GCM edge support ≥ 5 of 10 backbones | −1.90 | 17/4/7 | 0.004 | −2.38 (11/3/6, p = 0.027) | −0.03 (2/2/1) |

1. **The pre-registered primary endpoint fails.** `linsi|cons0.7` does not significantly improve MAGUS on held-out
   proteins (−0.59, p = 0.13). It replicates **only on new MAGUS draws of the BAliBASE sets it was found on**
   (−2.03, 7/1/0, p = 0.016). Outside BAliBASE its mean effect is zero (−0.01). On HomFam it hurts (+0.66, 3/2/5).
   On the simulated proteins it helps 5 of 8 times, but through SPFN, not SPFP.
2. **The reference-free gate does not work.** Its accuracy at predicting whether the primary method helps is
   15/28 = 54% on proteins, i.e. chance. It said "filter" on all 8 BAliBASE draws, which is correct, but also on 6
   of 10 HomFam families, where filtering hurt on 4. It said "keep" on all 8 simulated sets, where filtering helped
   on 5. It said "keep" correctly on all 5 held-out nucleotide sets (filtering hurts there by up to +10 points), but
   "support > 0.615" is also true of every simulated protein set. So it separates nucleotides from BAliBASE, not
   "help" from "hurt". The gated policy gives −0.35 on proteins (p = 0.51).
3. **The secondary method `L ∩ FFT-NS-2 --op 3` does generalize on proteins.**
   - Simulated: 8/0/0, −2.28.
   - BAliBASE: 8/0/0, −2.05.
   - Protein excluding BAliBASE: −0.97, p = 0.007.
   - HomFam: mixed (−0.15, 6/1/3).

   It is harmful on nucleotides (+4.76, 0/0/5), so it is a protein-only recipe.
4. **`(L + Clustal)|cons0.7` is the best recipe on BAliBASE** (−3.18, 8/0/0), and the worst elsewhere:
   - simulated high divergence: +8.93, 0/0/4;
   - HomFam: +1.34;
   - nucleotides: +24.

   This confirms protbench's finding that Clustal evidence does not travel outside BAliBASE.
5. **Mechanism differs by data.** The BAliBASE story (GCM is precision-limited; filtering buys SPFP) holds only on
   BAliBASE: ΔSPFP −5.5, ΔSPFN +1.5 for the primary method.
   - HomFam is recall-limited, consistent with published MAGUS numbers. Filtering raises SPFN (+2.4) more than it
     lowers SPFP (−1.1).
   - On the simulated proteins every gain comes from **lower SPFN**. This is a GCM/MCL clustering effect, not an
     evidence-precision effect (see the post-hoc result).
6. **Trees:** no recipe changes FastTree accuracy on the simulated data. Δ nRF vs MAGUS: primary +0.36, L∩F −0.11,
   (L+C) +0.19 (all p > 0.4, n = 8).
7. **Post hoc (not pre-registered; requested after the freeze):** the simplest competitor is a global edge-support
   threshold, deleting GCM graph edges that fewer than k of the 10 backbones contribute to. It beats every
   pre-registered recipe on non-BAliBASE protein.
   - Simulated: k = 5 gives −5.81, 8/0/0.
   - HomFam: neutral (−0.11).
   - BAliBASE: small gain (−0.71, 6/1/1).
   - Nucleotides: neutral (−0.03; 16S −0.12).

   It is the only recipe here that did not hurt any data type on average. It is post hoc, chosen from 3 values of k,
   on the same data, and needs its own held-out test.

**Verdict: "consistency-filtered evidence" as pre-registered (`linsi|cons0.7` + support gate) is not promising.
It does not generalize beyond BAliBASE, and its gate is at chance.** The broader direction, removing weakly
supported evidence before MCL, is **unclear-to-promising**. `L ∩ FFT-NS-2` generalizes moderately on proteins,
and a post-hoc global edge-support threshold gives large gains on simulated proteins without hurting nucleotides.
A 4-week project should pivot to that: an edge-support threshold or reweighting inside GCM, pre-registered on
fresh data, with a mechanism study of why it lowers SPFN.

## 1. Pre-registration and deviations

`PREREG.md` (commit "protcons: pre-registration …") fixed:
- the methods, exactly as in `cs581/bbevidence/code/bbe.py`;
- the gate: `support_a_only` of L-INS-i vs Clustal Omega, threshold τ = 0.615, the midpoint of 0.59/0.64;
- the endpoints: paired Δ error vs MAGUS's own merge, Wilcoxon, W/T/L with a 0.05 tie band;
- the held-out data in priority order.

An addendum, also committed before those runs, named the nucleotide replicates.

Deviations and notes:
- **Gate implementation.** It is recomputed without the reference (`code/pc.py gate`). It equals pairdiff.py's value
  exactly on 10AA coli_epi (0.6962).
- **Post-hoc analyses.** The orchestrating session asked for these after the freeze; they are labelled as post hoc:
  - the edge-support baseline (`code/run_edgesup.py`; k = 1 reproduces MAGUS's merge exactly);
  - SPFN/SPFP per family.

  The 20-backbone L-INS-i and T-Coffee TCS baselines were **not run** (time budget).
- **10AA coli_epi_100** was used as the pipeline test before the full run started. Its rows were produced by the
  same frozen code; it is included.
- **Contended timings.** From BB_BBA0101 on, the BAliBASE draws ran concurrently with the nucleotide lane and the
  post-hoc lane to fit the time budget, so their wall times are contended. FastTree jobs (niced) overlapped
  10AA/HomFam. Simulated-set timings are clean (one job at a time, 4 threads).
- **No dataset was dropped or rerun.** One MAGUS draw per dataset throughout.

## 2. Data

| family | datasets | reference | notes |
|---|---|---|---|
| SIMMOD / SIMHIGH | R1–R4 each, 1,000 seqs | true alignment and true tree | AliSim (IQ-TREE 3.1.4) LG+G4, root length 300, indels 0.05/0.05 POW(1.7, 40). Yule tree `-rlen 0.001 {0.06,0.10} 0.8`, seeds 100r+7 / 100r+13. Regenerated; SIMMOD_R1 has 7,348 columns, identical to protbench. |
| 10AA | 1GADBL_100 (561), coli_epi_100 (320) | structural | IDB-2567453 `salma_paper_datasets.zip` |
| HomFam | aat, p450, sdr, adh, blmb, rrm, PDZ, Acetyltransf, rvp, zf-CCHH | Homstrad seeds (5–20) | 2,000-sequence subsamples, `protbench/code/make_homfam.py`, seed 1; scored on the seeds |
| BAliBASE, fresh draws | 8 RV100 sets, `cs581/data/balibase_clean` | structural | one new MAGUS run each (new decomposition and backbones; the bbevidence draws were not reused) |
| nucleotide, held out | 1000M2_R1, 1000L1_R0, 1000S1_R0, 1000M3_R0, RNASim_R1 | true | cached MAGUS inputs from the worker branches |
| nucleotide, in sample | 16S.M_R0 | CRW | the gate was fit on this replicate |

MAGUS: `gcmx.bbtool_bench` with the paper's flags (25 subsets, 10 backbones × 200, MCL, minclusters), 4 threads. The
merge-only control (`linsi`) reproduces MAGUS's own output exactly (to 1e-6) on 23 of 34 datasets, within 0.03 points on 9 more, and
within 0.09 / 0.29 points on SIMHIGH_R1 / BBA0081 (trace non-determinism; all Δs use the merge control).

## 3. Results per family

Δ = method − MAGUS's own merge, in error points ((SPFN+SPFP)/2, %); negative = better.

| family | n | primary Δ (W/T/L, p) | L∩F Δ (W/T/L, p) | (L+C)\|cons Δ (W/T/L, p) | primary ΔSPFN / ΔSPFP |
|---|---|---|---|---|---|
| SIMMOD | 4 | −0.79 (3/0/1, 0.25) | −2.07 (4/0/0, 0.13) | +0.10 (2/0/2) | −1.56 / −0.03 |
| SIMHIGH | 4 | −0.90 (2/0/2, 0.63) | −2.48 (4/0/0, 0.13) | +8.93 (0/0/4) | −1.07 / −0.73 |
| 10AA | 2 | +0.02 (0/2/0) | +0.10 (1/0/1) | +0.02 (1/0/1) | +0.03 / +0.01 |
| HomFam | 10 | +0.66 (3/2/5, 0.25) | −0.15 (6/1/3, 0.43) | +1.34 (1/3/6, 0.055) | +2.38 / −1.06 |
| BAliBASE fresh | 8 | **−2.03 (7/1/0, 0.016)** | **−2.05 (8/0/0, 0.008)** | **−3.18 (8/0/0, 0.008)** | +1.47 / −5.54 |
| nucleotide held out | 5 | +4.66 (0/1/4, 0.13) | +4.76 (0/0/5, 0.06) | +24.17 (0/0/5, 0.06) | +10.73 / −1.40 |
| 16S.M (in sample) | 1 | −0.37 | −0.40 | +3.67 | +0.18 / −0.92 |

Per dataset (Δ error; the gate's support statistic and decision):

| dataset | MAGUS err | Δ L\|cons | Δ L∩F | Δ (L+C)\|cons | support | gate |
|---|---|---|---|---|---|---|
| SIMMOD_R1 | 13.13 | +0.14 | −0.97 | +1.51 | 0.825 | keep |
| SIMMOD_R2 | 11.27 | −1.63 | −3.79 | −2.75 | 0.874 | keep |
| SIMMOD_R3 | 11.16 | −1.32 | −2.51 | −0.52 | 0.842 | keep |
| SIMMOD_R4 | 10.26 | −0.36 | −1.01 | +2.16 | 0.834 | keep |
| SIMHIGH_R1 | 25.92 | +0.27 | −0.89 | +11.71 | 0.658 | keep |
| SIMHIGH_R2 | 22.83 | −1.99 | −4.17 | +3.57 | 0.757 | keep |
| SIMHIGH_R3 | 25.73 | −2.05 | −3.25 | +8.34 | 0.683 | keep |
| SIMHIGH_R4 | 22.70 | +0.16 | −1.60 | +12.08 | 0.668 | keep |
| 10AA_1GADBL | 3.21 | +0.01 | −0.14 | −0.08 | 0.875 | keep |
| 10AA_coliepi | 4.09 | +0.03 | +0.33 | +0.11 | 0.696 | keep |
| HF_Acetyltransf | 29.04 | −2.31 | −1.60 | −1.86 | 0.564 | filter |
| HF_PDZ | 14.64 | +3.54 | −0.65 | +4.11 | 0.537 | filter |
| HF_aat | 12.60 | +0.79 | −1.00 | +1.16 | 0.595 | filter |
| HF_adh | 1.03 | 0.00 | 0.00 | 0.00 | 0.588 | filter |
| HF_blmb | 20.31 | +3.42 | +2.09 | +7.69 | 0.472 | filter |
| HF_p450 | 20.31 | +1.10 | +0.43 | +1.87 | 0.506 | filter |
| HF_rrm | 19.76 | −0.09 | −0.33 | +0.01 | 0.727 | keep |
| HF_rvp | 17.80 | 0.00 | −0.09 | +0.15 | 0.639 | keep |
| HF_sdr | 23.17 | +0.46 | +0.13 | +0.05 | 0.729 | keep |
| HF_zf-CCHH | 12.83 | −0.33 | −0.44 | +0.23 | 0.800 | keep |
| BB_BBA0039 | 4.40 | +0.01 | −0.10 | −0.31 | 0.610 | filter |
| BB_BBA0067 | 25.97 | −1.87 | −0.38 | −2.17 | 0.474 | filter |
| BB_BBA0081 | 60.35 | −9.78 | −9.60 | −14.87 | 0.206 | filter |
| BB_BBA0101 | 28.71 | −2.49 | −1.24 | −3.50 | 0.431 | filter |
| BB_BBA0117 | 13.15 | −0.34 | −0.44 | −0.68 | 0.608 | filter |
| BB_BBA0134 | 18.61 | −0.60 | −0.81 | −0.44 | 0.364 | filter |
| BB_BBA0154 | 21.40 | −0.99 | −1.29 | −2.26 | 0.435 | filter |
| BB_BBA0190 | 23.32 | −0.19 | −2.53 | −1.19 | 0.600 | filter |
| 1000L1_R0 | 7.47 | +9.16 | +7.22 | +41.76 | 0.661 | keep |
| 1000M2_R1 | 11.01 | +10.25 | +11.76 | +38.53 | 0.630 | keep |
| 1000M3_R0 | 4.52 | +0.16 | +0.42 | +0.79 | 0.775 | keep |
| 1000S1_R0 | 9.78 | +3.77 | +3.89 | +39.61 | 0.706 | keep |
| RNASim_R1 | 8.93 | −0.02 | +0.52 | +0.15 | 0.775 | keep |
| 16S.M_R0 | 12.86 | −0.37 | −0.40 | +3.67 | 0.828 | keep |

**BAliBASE replicates well.** The fresh draws reproduce the bbevidence pilot almost exactly:
- primary −2.03 here vs −1.97 in the pilot;
- L∩F −2.05 vs −1.83;
- (L+C)|cons −3.18 vs −2.53;
- BBA0081 is again the biggest gain (−9.8).

So the earlier result was not draw noise. It is a property of BAliBASE RV100 (with MAGUS's L-INS-i evidence, many
confident cross-subset errors that the backbones disagree on), and that property is not shared by the simulated or
HomFam data.

## 4. Gate accuracy

The gate predicts "filter helps" when support < 0.615. Observed help means the primary method's Δ < 0.

| data | n | correct | filter & helped | filter & hurt | keep & would have hurt | keep & would have helped |
|---|---|---|---|---|---|---|
| protein | 28 | 15 (54%) | 8 | 6 | 7 | 7 |
| nucleotide (5 held out + 16S) | 6 | 4 (67%) | 0 | 0 | 4 | 2 |
| all | 34 | 19 (56%) | 8 | 6 | 11 | 9 |

- The gate says "filter" on every BAliBASE set. It does so even on BBA0039 (0.6095) and BBA0117 (0.608), where the
  pilot had them on either side of the threshold.
- It says "keep" on every simulated set. On HomFam it is wrong in both directions.
- Support is lowest on the hardest BAliBASE sets (BBA0081 0.21), and there it does track the size of the gain. But
  values of 0.47–0.60 occur both where filtering helps a lot (BBA0067) and where it hurts a lot (HF_blmb, HF_PDZ).
- **Gated policy** (filter if the gate says so, else MAGUS): proteins −0.35 (8/16/4, p = 0.51). Excluding BAliBASE
  it is +0.33 (1/15/4). It is worse than not filtering at all outside BAliBASE.
- The one thing the gate does do is say "keep" on all 5 held-out nucleotide sets. There, every filtering recipe is
  harmful (+0.2 to +10 for the primary method, up to +42 with Clustal). A cruder rule, "is it protein?", would do
  the same.

## 5. Tree accuracy (simulated proteins)

FastTree 2 `-lg -gamma`, nRF (%) to the true tree:

| dataset | true aln | MAGUS | L\|cons | L∩F | (L+C)\|cons |
|---|---|---|---|---|---|
| SIMMOD_R1 | 6.52 | 5.72 | 6.52 | 5.82 | 6.12 |
| SIMMOD_R2 | 6.52 | 6.92 | 6.92 | 6.92 | 7.22 |
| SIMMOD_R3 | 6.12 | 7.42 | 7.62 | 6.82 | 7.42 |
| SIMMOD_R4 | 5.52 | 6.82 | 6.52 | 6.82 | 7.12 |
| SIMHIGH_R1 | 6.12 | 11.23 | 13.54 | 10.63 | 12.34 |
| SIMHIGH_R2 | 5.22 | 10.33 | 10.53 | 9.43 | 8.93 |
| SIMHIGH_R3 | 6.42 | 8.53 | 8.43 | 9.13 | 8.53 |
| SIMHIGH_R4 | 6.62 | 9.03 | 8.83 | 9.53 | 9.83 |

Δ nRF vs MAGUS:
- primary: +0.36 (3/1/4, p = 0.47);
- L∩F: −0.11 (3/2/3, p = 0.63);
- (L+C)|cons: +0.19 (1/2/5, p = 0.41).

The 1–4-point SP-error gains of L∩F on the simulated data do not translate into better trees. The alignment-to-tree
gap is about 2–5 nRF points on SIMHIGH, and none of the recipes closes it.

## 6. Runtime

Seconds, 4 threads. "Extra" = prep (new backbone alignments + masking/intersection) + (merge − control merge). The
merge gets faster with filtered evidence, so "extra" can be negative.

| | MAGUS end to end | extra L\|cons | extra L∩F | extra (L+C)\|cons | gate (Clustal + statistic) |
|---|---|---|---|---|---|
| SIMMOD (4, clean) | 764–948 | −10 to −1 | −22 to −8 | +18 to +35 | 36–42 |
| SIMHIGH (4, clean) | 999–1,464 | −31 to +20 | −30 to +1 | +3 to +40 | 48–66 |
| HomFam (10; FastTree overlapping) | 203–1,734 | 0 to +6 | +2 to +13 | +1 to +51 | 1–59 |
| BAliBASE (8; contended from BBA0101) | 58–3,770 | −237 to +10 | −75 to +12 | −107 to +130 | 5–162 |
| nucleotide (6) | (cached, contended) | +3 to +32 | +30 to +84 | +171 to +675 | 210–744 |

- **All three recipes cost under 5% of MAGUS's end-to-end time on proteins.** `linsi|cons0.7` needs no new
  alignment. `L∩F` needs 10 FFT-NS-2 runs (3–4 s wall on the simulated sets).
- **Clustal is the expensive part on nucleotides.** The Clustal backbones that the gate and `(L+C)` need take 2–11
  minutes there.
- **The post-hoc edge-support merge** costs about the same as the control merge (Python graph build). It is often faster:
  e.g. 5 s at k = 5 vs 27 s for the control on SIMMOD_R4, because the pruned graph is smaller for MCL.

Per dataset: `results/summary.md`.

## 7. Post hoc: global edge-support threshold (not pre-registered)

The orchestrating session asked for this after the freeze. `code/run_edgesup.py` builds MAGUS's alignment graph
unchanged, then deletes every edge between two subset columns that fewer than k of the 10 backbones contribute a
residue pair to. With k = 1 it reproduces MAGUS's merge exactly.

| group | k = 2 | k = 3 | k = 5 |
|---|---|---|---|
| simulated (8) | −3.73 (8/0/0) | −5.18 (8/0/0) | **−5.81 (8/0/0)**, ΔSPFN −11.6, ΔSPFP −0.03 |
| HomFam (10) | +0.02 (3/2/5) | +0.08 (3/1/6) | −0.11 (3/1/6) |
| BAliBASE fresh (8) | +0.02 | −0.09 | −0.71 (6/1/1, p = 0.023) |
| 10AA (2) | 0.00 | 0.00 | +0.01 |
| all protein (28) | −1.05 (p = 0.07) | −1.48 (p = 0.07) | **−1.90 (17/4/7, p = 0.004)** |
| nucleotide held out (5) | 0.00 | −0.01 | −0.03 (2/2/1) |
| 16S.M | −0.01 | −0.02 | −0.12 |

- **A clustering effect, not an evidence-precision effect.** On the simulated proteins, removing rarely supported
  edges *lowers SPFN* by 8–15 points, with SPFP unchanged and TC up 4–18 points. MAGUS on these data is
  recall-limited, which suggests that sparse, weakly supported edges make MCL/minclusters split true columns.
  This is a hypothesis; the mechanism was not tested.
- **The primary method's simulated-data gains run through the same channel.** `linsi|cons0.7` and `L∩F` also gain
  through SPFN on the simulated data, which suggests they partly act as an edge-sparsifier too.
- **Per-family pattern.** BAliBASE gains much more from column filtering (−2.0) than from the edge threshold (−0.7):
  there the errors are confident, multi-backbone errors that a support count does not remove. The simulated
  proteins are the opposite. HomFam responds to neither.
- **Caveats.** k was chosen from three values on these same data, and the comparison of post-hoc and pre-registered
  methods is not fair. A combined recipe (edge threshold + L∩F) was not tried.

## 8. Verdict

**Does consistency-filtered evidence generalize? No, not as pre-registered.**
- The primary method `linsi|cons0.7` is a BAliBASE effect:
  - −2.0 on fresh BAliBASE draws;
  - −0.01 on the other 20 held-out protein datasets;
  - +0.66 on HomFam;
  - +4.7 on nucleotides.
- The frozen reference-free gate is at chance on proteins (54%).

That part of the project is **not promising**.

**What survives.**
- `L ∩ FFT-NS-2 --op 3` is a cheap, protein-only recipe that held up on held-out simulated and BAliBASE data
  (protein −1.28, 23/1/4, p = 0.0003). HomFam is neutral, trees are unchanged, and it must be switched off for
  nucleotides.
- The post-hoc finding that pruning weakly supported GCM edges helps where MAGUS is recall-limited is new: −5.8 on
  simulated proteins, neutral on nucleotides. It points to a mechanism (MCL fragmentation) rather than to evidence
  precision.

**For a 4-week CS581 project: unclear-to-promising, but only after a pivot.** The pivot is "support-aware GCM
graphs": edge-support thresholds or support-weighted edges.
1. Pre-register k (or a weighting) and test on fresh simulated/HomFam/ROSE draws.
2. Test the MCL-fragmentation mechanism directly: cluster counts and the trace's column splits.
3. Combine with L∩F on proteins.
4. Report trees, which none of the alignment gains here improved.

Honest limits:
- one MAGUS draw per dataset;
- HomFam is scored on 5–20 seeds;
- BAliBASE and nucleotide timings are contended;
- the 20-backbone and TCS baselines were not run.
