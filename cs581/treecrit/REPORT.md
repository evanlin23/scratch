AI-assisted (Claude), exploration code for CS581 project

# Which alignment errors cost FastTree accuracy? Why ~6.5 SP points buy only ~0.3 RF (treecrit)

CS581 (UIUC), Fall 2026. Branch `claude/cs581-treecrit`. Code in `code/`, data in `data/`, generated tables in
`results/partA_{simhigh,protein,all}.md`. Part B pre-registration: `PREREG.md`.

Course question (MSA lecture): *"Which alignment criteria are predictive of tree accuracy?"* Project question: the
GCM evidence filters (es4, recipe, vote hard-bb) lower SP error on hard simulated proteins (SIMHIGH) by about
6.5 points, yet FastTree nRF moves by only about 0.3 points. Why?

## Plain-language answer (Part A)

1. **Not all alignment errors cost the same.** Within a dataset, one point of SPFP (residues wrongly put in the
   same column) costs FastTree about **0.15 nRF points**. One point of SPFN (true homologies that are left in
   separate columns, mostly over-splitting) costs about **0.03**, five times less. A missed pair just becomes
   missing data. A false pair is a fake shared character that actively votes for wrong groupings.
2. **The filters remove the cheap kind of error.** es4, recipe and hard-bb remove about 13 SPFN points but only
   about 1 SPFP point. At these prices that is worth 13 × 0.03 + 1 × 0.15 ≈ **0.3–0.5 nRF points**, which is
   what we observe (−0.27 to −0.51 on SIMHIGH). So "6.5 SP points for 0.3 RF" is not a paradox: the SP gain is
   almost all in the cheap currency. (The true alignment is worth ~3.5 RF points over MAGUS because it also
   removes ~17 SPFP points.)
3. **The expensive errors are mostly out of reach of a support filter.** About half of the false pairs in a
   MAGUS-style merge come from errors already present in the 25 subset alignments: subset columns that mix
   residues of different true columns, which correct cross-subset edges then carry into the final column.
   Most of the rest ride on edges that nearly every backbone supports (posterior > 0.99): the backbones make the
   same mistake, so voting cannot see it. Only a few percent of false pairs sit on low-support edges, which is
   all a support filter can delete.
4. **The errors that predict tree error best are between close relatives, and the merge never touches them.**
   Across datasets, misalignment between sister taxa (true-tree cherries) is the strongest predictor of
   absolute nRF (adds 36% explained variance on top of the true-alignment tree's error; R² 0.24 → 0.60). 99% of
   cherries fall inside one MAGUS subset, so their alignment is fixed before the merge step. All merge variants
   leave it unchanged (Δ cherry SPFN ≈ −0.07 points, vs −30 points of SPFN between random taxon pairs).
5. **Within a dataset, tree error is mostly idiosyncratic.** No single measure (SP, TC, column precision,
   over-split/over-merge fractions, PI-restricted errors, gap placement, length, Fitch homoplasy on the true
   tree, distance distortion) explains more than ~6% of the within-dataset variance among estimated
   alignments; cross-validated R² is ≤ 0.04. Two different alignments of equal SP quality give trees that
   differ by SD 0.8–1.1 RF points, as much as any two alignments of the same data. Two independent MAGUS runs
   on the same data differ by SD 2.3 points. The tree depends on *which* residues are wrong, and no column
   summary captures that. Only the large contrast between estimated and true alignments is well predicted
   (within-draw R² ≈ 0.5 once the true alignment is included).

**Which measure best predicts tree error?**
- *Changes within a dataset* (the question for a method change): **SPFP** is the measure that matters per
  point. The SPFN + SPFP price model explains 6% of within-draw variance among estimated alignments, and 53%
  with the true alignment included. Single best measures here are column precision, TC and SP error, all with
  within-draw r ≈ 0.22 (R² ≈ 0.05).
- *Absolute nRF across datasets*: **alignment error between the closest relatives** (p-distance distortion
  or SPFN between cherry taxa), then misplaced residues and SPFP. SP error itself ranks in the middle.

## 1. Data (Part A)

- Every tree row with a paired alignment on `claude/cs581-gcmtrees{,-h1,-h2,-h3}` and `claude/cs581-gcmvote-*`
  (h1, h2, h8, h10/h10b/h10c, t03–t20; p21–p49 added later, see §6) whose replicate is in a rep bank
  (`code/reps.tsv`, `code/gather.py` → `data/trees.jsonl`). SIMHIGH_R13 (h7) has no bank and is excluded.
- 37 draws in the first pass: SIMHIGH 25 draws (R1–R8 from gcmtrees and R3–R20 from gcmvote; R3–R8 have two
  independent MAGUS draws), SIMMOD 8, ROSE 1000M1 R0–R3 (DNA) 4. 302 estimated alignments with a FastTree
  tree, plus 37 true-alignment trees and 10 oracle split(X) trees (gcmtrees, SIMHIGH R1–R5).
- **Alignments were regenerated** from the bank (MAGUS's subsets and backbones; `code/regen.py` runs
  `vote.py` and `gg.py` merge-only, one thread each, 4 lanes). All 283 alignments with a logged FastSP row
  reproduce it exactly (SPFN and SPFP to 1e-4; `data/aln_logged.jsonl`). The masked vote_hard copies are
  rebuilt by vote.py's own masking. SIMHIGH_R17's es4 tree used gg.py's `linsi#es4`, regenerated as such.
  Two bank reps (R4, R17 of gcmvote) lacked the gcmgen `aligned/linsi` layout; it was rebuilt from
  `inputs/backbones` (identical sequence sets checked).
- FastTree `-lg -gamma` (proteins) / `-nt -gtr -gamma` (DNA); nRF vs the true tree, in %.

## 2. Measures (`code/measures.py`, `code/measures2.py`)

Per alignment, against the true alignment (residue level, exact):
SPFN, SPFP, SP error, TC (column recall), column precision, fraction of true columns over-split, fraction of
estimated columns over-merged, split residues (outside the plurality estimated column of their true column),
misplaced residues (outside the plurality true column of their estimated column), the same restricted to
parsimony-informative (PI) columns, PI column count ratio, length ratio, internal gap-opening ratio, residue
error next to true indels vs away from them, per-taxon error (90th percentile), Fitch parsimony excess of the
alignment on the true tree (all / PI columns), p-distance distortion (random taxon pairs, true-tree cherries),
rank error of p-distances, compared-site ratio, SPFN/SPFP restricted to residue pairs between cherry taxa or
random taxon pairs. Our SPFN/SPFP/TC equal FastSP's exactly.

## 3. Results (SIMHIGH, 236 estimated alignments, 25 draws; full tables in `results/partA_simhigh.md`)

### 3.1 Within-draw association (draw fixed effects; 95% cluster-bootstrap CI over draws)

| measure | within r [95% CI] | within R² | paired ΔRF vs Δx (vs MAGUS) ρ (p) | across-draw ρ(RF − RF_true, x) |
|---|---|---|---|---|
| 1 − column precision | +0.23 [+0.03, +0.45] | 0.05 | +0.25 (3e-4) | +0.16 |
| 1 − TC | +0.22 [+0.02, +0.43] | 0.05 | +0.21 (0.002) | +0.15 |
| over-merged est. columns | +0.22 [+0.02, +0.44] | 0.05 | +0.23 (0.001) | +0.11 |
| Fitch excess on true tree | +0.22 [+0.00, +0.43] | 0.05 | +0.14 (0.04) | +0.24 |
| SP error | +0.21 [+0.03, +0.41] | 0.04 | +0.25 (3e-4) | +0.29 |
| misplaced residues | +0.20 [+0.09, +0.30] | 0.04 | +0.13 (0.07) | +0.42 |
| p-distance distortion, random pairs | +0.20 [+0.10, +0.29] | 0.04 | +0.18 (0.009) | +0.31 |
| SPFN | +0.20 [−0.00, +0.41] | 0.04 | +0.22 (0.001) | +0.21 |
| split residues | +0.20 [+0.02, +0.41] | 0.04 | +0.22 (0.002) | +0.23 |
| SPFP | +0.10 [−0.06, +0.24] | 0.01 | +0.09 (0.2) | +0.38 |
| over-split true columns | +0.06 [−0.09, +0.27] | 0.00 | −0.19 (0.006) | +0.33 |
| p-distance distortion, cherries | +0.07 [−0.12, +0.23] | 0.00 | −0.25 (2e-4) | +0.56 |

Leave-one-draw-out cross-validated within-draw R²: SP error 0.03, SPFN + SPFP 0.03, split + misplaced
residues 0.03, TC 0.03, Fitch 0.02, distance distortion 0.04; no model reaches 0.05. The proteins-only (33
draws) and all-data (37 draws) tables give the same picture (`results/partA_protein.md`, `results/partA_all.md`).

SPFP's raw within-draw correlation is weak (0.10) only because the variants barely change SPFP (SD across
variants ≈ 1 point vs ≈ 6 for SPFN). Per point it is the expensive error (§3.2).

### 3.2 Price of each error type (within-draw OLS on SPFN and SPFP together; nRF points per error point)

| rows | SPFN | SPFP | within R² |
|---|---|---|---|
| estimated alignments only (236, 25 draws) | +0.030 [+0.006, +0.052] | +0.175 [+0.025, +0.323] | 0.06 |
| + true alignment (261) | +0.028 [+0.004, +0.050] | +0.151 [+0.093, +0.209] | 0.53 |
| SIMHIGH R1–R5 incl. oracle split(X) (35) | +0.062 [+0.028, +0.110] | +0.106 [+0.093, +0.121] | 0.64 |
| residue-level: split vs misplaced residues (+ true) | +0.029 [−0.010, +0.069] | +0.277 [+0.179, +0.373] | 0.53 |

Proteins (33 draws): SPFN +0.029, SPFP +0.164; all data (37 draws): +0.029, +0.148. Same ratio everywhere: a
false pair costs 5–6× a missed pair (2–3× in the small oracle set).

### 3.3 Observed vs price-predicted ΔnRF (SIMHIGH, paired vs MAGUS)

| variant | n | ΔSP err | ΔSPFN | ΔSPFP | observed ΔnRF | predicted from prices |
|---|---|---|---|---|---|---|
| es4 | 17 | −6.69 | −13.72 | +0.33 | −0.27 | −0.33 |
| recipe | 25 | −7.30 | −13.60 | −1.00 | −0.51 | −0.53 |
| vote hard-bb | 17 | −6.69 | −12.49 | −0.88 | −0.43 | −0.48 |
| vote hard | 17 | −6.55 | −12.11 | −0.99 | −0.18 | −0.49 |
| es3 | 25 | −6.25 | −12.97 | +0.47 | −0.06 | −0.29 |
| hard filter (gg `linsi&fftns2-op3`) | 25 | −2.85 | −4.86 | −0.84 | −0.55 | −0.26 |
| vote soft4 | 17 | −2.97 | −4.51 | −1.43 | −0.08 | −0.34 |
| vote soft | 17 | +0.18 | +1.57 | −1.20 | +0.32 | −0.14 |
| oracle split(MAGUS) (SPFP → 0) | 5 | −8.50 | 0.00 | −17.01 | −1.10 | −2.57 |

The price model gets the size of the filter effects right (−0.3 to −0.5). Per-variant residuals are within the
noise (§3.5). It overpredicts the oracle that removes all false pairs (−1.1 observed vs −2.6), so the cost of
false pairs is not linear all the way down: the first false pairs removed are worth more than the last.

### 3.4 Where the false pairs come from (`code/fpsource.py`, SIMHIGH R3 and R10, gcmvote draws)

Every residue pair in a final column is classed by how the merge produced it.

| source of false pairs | MAGUS R10 | hard-bb R10 | MAGUS R3 | hard-bb R3 |
|---|---|---|---|---|
| within one subset alignment | 3.5% | 3.2% | 2.8% | 2.6% |
| correct edge (posterior > 0.99) joining impure subset columns | 48.0% | 53.4% | 32.0% | 34.5% |
| correct edge, lower posterior, impure columns | 5.3% | 6.0% | 6.8% | 8.1% |
| wrong edge, posterior > 0.99 | 16.0% | 16.7% | 11.1% | 12.4% |
| wrong edge, posterior 0.5–0.99 | 14.2% | 16.7% | 31.9% | 35.5% |
| wrong edge, posterior < 0.5 | 13.0% | 3.8% | 15.3% | 6.8% |
| transitive (no edge between the two nodes) | 0.1% | 0.2% | 0.1% | 0.2% |
| SPFP of the alignment | 13.5 | 12.8 | 19.3 | 18.3 |

"Impure" = the two subset columns' majority true columns agree, but the subset columns themselves mix
residues of several true columns (~8–9% of residues). Hard-bb deletes most low-posterior wrong edges (13–15% → 4–7% of false
pairs) but leaves everything else. Since hard-bb also keeps most true pairs, SPFP drops only ~1 point.

**Oracle check (true subset alignments, same subsets and backbones, GCM merge; `code/truesub.py`, SIMHIGH
R1–R10).** With correct subset alignments the merged SPFP falls from 13–19% to 4–9%, so about half of MAGUS's
false pairs originate in the subset alignments. FastTree on these oracles: *pending, see §3.6*.

### 3.5 Noise floor

| comparison | n | SD of ΔnRF |
|---|---|---|
| vote_hard vs its masked copy (3–8 of ~8,000 columns removed) | 17 | 0.45 |
| two different alignments of one draw within 0.25 SP points | 29 | 1.10 |
| all within-draw deviations among estimated alignments | 236 | 0.79 |
| variant − MAGUS | 219 | 1.24 |
| two independent MAGUS draws of the same dataset (R3–R8) | 6 | 2.33 |

FastTree is deterministic: re-running it on two MAGUS alignments reproduced the logged trees exactly. Tiny
alignment edits move nRF by SD ≈ 0.45. Different alignments of equal quality move it by about 1 point, as much
as alignments of very different quality. So most within-draw tree differences are not "noise" in the
sampling sense. They are a real but practically unpredictable response to *which* residues are misaligned.
The tiny-edit noise alone would cap within-draw R² at ~0.84; the measures reach 0.06.

### 3.6 Oracle trees (true subset alignments + GCM merge)

*Pending (FastTree running).*

## 4. Part B

See `PREREG.md` (pushed before any split tree). Results will be added here.
