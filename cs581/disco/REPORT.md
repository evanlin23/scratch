# DISCO-R pilot: does species-tree-guided rooting and tagging improve DISCO?

*CS581 (Fall 2026) project pilot, run 2026-10-09/10. This is evidence for a go/no-go decision, not a finished study.*

## 1. Question

**Source.** Lecture "Phylogenomics part 2", slide 35, "(Some) Open Questions": *"Can we improve DISCO?"* See also `literature/slide_open_problems.md` item 9 and `literature/wide_scan.md` §2.2 D.

**How DISCO works** (Willson et al., *Syst Biol* 2022, doi:10.1093/sysbio/syab070):
- It roots each gene-family tree in isolation, choosing the root with ASTRAL-Pro's MinDL score (fewest duplications from species-set overlap).
- It tags each node by species overlap.
- It then cuts off the smaller subtree at every duplication node to make single-copy trees.
- Those trees go to ASTRAL or ASTRID.

**The idea, DISCO-R:**
1. Estimate a first-pass species tree S0 with ASTRAL-Pro 2.
2. Re-root each gene-family tree by duplication–loss (DL) parsimony against S0 under LCA mapping (the NOTUNG approach), then re-tag it.
3. Run DISCO's decomposition.
4. Re-estimate the species tree. Optionally iterate.

**Questions tested:**
- (a) Does DISCO-R tag gene trees more accurately?
- (b) Does that make the species tree more accurate?
- (c) How much headroom is there, using an oracle that roots against the *true* species tree?

## 2. Prior art

Details and DOIs are in [`prior_art.md`](prior_art.md). Summary:

- **Not found published.** I found no pipeline that does all of: estimate S0, re-root/re-tag gene-family trees by DL parsimony against S0, decompose DISCO-style, re-estimate.
- **Closest prior art.**
  - OrthoFinder (Emms & Kelly 2019, doi:10.1186/s13059-019-1832-y) builds a STAG species tree, roots it with STRIDE, then roots gene trees against it. It never decomposes and re-estimates the species tree.
  - SpeciesRax and AleRax (doi:10.1093/molbev/msab365, doi:10.1093/bioinformatics/btae162) co-estimate gene and species trees under full DTL likelihood. They have no tagging or decomposition step.
- **ASTRAL-Pro 2 and 3 (ASTER).** I read the source. Rooting and tagging use only species overlap, with no species tree and no re-tagging rounds.
- **DISCO+QR and QR-STAR** (doi:10.1093/bioadv/vbad015, doi:10.1089/cmb.2023.0185) root the *species* tree, not the gene trees.
- **wQFM-DISCO** (doi:10.1093/bioadv/vbae189) keeps DISCO's tagging.
- **Must read before claiming novelty:** Parsons, Liu, Dua, Markin & Molloy, bioRxiv 2026 (doi:10.64898/2026.01.20.700722), which studies ASTRAL-Pro tagging correctness under DLCoal. I did not read the full text.
- **Reconciliation building blocks:** NOTUNG (doi:10.1093/bioinformatics/bts386), ecceTERA (doi:10.1093/bioinformatics/btw105) and GeneRax (doi:10.1093/molbev/msaa141).

## 3. Data, tools, baseline reproduction

**Data** (Illinois Data Bank, downloaded in full):
- **DISCO data**, doi:10.13012/B2IDB-4050038_V1 (`trees.tar.gz`, 5.7 GB).
  - 101 species (100 plus an outgroup), 1000 gene families per replicate, 10 replicates.
  - Gene trees are FastTree trees from 100 bp alignments.
  - SimPhy true gene trees and locus trees are included.
  - Conditions used: `default` (dup 5e-10, loss/dup 1, Ne 5e7, AD 20%), `ils_1e4` (AD 0%), `ils_2e8` (AD 50%), `gdl_1e-10_1`, `gdl_1e-9_1`.
- **DISCO+QR data**, doi:10.13012/B2IDB-5748609_V1.
  - 21 species, 1000 gene families, 10 replicates.
  - Dup rate 1e-10, 5e-10 and 1e-9; loss/dup 0 or 1; low ILS (~20% AD) or high ILS (~70% AD, `_hILS`).
  - Gene trees from 50 bp and 100 bp alignments, plus true gene trees and locus trees.
- **Data problems:**
  - `20_gdl_1e-12_1` has a 101-taxon `s_tree` with 21-taxon gene trees (a packaging error), so I skipped it.
  - In `20_gdl_1e-9_0` (no loss) gene trees average about 1800 leaves (up to 9500). Only 3 replicates were run there.

**Coverage actually run on the 21-species data:**
- 1000 genes, with 50 bp and 100 bp gene trees:
  - 10 replicates each for dup 1e-10 (loss 0 and 1), dup 5e-10 (loss 1) and dup 1e-9 (loss 1), at both ILS levels;
  - 10 replicates for 5e-10 with no loss at high ILS, but only 2 at low ILS;
  - 3 replicates for 1e-9 with no loss at high ILS.
- True gene trees: replicates 01–05.
- 100-gene subsets: 10 replicates in all 10 cells of dup {1e-10, 5e-10, 1e-9} × loss {0, 1} × ILS {low, high}, excluding 1e-9 with no loss.

**Tools** (built under `/opt/tools`; see `code/README.md`):
- DISCO v1.4.1.
- ASTRAL-Pro 2: ASTER `astral-pro` v1.16.2.4, the last commit before ASTRAL-Pro3.
- ASTRAL-IV: ASTER v1.25.4.8, used for "ASTRAL-DISCO".
- ASTRID: wASTRID/internode v0.0.7 `--preset vanilla`, which is ASTRID-2-style internode distances plus FastME. The DISCO paper used ASTRID-2 itself.
- These versions differ from the paper's, which used ASTRAL-III and older ASTRAL-Pro.

**Baseline reproduction** (DISCO data, 50 genes = the paper's default gene count, 100 bp, 10 replicates). Paper values are figure means read off the plots (±0.01; the paper has no RF tables).

| condition (paper figure) | ASTRAL-Pro 2: ours / paper | ASTRID-DISCO: ours / paper | ASTRAL-DISCO: ours / paper |
|---|---|---|---|
| default, 50 genes (Fig 3, Fig 5) | 0.087 / ≈0.08–0.083 | 0.082 / ≈0.08 | 0.100 / ≈0.105 |
| dup 1e-10, 50 genes (Fig 8) | 0.094 / ≈0.095 | 0.079 / ≈0.075 | 0.091 / – |
| dup 1e-9, 50 genes (Fig 8) | 0.084 / ≈0.085 | 0.090 / ≈0.083 | 0.106 / – |
| AD 0% (Ne 1e4), 50 genes (Fig 7) | 0.065 / ≈0.06 | 0.059 / ≈0.055 | 0.074 / – |
| AD 50% (Ne 2e8), 50 genes (Fig 7) | 0.174 / ≈0.165 | 0.135 / ≈0.14 | 0.183 / – |
| default, 1000 genes (Fig 3), reps 01–04 only | 0.028 / ≈0.05 | 0.023 / ≈0.045 | 0.033 / ≈0.05 |

**Reading the reproduction.**
- With 50 genes, every number is within about 0.01 of the paper's figure means.
- The paper's qualitative ranking holds: ASTRID-DISCO is better than or equal to ASTRAL-Pro, which is better than ASTRAL-DISCO, and the gap widens under high ILS.
- With 1000 genes, our errors are about half the paper's. Possible reasons:
  - only 4 of 10 replicates finished in time (about 20 min each);
  - newer ASTRAL tools;
  - I used the first 1000 gene trees, and the paper's gene subsets are unknown.
- The 50-gene condition is the paper's default, so I consider the baseline reproduced.


## 4. Method (what was implemented)

The code is `code/discor.py` and `code/run_rep.py`.

**1. S0.**
- S0 is the ASTRAL-Pro 2 tree on the multi-copy gene trees.
- It must be rooted for LCA mapping. I rooted it on the edge that minimises the total DL cost of 100–200 randomly sampled gene trees, so no outgroup is assumed.
- This gave the correct root in 240/243 runs on the 21-species data with 1000 genes, 92/100 with 100 genes, and 40/54 on the 101-species data. The oracle runs below show that this rooting error hardly matters.

**2. Root.**
- For each gene-family tree, every rooting is scored by DL parsimony against S0 with an O(n) rerooting dynamic program: LCA mapping, duplication cost 1.5, loss cost 1 (NOTUNG defaults).
- Ties are broken by DISCO's own MinDL score.
- Input trees are first re-rooted on a random edge so that no method inherits the input root.

**3. Tag.**
- `DISCO-R`: DISCO's overlap rule (a duplication iff the children's species sets overlap) at the new root.
- `DISCO-R-lca`: full LCA-reconciliation tags (a duplication iff M(v) = M(child)).

**4. Decompose and estimate.**
- DISCO v1.4.1 `decompose()`; trees with fewer than 4 leaves are dropped.
- The single-copy trees go to ASTRID and ASTRAL.
- `-it2` repeats steps 1–4 once, using the ASTRAL-DISCO-R tree as the new S0.
- `-oracle` uses the true species tree as S0. This bounds what better rooting against a species tree can give.

**Tagging truth.**
- A locus-tree node is a speciation iff its height equals a species-tree node height (relative tolerance 1e-6). SimPhy locus trees are ultrametric on the species-tree time scale.
- Two genes are orthologs iff the LCA of their loci in the locus tree is a speciation.
- **Pair accuracy** is the fraction of cross-species gene pairs whose ortholog/paralog call, from the LCA node's tag in the rooted, tagged gene tree, matches the truth.
- **Root accuracy** is defined only for true gene trees. It is restricted to multi-copy trees and compares the chosen root bipartition with SimPhy's.

**Statistics.**
- Species-tree error is normalised RF against the true species tree (all trees are binary, so FN = FP = RF).
- Comparisons are paired by replicate.
- W/T/L counts how often A is better / tied / worse. The tie band is "identical RF"; RF moves in steps of 1/18 at 21 species and 1/98 at 101 species.
- Two-sided Wilcoxon signed-rank test on the non-tied pairs.
- **Hold-out.** Replicates 06–10 were never looked at during development. Only replicate 01 of four conditions was used to debug, and no parameter was tuned. Results are reported for all replicates and for 06–10 alone.

## 5. Results

Full tables are in `results/`:
- `qr_final.md`: 21 species, 1000 genes, 50 bp / 100 bp / true gene trees;
- `qr_k100.md`: 21 species, first 100 genes, 50 bp;
- `disco_final.md`: 101 species;
- `pooled.md`: all experiments pooled;
- `*_rf.csv`: per-replicate RF values.

### 5.1 Tagging: DISCO-R tags more accurately, but only by about 1 point

Pair accuracy is the fraction of cross-species gene pairs whose ortholog/paralog call is correct. Each row pools replicates over all conditions of that experiment.

| experiment | DISCO | DISCO-R | DISCO-R, true S (oracle) | DISCO-R-lca | DISCO-R better in |
|---|---|---|---|---|---|
| 21 sp, true gene trees (low/high ILS, dup 1e-10…1e-9) | 0.934 | **0.948** | 0.948 | 0.872 | 51/53 reps |
| 21 sp, 1000 genes, 100 bp | 0.875 | **0.887** | 0.887 | 0.768 | 88/95 |
| 21 sp, 1000 genes, 50 bp | 0.851 | **0.860** | 0.860 | 0.727 | 83/95 |
| 21 sp, 100 genes, 50 bp | 0.852 | **0.863** | 0.863 | 0.736 | 83/100 |
| 101 sp, 50 genes, 100 bp (5 conditions) | 0.862 | **0.871** | 0.872 | 0.758 | 36/54 |

**Gene-tree root accuracy.**
- This uses true gene trees, multi-copy families only, with the input rooting randomised.
- DISCO picks the true root in **22%** of families and DISCO-R in **77%** (99% at low ILS with dup 1e-10).
- That large rooting gain becomes only +1.4 points of pair accuracy. Most of DISCO's wrong roots are *ties under MinDL*, so they produce the same overlap tags.

**Where DISCO-R loses.** It is slightly worse in one condition: 101 species, dup 1e-9, loss 1 (0.862 vs 0.870). There the higher ortholog recall is outweighed by lower precision: estimated-gene-tree error looks like extra duplications and losses to DL parsimony.

**LCA tagging is clearly worse** (−8 to −12 points). Under ILS and gene-tree error, LCA reconciliation calls many spurious duplications. Overlap tagging at the reconciled root is the right combination.

**The species tree is not the bottleneck.** Using the true species tree (oracle) instead of S0 changes pair accuracy by at most 0.001. This holds even though S0 was mis-rooted in 26% of the 101-species runs.

### 5.2 Species tree: no reliable improvement

Normalised RF against the true species tree, paired by replicate. W/T/L = DISCO-R better / tied / worse (tie = identical RF). Two-sided Wilcoxon signed-rank test.

| experiment | ASTRID-DISCO-R vs ASTRID-DISCO | ASTRAL-DISCO-R vs ASTRAL-DISCO | oracle ASTRID vs ASTRID-DISCO |
|---|---|---|---|
| 21 sp, 1000 genes, 100 bp (n=95) | +0.0006, 4/85/6, p=1 | +0.0012, 4/86/5, p=0.78 | +0.0006, 4/85/6 |
| 21 sp, 1000 genes, 50 bp (n=95) | −0.0023, 8/81/6, p=0.36 | **−0.0117, 20/72/3, p=0.0004** | −0.0023, 8/81/6 |
| 21 sp, 100 genes, 50 bp (n=100) | −0.0028, 20/63/17, p=0.77 | −0.0039, 22/65/13, p=0.20 | −0.0022, 19/63/18 |
| 101 sp, 50 genes (n=50) | **−0.0049, 24/16/10, p=0.041** | −0.0010, 21/10/19, p=0.67 | −0.0053, 25/13/12, p=0.038 |
| 101 sp, 1000 genes, default (n=4) | +0.0026, 0/3/1 | −0.0026, 1/3/0 | – |
| 21 sp, true gene trees (n=53) | +0.0021, 1/49/3, p=0.63 | −0.0031, 5/46/2, p=0.45 | – |
| **all estimated-gene-tree runs (n=344)** | −0.0020, 56/248/40, p=0.32 | **−0.0042, 68/236/40, p=0.004** | −0.0019, 56/245/43, p=0.40 |
| **held-out reps 06–10 only (n=165)** | −0.0038, 30/119/16, p=0.11 | −0.0046, 28/117/20, p=0.097 | −0.0033 |

**Caveats on these numbers:**
- The pooled rows are *not independent*. The 21-species "100 genes" and "1000 genes, 50 bp" rows share gene trees, so their pooled p-values are optimistic.
- Per condition (10 replicates each), no single condition reaches p < 0.05 for ASTRID-DISCO-R vs ASTRID-DISCO.
- **High-GDL/high-ILS cells, the case where we expected gains.** ASTRID: −0.006 at 1000 genes and 50 bp; **+0.011** (worse) with 100 genes. ASTRAL: −0.028 at 1000 genes and 50 bp; **+0.039** (worse) with 100 genes. In short, the gains do not concentrate there.
- **Second iteration (`-it2`).** It changes nothing: pooled −0.0015 vs DISCO.
- **The one consistent signal is ASTRAL-DISCO-R vs ASTRAL-DISCO on 50 bp gene trees** (−0.012 RF, 20/72/3, held-out p=0.035).
  - The oracle gives the *identical* result, so the cause is the re-rooting rule itself, not S0 quality.
  - DISCO-R emits slightly more single-copy trees: on average 1.0–2.6% more, and more in 56–89% of replicates. This may be what ASTRAL-IV benefits from.
  - This is a hypothesis; I did not test it.
- **Context.** ASTRID-DISCO, with or without R, beats ASTRAL-Pro 2 overall (−0.008 to −0.010, p<0.001), as in the DISCO paper. DISCO-R does not change which method wins.

**Why tagging gains don't transfer.**
- In these datasets, the pairs DISCO-R fixes are mostly near the root of multi-copy families. Re-rooting moves *which* subtree is cut off, but the large single-copy tree that dominates ASTRAL/ASTRID input is often the same.
- The quartet and internode summaries are already robust to the few percent of mis-tagged pairs. The DISCO paper said the same, and Zaman et al. 2026 found the same for rooting in STELAR.
- The oracle bound shows the headroom from better rooting *against a species tree* is at most about 0.005 RF here. Even a perfect S0 does not help.


## 6. Runtime

Wall-clock seconds per replicate, 1 thread, mean.

| setting | ASTRAL-Pro 2 | DISCO decomposition | DISCO-R extra (rooting S0 + re-root/decompose) | ASTRAL on decomposed trees | ASTRID |
|---|---|---|---|---|---|
| 21 sp, 100 genes | 0.4 | 0.3 | 4.9 | 0.5 | <0.1 |
| 21 sp, 1000 genes | 1–10 | 3–5 | 11–16 | 1–7 | <0.1 |
| 101 sp, 50 genes | 3–6 | 0.4 | 29 (S0 rooting is 16–33 of it) | 5–7 | <0.1 |
| 101 sp, 1000 genes | 252 | 10 | 73 | 258–308 | 0.1 |

- DISCO-R is pure Python: an O(n) DP per rooting pass, plus min-DL rooting of S0 over all 2n−3 edges using 100–200 sampled genes.
- The overhead is small next to ASTRAL / ASTRAL-Pro. It is large next to ASTRID-DISCO, which runs in seconds.
- A C++ or numba port would make it negligible.
- Gene families with about 1800 leaves (21 sp, dup 1e-9, no loss) took 10–30 min per replicate for the whole pipeline. That time is mostly the O(n²) tagging evaluation and ASTRAL-Pro, not DISCO-R.


## 7. Verdict

**Verdict: not promising as stated.**

- DISCO-R reliably improves *gene-tree* rooting (22% → 77% correct) and orthology tagging (+1 point of pair accuracy in 341/397 replicates).
- That does **not** become a species-tree improvement for ASTRID-DISCO, which is the strongest DISCO pipeline: pooled −0.002 RF, p=0.32; held out −0.004, p=0.11.
- The gain for ASTRAL-DISCO is small (−0.004 overall). It shows up mainly with 50 bp gene trees and 1000 genes, and is not explained by species-tree guidance (the oracle gives the same).
- With the true species tree as the guide, the ceiling is equally small. A better S0 or more iterations cannot rescue it.
- I would not pick this as a 4-week project with "beat DISCO on species-tree accuracy" as the success criterion.

**If it were still pursued, a 4-week plan that could produce a publishable-quality negative or positive result:**
1. **Week 1.** Read Parsons et al. 2026 in full. Add true-tag decomposition as an upper bound: decompose the estimated gene trees with *true* locus-tree tags. This answers whether *any* tagging improvement can help DISCO. If that bound is also flat, write the project up as "tagging is not DISCO's bottleneck".
2. **Week 2.** Isolate the ASTRAL-DISCO-R effect: the number and size of decomposed trees vs their orthology purity. Then try a decomposition that keeps *both* sides of a duplication when each is large, a "DISCO-max" variant. Decomposition is where the information is lost, not rooting.
3. **Week 3.** Harder regimes where tagging matters more: 1000-species DISCO data, missing data, fewer genes (10–50), and 50 bp. Use the ASTRAL-Pro3 and wQFM-DISCO baselines.
4. **Week 4.** Biological sanity check (DISCO's biological data, doi:10.13012/B2IDB-4050038_V1 `biological.tar.gz`) and the write-up.

**Better uses of the same infrastructure.** The code here (DL rooting DP, truth tagging from SimPhy locus trees, the paired batch harness) transfers directly:
- **wide_scan entry C, ASTRID-Pro: speciation-only internode distances.** It needs exactly this tagging. The finding that DISCO-R tagging raises orthology recall over MinDL tagging (pooled 0.60–0.69 → 0.65–0.74, with precision unchanged within ±0.015) is directly relevant there.
- **Measuring ASTRAL-Pro's consistency under tagging error** (slide 35). The pair-accuracy harness already quantifies tagging error per condition.

**Risks and limitations of this pilot.**
- **Tool versions.** ASTRAL-IV and wASTRID differ from the paper's ASTRAL-III and ASTRID-2. 1000-gene errors came out lower than the paper's.
- **Small n per condition.** 10 replicates; ties dominate at 21 species, where RF moves in steps of 0.056.
- **Fixed parameters.** One DL cost setting (1.5/1), with no support-based contraction of S0. Contraction could matter more at 1000 species.
- **Limited conditions.** The 101-species study used 50 genes and a single 1000-gene condition (4 replicates). dup 1e-9 with no loss was barely run (3 replicates, mean 1800 leaves per family).
- **Pooled p-values** combine non-independent experiments.
- **Prior art.** Not verified against Parsons et al. 2026 full text, or against RECOMB-CG/WABI 2025–26.

