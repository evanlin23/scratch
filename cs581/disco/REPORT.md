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

**Tools** (built under `/opt/tools`; see `code/README.md`):
- DISCO v1.4.1.
- ASTRAL-Pro 2: ASTER `astral-pro` v1.16.2.4, the last commit before ASTRAL-Pro3.
- ASTRAL-IV: ASTER v1.25.4.8, used for "ASTRAL-DISCO".
- ASTRID: wASTRID/internode v0.0.7 `--preset vanilla`, which is ASTRID-2-style internode distances plus FastME. The DISCO paper used ASTRID-2 itself.
- These versions differ from the paper's, which used ASTRAL-III and older ASTRAL-Pro.

**Baseline reproduction** (DISCO data, 50 genes = the paper's default gene count, 100 bp, 10 replicates). Paper values are figure means read off the plots (±0.01; the paper has no RF tables).

REPRO_TABLE

## 4. Method (what was implemented)

The code is `code/discor.py` and `code/run_rep.py`.

**1. S0.**
- S0 is the ASTRAL-Pro 2 tree on the multi-copy gene trees.
- It must be rooted for LCA mapping. I rooted it on the edge that minimises the total DL cost of 100–200 randomly sampled gene trees, so no outgroup is assumed.
- This gave the correct root in REPLACE_S0ROOT.

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

RESULTS

## 6. Runtime

RUNTIME

## 7. Verdict

VERDICT
