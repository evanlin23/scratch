# DISCO-R prior-art check (2026-10-09)

Checked by web search and by reading code. DOIs were resolved via Crossref. Anything I could not confirm is marked **(unverified)**.

## Baselines

- **DISCO.** Willson, Roddur, Liu, Zaharias & Warnow, *Syst Biol* 71:610 (2022), doi:10.1093/sysbio/syab070. Correction: doi:10.1093/sysbio/syad003.
  - Roots and tags each gene tree in isolation with ASTRAL-Pro's MinDL heuristic, then decomposes it.
  - Consistent only if rooting and tagging are correct.
  - The paper does not measure tagging error itself.
- **ASTRAL-Pro.** Zhang et al., *MBE* 2020, doi:10.1093/molbev/msaa139.
- **ASTRAL-Pro 2.** Zhang & Mirarab, *Bioinformatics* 2022, doi:10.1093/bioinformatics/btac620.
- **ASTER / ASTRAL-Pro3.** Zhang, Nielsen & Mirarab, *MBE* 2025, doi:10.1093/molbev/msaf172.
  - Read in `src/astral-pro.cpp`: rooting and tagging use only species-set overlap (MinDL), and no species tree is involved.
  - `-T` only prints the rooted, tagged gene trees.
  - `-g` and `-c` only constrain the species-tree search.
  - Rounds (`-R`/`-S`) are search iterations and do not re-tag.
  - **No species-tree-guided re-tagging and no iteration.**

## DISCO variants (none re-root the gene trees)

- **DISCO+QR.** Willson et al., *Bioinf Adv* 2023, doi:10.1093/bioadv/vbad015. Roots the *species tree* after DISCO.
- **QR-STAR.** Tabatabaee, Roch & Warnow, *J Comput Biol* 2023, doi:10.1089/cmb.2023.0185. Roots the species tree.
- **wQFM-DISCO.** Hakim, Ratul & Bayzid, *Bioinf Adv* 2024, doi:10.1093/bioadv/vbae189. Changes the summary method and keeps DISCO's tagging.
- **wQFM-GDL.** Rafi et al., bioRxiv 2025, doi:10.1101/2025.04.04.647228. Quartets from orthologs plus paralogs. No guided re-tagging found **(unverified: full text not read)**.

## Studies of tagging error

- **Parsons, Liu, Dua, Markin & Molloy**, bioRxiv 2026, doi:10.64898/2026.01.20.700722.
  - Defines correct tagging under DLCoal.
  - Gives statistical properties of ASTRAL-Pro's tagging objective.
  - Measures tagging accuracy in simulations.
  - **This is the closest motivation and must be read in full before claiming novelty.** Whether it tests tagging against a species tree is **(unverified)**.
- **Zaman, Momin & Bayzid**, *PLoS Comput Biol* 2026, doi:10.1371/journal.pcbi.1014782. Sensitivity of STELAR (single-copy) to gene-tree rooting. Tangential.

## Pipelines that root gene trees with a species tree

- **OrthoFinder.** Emms & Kelly, *Genome Biol* 2019, doi:10.1186/s13059-019-1832-y.
  - STAG species tree → STRIDE root → species-tree-guided gene-tree rooting → orthologs.
  - This is DISCO-R steps 1–2, the closest prior art.
  - It does not decompose and re-estimate the species tree.
- **Possvm.** Grau-Bové & Sebé-Pedrós, *MBE* 2021, doi:10.1093/molbev/msab234. Species-overlap-based orthology. Different.
- **ROADIES.** Gupta, Mirarab & Turakhia, *PNAS* 2025, doi:10.1073/pnas.2500553122. Iterates over the number of genes with ASTRAL-Pro2. No re-rooting.

## Reconciliation-based rooting and co-estimation

- **NOTUNG.** Stolzer et al., *Bioinformatics* 2012, doi:10.1093/bioinformatics/bts386. DL-parsimony rooting against a species tree. Building block for step 2.
- **ecceTERA.** Jacox et al., *Bioinformatics* 2016, doi:10.1093/bioinformatics/btw105.
- **GeneRax.** Morel et al., *MBE* 2020, doi:10.1093/molbev/msaa141. Needs a fixed species tree.
- **SpeciesRax.** Morel et al., *MBE* 2022, doi:10.1093/molbev/msab365.
  - Species-tree search under DTL likelihood, with gene roots handled implicitly.
  - Overlapping in spirit (species tree ↔ gene tree reconciliation), but it is full likelihood, with no tagging, decomposition or ASTRAL/ASTRID.
  - DISCO's paper compared against it.
- **AleRax.** Morel et al., *Bioinformatics* 2024, doi:10.1093/bioinformatics/btae162. Co-estimation. Different method.

## Verdict

- The two-pass DISCO-R pipeline does not appear to be published. The pipeline is:
  1. estimate S0;
  2. re-root and re-tag the gene-family trees by DL parsimony against S0;
  3. DISCO-decompose;
  4. re-estimate the species tree.
- Each ingredient is classic (NOTUNG 2012, OrthoFinder 2019), so the novelty is the combination plus the evaluation.
- Coverage gap: RECOMB-CG and WABI 2025–26 proceedings and theses were not searched in depth.
