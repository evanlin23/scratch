# Prior art (search done 2026-10-09; DOIs checked against Crossref)

[A] = read in the abstract or full text; [M] = from memory or a secondary source (unverified).

| Paper | Venue/yr | DOI | Relevance |
|---|---|---|---|
| Legried, Molloy, Warnow, Roch, Polynomial-time statistical estimation of species trees under GDL | JCB 2021 / RECOMB 2020 | 10.1089/cmb.2020.0424 | ASTRAL-one and ASTRAL-multi are consistent under GDL (true gene trees, no ILS) [A]. ASTRID-multi, STAG and NJst appear only in simulations. |
| Markin & Eulenstein, Quartet-based inference is statistically consistent under DLCoal | Bioinformatics 2021 | 10.1093/bioinformatics/btab414 | Quartet methods are consistent under DLCoal [A]. Nothing on distance methods. |
| Hill, Legried, Roch, sample complexity of quartet methods under ILS+GDL | Ann Appl Probab 2022 | 10.1214/22-aap1799 | Quartet methods only [A]. |
| Willson, Roddur, Warnow, Comparing methods for species tree estimation with GDL | AlCoB 2021 | 10.1007/978-3-030-74432-8_8 | Simulation study; ASTRID-multi was competitive [M]. No theory. |
| Willson et al., DISCO | Syst Biol 2022 | 10.1093/sysbio/syab070 | ASTRAL-DISCO is consistent if rooting and tagging are correct [A]. ASTRID-DISCO is empirical only. |
| Zhang, Scornavacca, Molloy, Mirarab, ASTRAL-Pro | MBE 2020 | 10.1093/molbev/msaa139 | Consistent under GDL for correctly rooted and tagged trees. Imperfect rooting and tagging is left to future work [A]. |
| Zhang & Mirarab, ASTRAL-Pro 2 | Bioinformatics 2022 | 10.1093/bioinformatics/btac620 | Speed-up only [M]. |
| Emms & Kelly, STAG | bioRxiv 2018 | 10.1101/267914 | Closest-copy distance, then FastME. No proof [A/M]. |
| Molloy & Warnow, FastMulRFS | Bioinformatics 2020 | 10.1093/bioinformatics/btaa444 | Consistent under GDL when GDL is not adversarial [A]. Uses an RF-type criterion, not internode distances. |
| Vachaspati & Warnow, ASTRID | BMC Genomics 2015 | 10.1186/1471-2164-16-s10-s3 | Consistent under the MSC, single-copy genes only. |
| Allman, Degnan, Rhodes, unrooted STAR / NJst | IEEE/ACM TCBB 2018 | 10.1109/tcbb.2016.2604812 | Consistent under the MSC only [M]. |
| Rhodes, Nute, Warnow, NJst and ASTRID inconsistent under random missing data | arXiv 2001.07844 (2020) | none | Negative result for i.i.d. taxon deletion. Gene loss is a similar mechanism. |
| Morel et al., SpeciesRax (MiniNJ) | MBE 2022 | 10.1093/molbev/msab365 | Minimum-distance (STAG-like) starting tree. No consistency proof [A/M]. |
| Legried, Anomaly zones under GDL | arXiv 2023 | 10.48550/arXiv.2309.01663 | No anomaly zone for balanced 4-taxon trees, but caterpillars have one [A]. |

**Bottom line.** No paper found proves or disproves consistency under GDL or DLCoal for ASTRID-multi, NJst, STAG or any internode-distance method. No paper found defines a GDL-corrected internode distance (ortholog pairs only, counting only speciation nodes) with a proof. The search covered Crossref plus web searches but was not exhaustive; recent Warnow-lab and Mirarab-lab preprints still need a check.
