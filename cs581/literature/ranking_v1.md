# Ranked comparison of explored project ideas (first pass, 2026-10-10 03:05 UTC)

Sources: each session's `cs581/<dir>/REPORT.md` on its `claude/cs581-<dir>` branch (verdicts
extracted verbatim), the MAGUS-variant fan-out (50 held-out replicates, `variants/SUMMARY.md`) and the
measured end-to-end benchmark (`claude/cs581-e2e-*`, 4 of 14 datasets so far). Nothing is chosen.
Pending: fastmagus, sota2026, dldtm, models, rogue, lenhet (launched 01:20-01:50), gdl follow-up
(sec. 6.4), gdlcons atlas tables, camus and ml (reports still have placeholders).

| rank | idea (branch) | slide / instructor question | strongest result | beats strongest baseline? | novelty | verdict |
|---|---|---|---|---|---|---|
| 1 | **ASTRAL-Pro under GDL: consistency vs tagging error** (gdlcons) | Slide 35, verbatim: "Is ASTRAL-Pro consistent for GDL under ... error for rooting and tagging?"; "Which other methods are consistent under GDL?" | Exact 4-taxon formula: stock ASTRAL-Pro3 is **inconsistent under pure GDL even with true gene trees** (wrong 40/40, 20/20, 4/4 at 500/2k/10k families; margin -0.136 +- 0.004); consistent with true tags; tag-error threshold predicted 0.079, observed 0.10-0.15 | n/a (theory) | none found (Parsons et al. 2026 unread) | promising (Q1) |
| 2 | **ASTRID-Pro: GDL-corrected distance** (gdl) | Slide 35, verbatim: "distance correction for GDL so ASTRID/NJst are consistent?" | Theorem: ortholog-only internode distance is a tree metric with the species-tree topology iff no supercritical branch; counterexample matches prediction to 3 decimals | no: ties ASTRID-multi (p=0.90, pre-registered) | open | unclear: theory promising, method not |
| 3 | **Sample complexity of ASTRID/NJst vs ASTRAL** (samplecx) | Slide 35, verbatim | Caterpillar n=64, f=0.1: ASTRID needs 1189 genes vs ASTRAL 492; constant grows 4.9->15.3 vs 4.0->7.6 (n=8->64); suggestive, CIs overlap | n/a (evaluation) | open (Roch 2018 conjecture, untested beyond n=8) | promising as evaluation |
| 4 | **GTM blending** (gtm) | D&C deck: "develop a better DTM that allows blending"; "find a condition where GTM < TreeMerge" | GTM-Blend-ML vs GTM: -14.3 / -6.1 / -10.1 FN pts (20/0/0, p=9e-5) in simulation | in simulation only; loses to full IQ-TREE (+11 pts); n.s. on published data | blending DTMs exist; ML-blending not found | yes, if framed "blending fixes GTM when guide is poor" |
| 5 | **Soft-constraint MAGUS** (main track) | project list: improve MAGUS / merging alignments | self-soft vs MAGUS -0.78 pts (42/4/1, p=5e-13, n=47); vs MAGUS(Slow) -0.48 (32/8/7, p=1e-6); measured runtime +6-9% on 1000-seq data | yes for alignment; **no for trees** (FastTree tie, p=0.90) | new | modest (about 1/3 of MAGUS's gain over PASTA) |
| 6 | **Which alignment errors matter for trees** (alncrit) | MSA deck, verbatim: "Which alignment criteria are predictive of tree accuracy?" | at equal SPFN, real aligner errors cost 2-9x more tree accuracy than random errors (p=0.020) | n/a (evaluation) | open | unclear, promising angle |
| 7 | DISCO-R (disco) | Slide 35: "Can we improve DISCO?" | vs ASTRID-DISCO -0.002 RF, p=0.32 (n.s.); only vs ASTRAL-DISCO p=0.0004 | no | not found | not promising as stated |
| 8 | Consensus alignment (main track) | First day: "consensus alignments" | -0.78 vs chosen input, tie with MAGUS(Slow), worse than best input | no | partly | modest |
| - | Pairwise merging in PASTA (pairmerge) | First day: "merging two alignments" | exact DP vs OPAL n.s. either way | no | algorithm published | not promising |
| - | Quartet amalgamation (quartets) | Slide 35: "better quartet amalgamation" | ASTRAL-IV already reaches the best quartet score in 377/380 | no | (a),(b) published | not promising |
| - | Supertrees at scale (supertree) | First day: supertrees at 100k+ | ties ASTRAL-III at ~1/30 memory; SCS dominates | no | low | not promising |
| - | Forest+DTM (forest) | D&C deck | gains come from the decomposition; FastTree beats all | no | new but not useful | not promising |
| - | Linguistic models (ling) | Linguistics deck | new method loses; idea published (Canby 2024) | no | low | leaning not |
| - | Learned edge weights (local) | First day: deep learning | -0.08 on simulated, worse on proteins | no | new | not promising |
| - | MAGUS-lite / small backbones (local) | goal: faster | 15%+ error / +1.5-2 pts | no | - | not promising |

## Measured end-to-end runtime (same 4-core machine; 4 of 14 datasets)

| dataset | PASTA 1.8.3 | MAGUS | MAGUS(Slow) | self-soft | slow-soft |
|---|---|---|---|---|---|
| 1000L1 R0 | 23.4 min / 12.08% | 23.2 / 7.73 | 25.7 / 7.34 | 25.3 / **7.17** | 26.8 / 7.88 |
| 1000S3 R0 | 27.5 / 4.84 | 22.9 / 4.52 | 24.8 / 4.07 | 24.4 / **3.63** | 25.2 / 3.71 |
| 1000L2 R0 | 28.1 / 7.43 | 23.8 / 5.12 | 25.9 / **3.47** | 25.7 / 3.72 | 27.5 / 4.16 |
| BAliBASE 0101 (322 seqs) | 5.9 / 29.68 | 14.9 / 27.98 | 15.0 / 29.05 | 21.0 / **27.13** | 20.2 / 27.30 |

On 4 cores PASTA is about as fast as MAGUS on 1000 sequences (the paper's 2.5x speedup used 16-core
nodes); self-soft costs +6-9% over MAGUS there (about the cost of MAGUS(Slow)), +41% on the small protein set.
