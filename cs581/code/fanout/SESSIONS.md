# Fan-out worker sessions (launched 2026-10-09 19:13 UTC)

Each runs `worker.sh jobs_wN.txt claude/cs581-worker-N` in its own cloud
container (Sonnet; prompt = WORKER_PROMPT.md) and pushes results to its branch.

| worker | jobs | branch | session |
|---|---|---|---|
| 1 | jobs_w1.txt | claude/cs581-worker-1 | session_01RenQV9joJs7AwJubxkMLZN |
| 2 | jobs_w2.txt | claude/cs581-worker-2 | session_012Gam4VVwtNh7DoKRvx3B6E |
| 3 | jobs_w3.txt | claude/cs581-worker-3 | session_01F2kyZ7CpF1RoQzXk2JvS5C |
| 4 | jobs_w4.txt | claude/cs581-worker-4 | session_01McxGQrNS5jHLZ6U71i33bE |
| 5 | jobs_w5.txt | claude/cs581-worker-5 | session_01MQJbogSY6hT4UwkSokSmfx |
| 6 | jobs_w6.txt | claude/cs581-worker-6 | session_017HQsL4VsUsh8hRcJ25nqjw |
| 7 | jobs_w7.txt | claude/cs581-worker-7 | session_01Lg23wQ5RtFcNQPgejDFmss |
| 8 | jobs_w8.txt | claude/cs581-worker-8 | session_0184agjhEq4MAHXx6fyEVDVr |
| 9 | jobs_w9.txt | claude/cs581-worker-9 | session_019JV6jVL8RzRLn1VoKrDTQA |
| 10 | jobs_w10.txt | claude/cs581-worker-10 | session_01YYR7PKwA5uZWeB6dEjoDMC |
| 11 | jobs_w11.txt | claude/cs581-worker-11 | session_01MZ33magJd3UFGYJvqxc2Za |
| 12 | jobs_w12.txt | claude/cs581-worker-12 | session_01RmGYN29ByRYkDPbCFQFBhh |
| local | jobs_local_large.txt | claude/charming-pasteur-yl6k2v | orchestrating session |

Workers 7-12 were added at 22:12 UTC; the back halves of lists 1-5 were moved to them
(workers re-read their job list from the base branch before each job).

Collect results: `git fetch origin 'refs/heads/claude/cs581-worker-*:refs/remotes/origin/claude/cs581-worker-*'`
then copy each branch's `cs581/experiments/runs/*` (see aggregate.sh).

## Other sessions

| purpose | branch | session |
|---|---|---|
| maximum-likelihood exploration (literature, ML tree harness on the MAGUS datasets, pilots, report in `cs581/ml/REPORT.md`) | claude/cs581-ml | session_015LodS6GpCKvPXaDh4pLyME |
| consensus-MSA significance test, worker 1 (consensus_reps_w1.txt) | claude/cs581-consensus-1 | session_01Mb4gbSBYyvJeZFhya7BYLf |
| consensus-MSA significance test, worker 2 (consensus_reps_w2.txt) | claude/cs581-consensus-2 | session_019ZJfBbbRnLqeWx43FfaGmf |
| consensus-MSA significance test, worker 3 (consensus_reps_w3.txt) | claude/cs581-consensus-3 | session_01PD6mMSCA1xbEn3wPJXiJrm |
| consensus-MSA significance test, worker 4 (consensus_reps_w4.txt) | claude/cs581-consensus-4 | session_015dxKxLgEdAbLqs172yF3x3 |
| GTM-Blend exploration: validate GTM on the published DTM data, pilot a blending merger (`cs581/gtm/REPORT.md`) | claude/cs581-gtm | session_01Q8g6VEyCooDvxCvUjRZhMd |
| controlled runtime benchmark A (1000L3_R0, 1000M2_R0) | claude/cs581-timing-a | session_01BSVu1PkMDdK1YPBRF8CVrC |
| controlled runtime benchmark B (1000S2_R0, BBA0101_R0, RNASim_R0) | claude/cs581-timing-b | session_012V4tkNQ9BcgKy53szBAjDi |

## Exploration sessions (all slide open questions; launched 2026-10-09 22:46 UTC)

Each pilots one idea from `cs581/literature/slide_open_problems.md`: prior art, baseline
reproduction, paired pilot with runtime, and `cs581/<dir>/REPORT.md` with a verdict.

| idea | dir / branch | session |
|---|---|---|
| ASTRID under GDL: consistency probe + ASTRID-Pro | cs581/gdl, claude/cs581-gdl | session_01XKSHzvHypk3AcgwZBnTa2Y |
| ASTRID/NJst sample complexity + missing-data correction | cs581/samplecx, claude/cs581-samplecx | session_01G7cWtFUbRG4EbXXYgioKNN |
| DISCO-R (species-tree-guided root and tag) | cs581/disco, claude/cs581-disco | session_01Xmq8bA62x69pyxo7ZqohoW |
| CAMUS base tree / quartet filter | cs581/camus, claude/cs581-camus | session_019m2PfNkuX9mdkgzEKuJPbZ |
| Quartet amalgamation (ILS + HGT) | cs581/quartets, claude/cs581-quartets | session_019bppS3hdNG9kPeDaumUP6z |
| Linguistic evolution models and methods | cs581/ling, claude/cs581-ling | session_01SHe3R3FJaRhu1ZfHVdtBCg |
| Alignment criteria vs ML tree accuracy (incl. soft-MAGUS) | cs581/alncrit, claude/cs581-alncrit | session_01NgvCDj8KAVSiKtmcxkoUwQ |
| Merging two alignments (PASTA pairwise merger) | cs581/pairmerge, claude/cs581-pairmerge | session_01Aw3VmUYbcCjZkwsBjhiVCU |
| Supertrees at scale | cs581/supertree, claude/cs581-supertree | session_01RFzDUN2Z8wuQKvjvTBz7qB |
| Forest+DTM | cs581/forest, claude/cs581-forest | session_01QJAxbbXp6GDNfvQ1ZXGxJP |
| Learned evidence weights for GCM (deep learning x merging) | cs581/code/gcmx (local, this session) | orchestrating session |
| Faster MAGUS at equal accuracy (cheaper guide tree / fewer backbones + self-soft) | cs581/fastmagus, claude/cs581-fastmagus | session_016eRY2Kh3vGd3hmAfDk9d2R |
| Is MAGUS still SOTA in 2026? (TWILIGHT, FAMSA2, MUSCLE5, ...) | cs581/sota2026, claude/cs581-sota2026 | session_012Rgpxy2uEEPz2XSy4Tz17n |
| GDL consistency atlas (ASTRAL-Pro under rooting/tagging error; other methods under GDL/DLCOAL) | cs581/gdlcons, claude/cs581-gdlcons | session_017ko431PE4LESFrE9aX84WU |
| Deep-learning subset trees merged by GTM (DL for large-scale trees) | cs581/dldtm, claude/cs581-dldtm | session_01NgfgtKvrYSpihCVPR83VUo |
| Methods under new models (clock violation, realistic indels, tree shape) | cs581/models, claude/cs581-models | session_016r24PFxd7c5QBpCEq4hEba |
| Rogue taxa in alignment and tree estimation | cs581/rogue, claude/cs581-rogue | session_017Wru2Wv82AvDP88BymCMXk |
| Adding sequences with length heterogeneity (UPP/WITCH/EMMA) | cs581/lenhet, claude/cs581-lenhet | session_01RGx1ww8rFxYptzNLytBfNc |

The ML session (claude/cs581-ml) was asked to also cover "better ML heuristics" and
"scaling concatenation"; the GTM session covers DTM blending; consensus alignments run
on claude/cs581-consensus-1..4.

## Measured end-to-end benchmark (launched 2026-10-10 01:19 UTC)

PASTA 1.8.3 (paper's version) vs MAGUS vs MAGUS(Slow) vs self-soft / slow-soft MAGUS, every
pipeline run for real from unaligned sequences on the same idle 4-core machine
(`gcmx.e2e_bench`, `fanout/e2e_worker.sh`; results in cs581/experiments/e2e/<branch>.jsonl).

| jobs | branch | session |
|---|---|---|
| 1000L3_R0, 1000S1_R0 | claude/cs581-e2e-1 | session_0183ArYxNTput6waSGiaMsqv |
| 1000M1_R0, 1000M2_R0 | claude/cs581-e2e-2 | session_01VnRadyNWjuGQ3nVtXSXBCM |
| 1000L1_R0, 1000S2_R0 | claude/cs581-e2e-3 | session_01Bk4wiTG4VFFMZD88ta8BnZ |
| 1000S3_R0, 1000M3_R0 | claude/cs581-e2e-4 | session_0128Qxsj4s6xSbnPqBneEq9v |
| 1000L2_R0, 1000M4_R0 | claude/cs581-e2e-5 | session_016uB1sp1SwXKKqGhDJC4EUz |
| RNASim_R0, 16S.M_R0 | claude/cs581-e2e-6 | session_01XgjLDdFFSFeyfAcVFCo5LL |
| BBA0101_R0, BBA0190_R0 | claude/cs581-e2e-7 | session_01QLmz7p1YWy9UFKDR1QP8gy |

## Backbone-aligner benchmark (launched 2026-10-10 ~06:05 UTC)

Lead from cs581/sota2026 (merge pilot, n = 5): MAGUS with its GCM backbones aligned by Clustal
Omega instead of MAFFT L-INS-i had 1.6-1.9 points lower error on BAliBASE (5/5). Confirmation:
3 fresh MAGUS draws per dataset, paired merge-only comparison on each draw (clustalo, mafft-auto;
phase 2 adds linsi-noep, ginsi) plus measured end-to-end MAGUS with Clustal backbones
(`gcmx.bbtool_bench`, `fanout/bb_worker.sh`, `run_magus.py --gcmx-backbonetool`; results in
cs581/experiments/bbtool/<branch>.jsonl).

| jobs | branch | session |
|---|---|---|
| BBA0039_R0, BBA0067_R0 | claude/cs581-bbtool-1 | session_01PMMNFk2MXnWqcUH1VEDcua |
| BBA0081_R0, BBA0101_R0 | claude/cs581-bbtool-2 | session_01TWxgcSvAcr1WhRTf4ys2Hf |
| BBA0117_R0, BBA0134_R0 | claude/cs581-bbtool-3 | session_01LGeoTwv7CU7dgtZqEJYNbE |
| BBA0154_R0, BBA0190_R0 | claude/cs581-bbtool-4 | session_011WVrS4EKz2bB1BACnFCsv9 |
| RNASim_R0, 1000M2_R0 (nucleotide control) | claude/cs581-bbtool-5 | session_01TuxaQAE6YiMYWb4HzCi3t7 |
| 16S.M_R0, 1000L1_R0 (nucleotide control) | claude/cs581-bbtool-6 | session_01UZfnXAvSJ4zJVxqVt8DHJU |

## Overnight exploration (launched 2026-10-10 06:16 UTC onward)

| idea | dir / branch | session |
|---|---|---|
| What makes good GCM evidence: mechanism of the Clustal-backbone effect, mixed/more/cheaper backbones, MAFFT-only recipes | cs581/bbevidence, claude/cs581-bbevidence | session_0153wqpgcxzhkqFkK245ZgRF |
| Does the Clustal-backbone gain generalize: 10AA, HomFam, simulated proteins with tree accuracy | cs581/protbench, claude/cs581-protbench | session_01Bq4K1Uniua56tnkwkaNipt |
| EPA-ng accuracy drop on > 2,000-leaf placement subtrees (BSCAMPP open problem) | cs581/epang, claude/cs581-epang | session_01M5JKGxByGqhBP4rkmE4AMV |
| WITCH-lite: faster query alignment for TIPP3 profiling | cs581/witchlite, claude/cs581-witchlite | session_019bGdVUa8cxTYZFxv5vLeYW |
| Predicted-3Di (ProstT5) protein MSA, and as MAGUS evidence | cs581/prost3di, claude/cs581-prost3di | session_014DeZx1BZSUWuGQZsibVsdX |
| KH-test early stopping for IQ-TREE 3 | cs581/iqstop, claude/cs581-iqstop | session_018wqFifessu4UYH14WK1CxP |
| Identifiability and choosing k in phylogenetic distance deconvolution (theory) | cs581/decodiphy, claude/cs581-decodiphy | session_01D8vRCMw44Lw5ADdqt4qm6W |
| Downstream payoff of the EPA-ng fix: BSCAMPP with larger subtrees, SIMD speed, TIPP3/PICRUSt2 | cs581/epangdown, claude/cs581-epangdown | session_01XUxoJUYQhJFfVzLYXzkjTJ |
| Does the EPA-ng bug change PICRUSt2 placements and predictions? | cs581/picrust, claude/cs581-picrust | session_01YN3K4y9shtvUweCRkCx7DK |
| PASTA rerun on RNASim R0 (-d rna) and 16S.M R0 (-d dna) after the e2e datatype fix | claude/cs581-pastafix | session_01QMhfydVNQZ7wmwthEnGdTM |
| Held-out test of consistency-filtered GCM evidence (MAFFT-only), pre-registered | cs581/protcons, claude/cs581-protcons | session_01EFnDHNVbHHWbWNFzR63MHi |
| Generalize consensus GCM evidence to DNA/RNA (reliable second opinions, soft weights, self-consistency, agreement switch) | cs581/gcmgen, claude/cs581-gcmgen | session_01Mde1BAA2LWYXEvQGVhUdyE |

## Round 2: general directions in parallel (launched 2026-10-10 15:10 UTC)

| idea | dir / branch | session |
|---|---|---|
| Blending disjoint tree merger at scale (GTM-Blend round 2) | cs581/gtmscale, claude/cs581-gtmscale | session_01YDdCNu6WP2RqDSJGtsTMuH |
| ASTRID-Pro: theorem write-up + broader GDL evaluation | cs581/astridpro2, claude/cs581-astridpro2 | session_01SA9eLis44FAe1cbduxbZqp |
| ASTRAL-Pro GDL inconsistency: proof + practical reach | cs581/astralpro2, claude/cs581-astralpro2 | session_01VjGVSAAam7ugjhXoweAbar |
| Distance-mixture identifiability: proofs + identifiability-aware output | cs581/decodiphy2, claude/cs581-decodiphy2 | session_01MHDF53LS457HdYpmfTTir4 |
| One general improved MAGUS (self-soft + consensus evidence) for DNA, RNA and proteins | cs581/magusgen, claude/cs581-magusgen | session_013Ua7rdYjAfh7cX812BxYh4 |
| gcmtrees: does the GCM evidence recipe improve ML trees? (FastTree/IQ-TREE on 16 simulated protein sets) | claude/cs581-gcmtrees | session_012jUiaZ8rzSFs56nuL5D3Qt |
| fragml2: fast ML trees with fragmentary sequences (backbone + fixed EPA-ng / constrained ML) | claude/cs581-fragml2 | session_017GJCTSw6j4rhjwWQvYujTi |
| wapro: weighted ASTRID-Pro (support/length weighting, contraction, missing-data norm) vs best GDL methods | claude/cs581-wapro | session_013LfRmp3zorCpsEbqRSCdPT |
| aproroom: headroom for GDL species-tree methods (true vs estimated trees/tags, 1000 species, fair speed, hybrid) | claude/cs581-aproroom | session_01S5dW3pZ3M1RoZppz3YMxpv |
| gcmclust: different clustering method inside GCM (MCL inflation, Leiden/Louvain, CC, agglomerative), with/without support filter | claude/cs581-gcmclust | session_01NkgLadfKkHYrR3iaEo1CYT |
| basemeth: MAGUS/PASTA with MUSCLE5/FAMSA2/ProbCons/Prank/Clustal subset aligners + Regressive comparison | claude/cs581-basemeth | session_01UZ9gUSWa5AGU8MVZgk2wfJ |
| gcmtrees-h1 helper: SIMHIGH_R6-R8 | claude/cs581-gcmtrees-h1 | session_01CVUHoeHrzN39vhfp85zAWp |
| gcmtrees-h2 helper: SIMMOD_R1-R4 | claude/cs581-gcmtrees-h2 | session_01LrNTfKjroQAetYhP1jYT7T |
| gcmtrees-h3 helper: SIMMOD_R5-R8 | claude/cs581-gcmtrees-h3 | session_01RvNzN9F7tm6AN2BRzhv4AT |
