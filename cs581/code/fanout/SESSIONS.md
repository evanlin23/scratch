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
| fragscale: two-step ML vs RAxML-NG at equal time, RNASim 10K fragmentary, estimated alignments, stock vs patched EPA-ng | claude/cs581-fragscale | session_01XZ4oCa5kuuTD71UBphf7yY |
| magusfast: profile MAGUS; fewer/smaller backbones with support pruning; >=2x at equal accuracy | claude/cs581-magusfast | session_01ALcCHpGgfjs3X2Uk2qc2kr |
| gcmvote: reference-free vote (binomial-mixture) model for GCM edges; threshold x backbone-count sensitivity; trees SIMHIGH R1-R4 | claude/cs581-gcmvote | session_01HVdpQL7fYBwJKGSwChipPv |
| gcmvote-h1 helper: SIMHIGH R5-R8 MAGUS runs + trees | claude/cs581-gcmvote-h1 | session_01CmMHK7h7WnNYfvGU2r1WB4 |
| gcmvote-h2 helper: SIMHIGH R9-R12 MAGUS runs + trees | claude/cs581-gcmvote-h2 | session_01GoHVT6GpRSf3DY6H3qTo2H |
| gcmvote-h3: held-out BAliBASE BBA0154/0190/0081/0117 + HomFam 4 | claude/cs581-gcmvote-h3 | session_01SrCc4mFPMZLbGNWLqfWBw7 |
| gcmvote-h4: held-out ROSE 1000L3/M3/S1/S2/M4/S3 | claude/cs581-gcmvote-h4 | session_01UKXEnuJh7YdekbkauWYR2B |
| gcmvote-h5: held-out SIMMOD_R2, SIMHIGH_R2, 1000M2_R1, 1000L1_R1, RNASim R0/R1 | claude/cs581-gcmvote-h5 | session_01KKVjwN8B5Zvb1erTmsSVPA |
| gcmvote-h6: B=5/10/20 sensitivity, proteins | claude/cs581-gcmvote-h6 | session_017cEpz8mp13zBM3tp1rkXBS |
| gcmvote-h9: B=5/10/20 sensitivity, DNA/RNA | claude/cs581-gcmvote-h9 | session_01VDyiALio6XXjAwodF381Vu |
| gcmvote-h7: trees SIMHIGH R13-R16 | claude/cs581-gcmvote-h7 | session_017NndPrbeyffzkvavo85KHn |
| gcmvote-h8: trees SIMHIGH R17-R20 | claude/cs581-gcmvote-h8 | session_01KnTMhP5zb9H1gTX9otMqa2 |
| gcmvote-h10: MAGUS-paper ROSE 1000M1 R0-R3 (aln + trees) | claude/cs581-gcmvote-h10 | session_0183UYDvTv7komXdMwUC7GQ6 |
| gcmvote-h11: MAGUS-paper 16S.3 | claude/cs581-gcmvote-h11 | session_01VTnVnxuXg4odXiiianLbP1 |
| gcmvote-h12: MAGUS-paper 16S.T | claude/cs581-gcmvote-h12 | session_012UsuJKCUwQfEwd4ZNw1RsK |
| gcmvote-h13: MAGUS-paper RNASim 10K (aln + trees) | claude/cs581-gcmvote-h13 | session_01XNEer5rT7xkcSZPzDWcfPG |

### gcmvote fan-out (02:27–02:33 UTC, Oct 11): one dataset or replicate per machine

Existing helpers were cut down to the unit already running: h1 SIMHIGH_R5, h2 R9, h7 R13, h8 R17, h3 BBA0154/0190, h4 1000L3/M3,
h5 SIMMOD_R2/SIMHIGH_R2, h6 BBA0101/0067, h9 1000M2/1000L1, h10 1000M1 R0/R1. Main dropped the SIMMOD_R2 and SIMHIGH_R3/R4 draws.

| session | branch | session id |
|---|---|---|
| gcmvote-t03: SIMHIGH_R3 trees | claude/cs581-gcmvote-t03 | session_01UdLszXi3YRBJB5MsHefVQ6 |
| gcmvote-t04: SIMHIGH_R4 trees | claude/cs581-gcmvote-t04 | session_01AvPG43SGyEiChAUUej5ehA |
| gcmvote-t06: SIMHIGH_R6 trees | claude/cs581-gcmvote-t06 | session_019ZwJtiLC2S7rUetoD9XStz |
| gcmvote-t07: SIMHIGH_R7 trees | claude/cs581-gcmvote-t07 | session_01RAHZdFFx8yVj3XwYhD6Xa3 |
| gcmvote-t08: SIMHIGH_R8 trees | claude/cs581-gcmvote-t08 | session_01XhSC4usqEwJVQx86TFraMP |
| gcmvote-t10: SIMHIGH_R10 trees | claude/cs581-gcmvote-t10 | session_01MiwjsURTu5gX3EQYDGvWCC |
| gcmvote-t11: SIMHIGH_R11 trees | claude/cs581-gcmvote-t11 | session_01PWmPsqizmpvkd5nETirMEh |
| gcmvote-t12: SIMHIGH_R12 trees | claude/cs581-gcmvote-t12 | session_017Jwkz9heH9WUPFbFoUZNUV |
| gcmvote-t14: SIMHIGH_R14 trees | claude/cs581-gcmvote-t14 | session_01VbuGGfY4igvE3DjEcRECC6 |
| gcmvote-t15: SIMHIGH_R15 trees | claude/cs581-gcmvote-t15 | session_01GdEV8KfY1cs2MQ66xMojpb |
| gcmvote-t16: SIMHIGH_R16 trees | claude/cs581-gcmvote-t16 | session_015kjcJVfCE4D8LTLeS51mn6 |
| gcmvote-t18: SIMHIGH_R18 trees | claude/cs581-gcmvote-t18 | session_011sDrsbPqFgwo4fXGbr7rK9 |
| gcmvote-t19: SIMHIGH_R19 trees | claude/cs581-gcmvote-t19 | session_0142FFeMDqHi17MRNZFEqBHP |
| gcmvote-t20: SIMHIGH_R20 trees | claude/cs581-gcmvote-t20 | session_01Vu6yW4PFqYYxozZN1zwbnM |
| gcmvote-h3b: BBA0081 + BBA0117 | claude/cs581-gcmvote-h3b | session_01XxU57Ta4e6gMVwVa4ELGdP |
| gcmvote-h3c: HomFam Acetyltransf + PDZ | claude/cs581-gcmvote-h3c | session_01EUax5HQBk5gemvbDF8zivF |
| gcmvote-h3d: HomFam blmb + aat | claude/cs581-gcmvote-h3d | session_01Pd4qkoCN5ekJuPFR7GtrW2 |
| gcmvote-h4b: ROSE 1000S1 + 1000S2 | claude/cs581-gcmvote-h4b | session_01X79TL2RNLbxqWSAqBRsFFX |
| gcmvote-h4c: ROSE 1000M4 + 1000S3 | claude/cs581-gcmvote-h4c | session_01Kus3fzHct4cBXL9LDQbS7N |
| gcmvote-h5b: 1000M2_R1 + 1000L1_R1 | claude/cs581-gcmvote-h5b | session_016CZioNZLM2aK8TNKi2pVJE |
| gcmvote-h5c: RNASim 1000 R0 | claude/cs581-gcmvote-h5c | session_01DewgomzXpRiVTssysJJTAB |
| gcmvote-h5d: RNASim 1000 R1 | claude/cs581-gcmvote-h5d | session_01WzUVeaJp2iiBBzBq7Lko3R |
| gcmvote-h6b: B=5/10/20, SIMHIGH_R1 + SIMMOD_R1 | claude/cs581-gcmvote-h6b | session_019WaEjgwQKzC1XMbVVAGCM5 |
| gcmvote-h9b: B=5/10/20, 16S.M | claude/cs581-gcmvote-h9b | session_01HmXDzVoAJo7yo2KqBqR9dQ |
| gcmvote-h10b: 1000M1 R2 (aln + trees) | claude/cs581-gcmvote-h10b | session_01RAcTHCaeXkpEmGH4u7ehmm |
| gcmvote-h10c: 1000M1 R3 (aln + trees) | claude/cs581-gcmvote-h10c | session_01WzqN86jpzzUL1Y9hYdqhNs |

### Round 3 (04:07 UTC, Oct 11): open questions the student raised

| session | branch | session id |
|---|---|---|
| gcmwhy: why filtering helps proteins but not DNA/RNA (headroom, near-miss votes, difficulty, MCL/trace) | claude/cs581-gcmwhy | session_016o1oVimWEtaon5tqZPMKHy |
| gcmvote2: better voting models (strength-weighted, Dawid–Skene, overlap-aware, combined), pre-registered | claude/cs581-gcmvote2 | session_01UrxbGvxw26RRQwfdbbwTT7 |
| bbsize: backbone size (100/200/400) and count (5–40) × filtering | claude/cs581-bbsize | session_01HKPtCJnVFneBWA5FthQfxZ |

All gcmvote helpers were asked (03:58–04:05 UTC) to push a "rep bank" (cs581/gcmvote/bank/<rep>.tar.gz + MANIFEST.tsv) so new
merge-step variants can be tested without re-running MAGUS.

### Round 4 (04:14–04:18 UTC, Oct 11): "explore everything"

| session | branch | session id |
|---|---|---|
| gcmvote-p21: tree power SIMHIGH_R21+R22 | claude/cs581-gcmvote-p21 | session_01Qr2UcPghiNLnqRbRxPSPzJ |
| gcmvote-p23: R23+R24 | claude/cs581-gcmvote-p23 | session_01GE7ZUFGcGMPpe3Wgjv7B5r |
| gcmvote-p25: R25+R26 | claude/cs581-gcmvote-p25 | session_015jpoNBBfrrDwyZoMCGmH75 |
| gcmvote-p27: R27+R28 | claude/cs581-gcmvote-p27 | session_01GhawMA68QavRQKdYN99EUg |
| gcmvote-p29: R29+R30 | claude/cs581-gcmvote-p29 | session_011xHKkFnFr4KVNw9ZCcnXuD |
| gcmvote-p31: R31+R32 | claude/cs581-gcmvote-p31 | session_01JpHuMWw7NzUfcnuUEVNLUy |
| gcmvote-p33: R33+R34 | claude/cs581-gcmvote-p33 | session_014RRLjAjqtFpcPRscdSxZV6 |
| gcmvote-p35: R35+R36 | claude/cs581-gcmvote-p35 | session_01McQ9xFQ4VboVvMuGTK9ub3 |
| gcmvote-p37: R37+R38 | claude/cs581-gcmvote-p37 | session_01VPw5GucbxB3mwZCgcovxT4 |
| gcmvote-p39: R39+R40 | claude/cs581-gcmvote-p39 | session_01TLcyu34FPWXn6wX3w6o7tn |
| gcmvote-p41: R41+R42 | claude/cs581-gcmvote-p41 | session_01PZxSskpymjxgAohzZKAVZD |
| gcmvote-p43: R43+R44 | claude/cs581-gcmvote-p43 | session_018hFbZBoC2CGFoX9TRs6Yfz |
| gcmvote-p45: R45+R46 | claude/cs581-gcmvote-p45 | session_01FmXx5t1KNYq7Gvc9tR8HXC |
| gcmvote-p47: R47+R48 | claude/cs581-gcmvote-p47 | session_01RdobQvx2BtrKWbkDRPjdBe |
| gcmvote-p49: R49+R50 | claude/cs581-gcmvote-p49 | session_01U51E5YWrT9m5dVhR1SyZoR |
| treecrit: which alignment errors predict tree error + one tree-targeted fix | claude/cs581-treecrit | session_01L2ecceKe8XyMBbwKuYDjWe |

Also: gcmvote2 asked to add M5 (reference-free dataset gate) and M6 (CPM clustering combo); local novelty audit of the voting
ideas → cs581/literature/novelty_gcmvote.md.
