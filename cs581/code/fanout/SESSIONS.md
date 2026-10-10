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
