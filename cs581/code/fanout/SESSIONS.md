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
