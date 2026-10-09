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
| local | jobs_local_large.txt | claude/charming-pasteur-yl6k2v | orchestrating session |

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
