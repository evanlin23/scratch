You are an unattended compute worker for a CS581 research project (multiple
sequence alignment; reproducing the MAGUS paper's benchmarks and testing
merge-step variants). The repository is already checked out on branch
`claude/charming-pasteur-yl6k2v`. Do exactly the following and nothing else:
do not modify any code or job files, do not open pull requests, and do not
push to any branch other than `BRANCH`.

1. From the repository root run `bash cs581/code/setup.sh` (installs aligners
   and downloads ~390 MB of benchmark data; takes ~5-10 minutes). If a step
   fails because of a transient network error, run it again.

2. Run this with the Bash tool using `run_in_background: true` and
   `timeout: 7200000`:

       bash cs581/code/fanout/worker.sh cs581/code/fanout/JOBFILE BRANCH

   It runs MAGUS on benchmark replicates (about 30-40 minutes per job, ~10 jobs)
   and commits + pushes small result files to `BRANCH` after every job.

3. Do not poll and do not use foreground `sleep`; you are notified when the
   background command stops. When it stops, read the last lines of its output:
   - if they contain `ALL JOBS DONE`, you are finished;
   - otherwise (for example the 2-hour background limit was reached), run the
     exact same command again in the background. It skips finished jobs and
     cleans up after the interrupted one. Repeat until `ALL JOBS DONE`.

4. Finish with a short message listing which job names completed and which
   (if any) printed `FAILED`.
