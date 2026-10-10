#!/bin/bash
# Restartable driver: rerun after an interruption; bbtool_bench skips finished (dataset, draw) rows.
cd /home/user/scratch/cs581/code
P=/home/user/scratch/cs581/protbench
python3 -m gcmx.bbtool_bench $P/jobs.txt $P/results/bbtool.jsonl /opt/work/protbench --draws 0 --tools clustalo --union '' --e2e clustalo --threads 4
