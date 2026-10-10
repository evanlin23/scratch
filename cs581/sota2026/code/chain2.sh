#!/bin/bash
# Sequential tail of the benchmark: waits for 1000S3 from chain 1, stops chain 1 before its
# clustalo step (dropping any partial clustalo rows), then runs the remaining jobs in priority order.
cd "$(dirname "$0")"
R=../results/full.jsonl
until grep -q '"1000S3_R0", "tool": "mafft-auto"' $R; do sleep 15; done
kill 21230 2>/dev/null
sleep 1
pkill -f "jobs_clustalo.txt 900"
pkill -x clustalo
sleep 3
python3 drop_rows.py "$R" clustalo
bash run_full.sh twilight-truetree,famsa-truetree jobs_diag.txt 900
bash run_full.sh twilight jobs_16st.txt 2400
bash run_full.sh twilight-1,mafft-parttree,twilight jobs_10k.txt 2400
bash run_full.sh mafft-linsi,muscle5 jobs_bba.txt 900
bash run_full.sh clustalo jobs_clustalo.txt 600
echo CHAIN2 DONE
