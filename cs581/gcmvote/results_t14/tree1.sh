# AI-assisted (Claude), exploration code for CS581 project
#!/bin/bash
# one tree: METHOD SRC RESULTS_FILE VARIANT  (waits until VARIANT's row is in RESULTS_FILE)
m=$1; src=$2; res=$3; v=$4
T=/opt/work/gcmvote/trees/SIMHIGH_R14_d0
until grep -qF "\"variant\": \"$v\"" "$res" 2>/dev/null; do sleep 20; done
ln -sf "$src" $T/$m.fasta
nice -n 10 /opt/mm/root/envs/pasta183/bin/python /home/user/scratch/cs581/protbench/code/trees.py $T \
  /opt/data/sim/SIMHIGH/R14/tree.nwk $T/trees.jsonl $m
