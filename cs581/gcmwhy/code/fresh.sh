#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Fresh MAGUS draws (paper flags) for training reps missing from the bank (BBA0134, BBA0039, 1000L2), as gcmgen fresh.sh.
W=/tmp/claude-0/fresh; mkdir -p $W/reps; C=/home/user/scratch/cs581/code
cd $C
while read -r name src k; do
  rep=$W/reps/$name; [ -d $rep/inputs ] && { [ -e /tmp/claude-0/bank/$name ] || ln -s $rep /tmp/claude-0/bank/$name; continue; }
  echo "$name $src $k" > $W/job_$name.txt
  nice -n 5 python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/draws --draws 0 --tools '' --e2e '' --union '' --threads 4 || continue
  python3 /home/user/scratch/cs581/protcons/code/pc.py rep $W/draws/${name}_d0 $rep && ln -s $rep /tmp/claude-0/bank/$name
done <<JOBS
BBA0134 /home/user/scratch/cs581/data/balibase_clean/RV100_BBA0134.fasta 25
BBA0039 /home/user/scratch/cs581/data/balibase_clean/RV100_BBA0039.fasta 25
1000L2 /opt/data/Datasets/ROSE/1000L2/R0/rose.aln.true.fasta 25
JOBS
echo FRESH_DONE
