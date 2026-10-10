#!/bin/bash
# Fresh MAGUS draws (paper flags) for sets without cached inputs, then a gg replicate. Restartable.
W=/opt/work/gcmgen; G=/home/user/scratch/cs581/gcmgen/code
cd /home/user/scratch/cs581/code
while read -r name src k; do
  [ -z "$name" ] && continue
  echo "$name $src $k" > $W/job_$name.txt
  python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/fresh --draws 0 --tools '' --e2e '' --union '' --threads 4 || continue
  python3 /home/user/scratch/cs581/protcons/code/pc.py rep $W/fresh/${name}_d0 $W/reps/$name
done < "$1"
cp $W/magus.jsonl /home/user/scratch/cs581/gcmgen/results/fresh_magus.jsonl
echo FRESH_DONE
