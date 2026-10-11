#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# h4: fresh MAGUS draw (paper flags) + rep dir (as gcmgen/code/fresh.sh) + gg.py baselines, per dataset. Restartable.
W=/opt/work/h4; R=/home/user/scratch/cs581
mkdir -p $W/reps
cd $R/code
while read -r name src k; do
  [ -z "$name" ] && continue
  if [ ! -d $W/reps/$name/sets/s0 ]; then
    echo "$name $src $k" > $W/job_$name.txt
    python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/fresh --draws 0 --tools '' --e2e '' --union '' --threads 4 || continue
    python3 $R/protcons/code/pc.py rep $W/fresh/${name}_d0 $W/reps/$name || continue
  fi
  python3 $R/gcmgen/code/gg.py run $W/reps/$name linsi 'linsi#es3' 'linsi#es4' 'linsi#es5'
  echo "BASE_DONE $name"
done < "$1"
echo ALL_BASE_DONE
