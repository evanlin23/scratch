#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Restartable: one fresh MAGUS draw (paper flags, as gcmgen fresh.sh) -> rep dir -> gg baselines + vote variants.
#   bash pipe.sh NAME TRUE_ALIGNMENT DRAW_THREADS
name=$1; src=$2; th=$3
C=/home/user/scratch/cs581/code; W=/opt/work/gcmvote; H=/home/user/scratch/cs581/gcmvote/code_h5b
rep=$W/reps/$name; mkdir -p $W/reps
cd $C
if [ ! -d $rep/sets/s0 ]; then
  echo "$name $src 25" > $W/job_$name.txt
  python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/draws --draws 0 --tools '' --e2e '' --union '' --threads $th || exit 1
  python3 /home/user/scratch/cs581/protcons/code/pc.py rep $W/draws/${name}_d0 $rep || exit 1
fi
python3 /home/user/scratch/cs581/gcmgen/code/gg.py run $rep linsi 'linsi#es3' 'linsi#es4' 'linsi#es5' > $W/gg_$name.log 2>&1 &
python3 $H/par.py $rep 3 magus es4 hard soft soft2 soft4 hard-bb soft-bb
wait
echo PIPE_DONE $name
