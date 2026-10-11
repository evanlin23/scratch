#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Helper h3c, restartable: HomFam NAME -> fresh MAGUS draw (paper flags, as gcmgen fresh.sh), rep dir,
# gg.py baselines, vote.py variants (run.py; scored vs the Homstrad seeds = REP/true.fasta).
#   bash h3c.sh NAME THREADS     e.g. bash h3c.sh Acetyltransf 2
name=$1; th=${2:-4}; W=/opt/work/gcmvote_h3c; R=/home/user/scratch/cs581
rep=$W/reps/HF_$name; mkdir -p $W
cd $R/code
if [ ! -d $rep/sets/s0 ]; then
  echo "HF_$name /opt/data/homfam2k/$name/ref.fasta 25 /opt/data/homfam2k/$name/unaln.fasta" > $W/job_$name.txt
  python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus_$name.jsonl $W/fresh --draws 0 --tools '' --e2e '' --union '' --threads $th || exit 1
  python3 $R/protcons/code/pc.py rep $W/fresh/HF_${name}_d0 $rep || exit 1
fi
python3 $R/gcmgen/code/gg.py run $rep linsi 'linsi#es3' 'linsi#es4' 'linsi#es5' &
python3 $R/gcmvote/code/run.py $rep magus es4 hard &
python3 $R/gcmvote/code/run.py $rep soft soft2 soft4 &
python3 $R/gcmvote/code/run.py $rep hard-bb soft-bb &
wait
echo H3C_DONE $name
