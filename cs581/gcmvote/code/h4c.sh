#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Helper h4c, restartable: one fresh MAGUS draw (paper flags, as gcmgen/fresh.sh) per ROSE set, a replicate dir,
# then gg.py baselines (lane 1) and vote.py variants (lanes 2-3) in parallel while the next draw runs.
#   bash h4c.sh "1000M4 /path/true.fasta" "1000S3 /path/true.fasta"
C=/home/user/scratch/cs581/code; W=/opt/work/gcmvote; G=/home/user/scratch/cs581
mkdir -p $W/reps $W/logs; cd $C
for job in "$@"; do
  set -- $job; name=$1; src=$2; rep=$W/reps/$name
  if [ ! -d $rep/sets/s0 ]; then
    echo "$name $src 25" > $W/job_$name.txt
    python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/draws --draws 0 --tools '' --e2e '' --union '' --threads 4 || continue
    python3 $G/protcons/code/pc.py rep $W/draws/${name}_d0 $rep
  fi
  # vote rows go to a twin dir (symlinks to the same inputs) so the two results.jsonl formats stay apart
  v=$W/reps/${name}_vote; mkdir -p $v
  for f in inputs true.fasta sets; do [ -e $v/$f ] || ln -s $rep/$f $v/$f; done
  ( python3 $G/gcmgen/code/gg.py run $rep linsi 'linsi#es3' 'linsi#es4' 'linsi#es5' > $W/logs/$name.gg.log 2>&1; echo GG_DONE >> $W/logs/$name.gg.log ) &
  ( python3 $G/gcmvote/code/run.py $v magus es4 hard soft > $W/logs/$name.vote1.log 2>&1; echo VOTE_DONE >> $W/logs/$name.vote1.log ) &
  ( python3 $G/gcmvote/code/run.py $v hard-bb soft-bb soft2 soft4 > $W/logs/$name.vote2.log 2>&1; echo VOTE_DONE >> $W/logs/$name.vote2.log ) &
done
wait
echo H4C_DONE
