#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# One dataset end to end (helper h9b, reusing h9 scripts): MAGUS draw with 20 backbones, merge-only reps with B = 5, 10, 20
# (backbones 1..B of that draw, same subsets), then every merge variant (MAGUS's merge, edge-support
# fractions 0.2-0.5 of B, the pre-registered vote variants) on LANES parallel lanes. Restartable.
#   bash dataset.sh NAME [LANES]
name=$1; LANES=${2:-4}
W=/opt/work/h9b; R=/home/user/scratch; C=$R/cs581/code; G=$R/cs581/gcmgen/code; H=$R/cs581/gcmvote/code_h9
mkdir -p $W; cd $C
while pgrep -f "bbtool_bench $W/job_$name.txt" >/dev/null; do sleep 30; done  # a draw already running
grep -h "^$name " fanout/jobs_*.txt | head -1 | awk '{print $1, $2, $3}' > $W/job_$name.txt
GCMX_NUM_BACKBONES=20 python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/fresh --draws 0 \
  --tools '' --e2e '' --union '' --threads 4 || exit 1
python3 $R/cs581/protcons/code/pc.py rep $W/fresh/${name}_d0 $W/reps/${name}_B20 || exit 1
for B in 5 10; do
  d=$W/reps/${name}_B$B
  if [ ! -d $d/sets/s0 ]; then
    rm -rf $d; mkdir -p $d/inputs/backbones
    cp -r $W/reps/${name}_B20/inputs/subalignments $d/inputs/
    for i in $(seq 1 $B); do cp $W/reps/${name}_B20/inputs/backbones/backbone_${i}_mafft.txt $d/inputs/backbones/; done
    cp $W/reps/${name}_B20/{unaligned.fasta,magus.json} $d/
    python3 -c "import sys; sys.path.insert(0,'$R/cs581/bbevidence/code'); import bbe; bbe.prep('$d', '$W/reps/${name}_B20/true.fasta')"
  fi
done
tasks=$W/tasks_$name.txt; : > $tasks
for B in 20 10 5; do
  for v in $(for f in 0.2 0.3 0.4 0.5; do python3 -c "import math; print('linsi#es%d' % math.ceil($f*$B-1e-9))"; done | awk '!s[$0]++') linsi; do
    echo "python3 $G/gg.py run $W/reps/${name}_B$B $v" >> $tasks
  done
  for v in magus hard hard-bb soft soft-bb soft2 soft4; do
    echo "python3 $H/vote_run.py $W/reps/${name}_B$B $v" >> $tasks
  done
done
xargs -P $LANES -I{} sh -c '{}' < $tasks
echo DATASET_DONE $name
