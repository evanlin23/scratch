#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Helper h4b, restartable: fresh MAGUS draw (paper flags) -> rep, gg.py baselines, vote.py variants, copy rows.
#   bash h4b.sh NAME   (e.g. 1000S1; reference /opt/data/Datasets/ROSE/NAME/R0/rose.aln.true.fasta)
name=$1; C=/home/user/scratch/cs581; W=/opt/work/gcmvote; rep=$W/reps/$name; vrep=$W/vreps/$name
out=$C/gcmvote/results_h4b; mkdir -p $W/reps $W/vreps $out
cd $C/code
if [ ! -d $rep/sets/s0 ]; then
  echo "$name /opt/data/Datasets/ROSE/$name/R0/rose.aln.true.fasta 25" > $W/job_$name.txt
  python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus_$name.jsonl $W/fresh --draws 0 --tools '' --e2e '' --union '' --threads 4 || exit 1
  python3 $C/protcons/code/pc.py rep $W/fresh/${name}_d0 $rep || exit 1
fi
# vote rows go to a separate rep dir (run.py's results.jsonl rows carry "B", gg.py's do not)
mkdir -p $vrep
for f in $rep/*; do b=$(basename $f); case $b in results.jsonl|variants|vote) ;; *) [ -e $vrep/$b ] || ln -s $f $vrep/$b;; esac; done
G="python3 $C/gcmgen/code/gg.py run $rep"; V="python3 $C/gcmvote/code/run.py $vrep"
# 4 lanes; gg.py locks variants, run.py lanes get disjoint variants
( $G linsi 'linsi#es3' 'linsi#es4' 'linsi#es5'; $V hard-bb ) > $W/lane1_$name.log 2>&1 &
( $G 'linsi#es5' 'linsi#es4' 'linsi#es3' linsi; $V soft-bb ) > $W/lane2_$name.log 2>&1 &
( $V magus es4 hard ) > $W/lane3_$name.log 2>&1 &
( $V soft soft2 soft4 ) > $W/lane4_$name.log 2>&1 &
wait
cp $rep/results.jsonl $out/$name.gg.results.jsonl
cp $vrep/results.jsonl $out/$name.results.jsonl
for d in $vrep/vote/*/; do t=$(basename $d); cp $d/model.json $out/$name.$t.model.json; done
cp $W/magus_$name.jsonl $out/$name.fresh_magus.jsonl 2>/dev/null
echo H4B_DONE $name
