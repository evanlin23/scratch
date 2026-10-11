#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Restartable: HomFam held-out sets (protbench subsamples, 2,000 seqs, seed 1, scored on the Homstrad seeds).
# One fresh MAGUS draw (paper flags) -> rep dir (pc.py rep) -> gg.py L-INS-i baselines and vote.py variants.
# vote.py runs on a symlinked copy of the rep (VREP) because run.py and gg.py both append to REP/results.jsonl.
#   bash h3d.sh blmb aat
C=/home/user/scratch/cs581; W=/opt/work/gcmvote; R=$C/gcmvote/results_h3d
mkdir -p $W/reps $W/vreps $R
cd $C/code
for fam in "$@"; do
  name=HF_$fam; rep=$W/reps/$name; vrep=$W/vreps/$name
  if [ ! -d $rep/sets/s0 ]; then
    echo "$name /opt/data/homfam2k/$fam/ref.fasta 25 /opt/data/homfam2k/$fam/unaln.fasta" > $W/job_$name.txt
    python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W/draws --draws 0 --tools '' --e2e '' --union '' --threads 4 || continue
    python3 $C/protcons/code/pc.py rep $W/draws/${name}_d0 $rep || continue
  fi
  if [ ! -d $vrep ]; then
    mkdir -p $vrep.tmp
    for f in inputs true.fasta unaligned.fasta magus.json sets; do ln -s $rep/$f $vrep.tmp/$f; done
    mv $vrep.tmp $vrep
  fi
  python3 $C/gcmgen/code/gg.py run $rep linsi 'linsi#es3' 'linsi#es4' 'linsi#es5' > $W/gg_$name.log 2>&1 &
  python3 $C/gcmvote/code/run.py $vrep magus hard soft2 hard-bb > $W/vote1_$name.log 2>&1 &
  python3 $C/gcmvote/code/run.py $vrep es4 soft soft4 soft-bb > $W/vote2_$name.log 2>&1 &
  wait
  cp $rep/results.jsonl $R/$name.gg.results.jsonl
  cp $vrep/results.jsonl $R/$name.results.jsonl
  cp $rep/magus.json $R/$name.magus.json
  grep "\"$name\"" $W/magus.jsonl > $R/$name.bbtool.jsonl
  for d in $vrep/vote/*/; do cp $d/model.json $R/$name.$(basename $d).model.json; done
  echo "DONE $name"
done
echo H3D_DONE
