#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Helper h5: every pre-declared vote variant (PREREG.md; hard+mask skipped, its alignment equals hard) on a rep
# built by run_h5.sh, scored by gcmvote/code/run.py. Vote rows go to a separate view of the rep (symlinks) so
# gg.py's results.jsonl is untouched. 3 parallel lanes. Restartable.
#   bash vote_h5.sh NAME
W=/opt/work/h5; V=/home/user/scratch/cs581/gcmvote/code; name=$1
src=$W/reps/$name; rep=$W/vreps/$name
[ -d $src/sets/s0 ] || { echo "no rep $name"; exit 1; }
mkdir -p $rep
for f in $src/*; do b=$(basename $f); case $b in results.jsonl|variants|vote) ;; *) [ -e $rep/$b ] || ln -s $f $rep/$b;; esac; done
python3 $V/run.py $rep magus es4 hard frac0.2 > $W/vlane1_$name.log 2>&1 &
python3 $V/run.py $rep hard-bb soft soft-bb frac0.3 > $W/vlane2_$name.log 2>&1 &
python3 $V/run.py $rep soft2 soft4 frac0.4 frac0.5 > $W/vlane3_$name.log 2>&1 &
wait
python3 $V/calib.py $rep magus > $W/calib_$name.log 2>&1
grep -h FAILED $W/vlane?_$name.log
echo "VOTE_DONE $name"
