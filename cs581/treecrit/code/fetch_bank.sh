#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# Extract every bank replicate in reps.tsv to /opt/work/treecrit/reps/<key> (restartable).
set -u
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../../.." && pwd); W=/opt/work/treecrit
mkdir -p $W/reps $W/bank
grep -v '^#' $HERE/reps.tsv | while IFS=$'\t' read key br path res src; do
  [ -f $W/reps/$key/true.fasta ] && continue
  git -C $REPO fetch -q origin $br 2>/dev/null
  git -C $REPO show origin/$br:$path > $W/bank/$key.tgz || { echo "no bank $key"; continue; }
  d=$W/bank/x_$key; rm -rf $d; mkdir -p $d; tar xzf $W/bank/$key.tgz -C $d
  sub=$(ls $d); rm -rf $W/reps/$key; mv $d/$sub $W/reps/$key; rmdir $d
  echo "ok $key $(ls $W/reps/$key | tr '\n' ' ')"
done
