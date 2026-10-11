#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# aln.sh NAME : gcmtrees run.sh + gg es4 + vote variants
n=$1; cd /home/user/scratch
bash cs581/gcmtrees/code/run.sh $n
python3 cs581/gcmgen/code/gg.py run /opt/work/gcmtrees/reps/$n 'linsi#es4'
G=/opt/work/gcmtrees/reps/$n; V=/opt/work/gcmvote/reps/$n; mkdir -p $V
for f in inputs true.fasta unaligned.fasta; do [ -e $G/$f ] && ln -sfn $G/$f $V/$f; done
python3 cs581/gcmvote/code/run.py $V magus hard hard-bb soft4
echo "ALNSTAGE_DONE $n"
