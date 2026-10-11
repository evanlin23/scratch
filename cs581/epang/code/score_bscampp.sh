#!/bin/bash
# Score all finished BSCAMPP runs of one replicate: score_bscampp.sh <datadir>
D=$1; args=()
for d in $D/bscampp/*/; do n=$(basename $d); [ -s $d/result.jplace ] && args+=("$n=$d/result.jplace"); done
rm -f $D/bscampp/scores.tsv
/opt/mm/root/envs/place/bin/python -W ignore /home/user/scratch/cs581/epang/code/score.py $D $D/bscampp/scores.tsv "${args[@]}"
