#!/bin/bash
# re-run hsp for the leave-out set (after dropping all-zero KO columns), score it, then rerun hmp
E=/opt/mm/root/envs/picrust2/bin; S=/opt/work/loo/s1; C=/home/user/scratch/cs581/picrust/code
cd $S
for v in stock fix; do
  s=$(date +%s); $E/hsp.py -t placed_$v.tre -o ko_$v.tsv.gz --observed_trait_table ko_reduced.txt.gz -n -p 4 > hsp_$v.log 2>&1
  echo -e "$S\t$v\thsp_rerun=$(( $(date +%s)-s ))" >> /opt/work/loo/times.tsv
done
$E/python $C/score_loo.py $S > $S/summary.json 2> $S/score.err
for v in stock fix; do $C/run_pipeline2.sh hmp $v; done
