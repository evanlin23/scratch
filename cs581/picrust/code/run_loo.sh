#!/bin/bash
# run_loo.sh <setdir> <stock|fix> : place held-out V4 amplicons into the pruned reference with
# PICRUSt2's place_seqs.py (default options) and predict KOs + NSTI with hsp.py (default mp).
S=$1; v=$2; E=/opt/mm/root/envs/picrust2
export PATH=/opt/epa/$v:$E/bin:$PATH EPA_LOG=$S/epa_times.log
cd $S; rm -rf int_$v placed_$v.tre ko_$v.tsv.gz
s=$(date +%s)
place_seqs.py -s queries.fna -r ref -o placed_$v.tre -p 4 --intermediate int_$v > place_$v.log 2>&1
m=$(date +%s)
hsp.py -t placed_$v.tre -o ko_$v.tsv.gz --observed_trait_table ko_reduced.txt.gz -n -p 4 > hsp_$v.log 2>&1
echo -e "$S\t$v\tplace=$(( m-s ))\thsp=$(( $(date +%s)-m ))" >> /opt/work/loo/times.tsv
