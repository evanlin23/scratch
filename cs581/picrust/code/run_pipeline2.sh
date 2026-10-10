#!/bin/bash
# run_pipeline.sh <dataset> <stock|fix|conda>
# Full PICRUSt2 2.6.3 default pipeline (picrust2_pipeline.py, default options) with the chosen EPA-ng on PATH.
ds=$1; v=$2; D=/home/user/picrust/ms/data/16S_datasets/$ds
E=/opt/mm/root/envs/picrust2
export PATH=/opt/epa/$v:$E/bin:$PATH EPA_LOG=/opt/work/epa_times.log
fa=$(ls $D/*rep_seqs.f* $D/*asvs.fna 2>/dev/null | head -1)
tab=/opt/work/$ds/table.tsv   # biom-converted table with the "# Constructed from biom file" line removed
out=/opt/work/$ds/pipe_$v; rm -rf $out
s=$(date +%s)
/home/user/scratch/cs581/picrust/code/picrust2_pipeline_nocheck.py -s $fa -i $tab -o $out -p 4 --verbose > /opt/work/$ds/pipe_$v.log 2>&1
rc=$?
echo -e "$ds\t$v\t$(( $(date +%s)-s ))\t$rc" >> /opt/work/pipe_times.tsv
