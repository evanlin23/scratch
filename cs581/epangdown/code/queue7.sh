#!/bin/bash
# Exp.-5 analogue on a second dataset: RNASim 8K-leaf random sub-backbone (10K does not fit at 1,624 sites).
C=/home/user/scratch/cs581/epangdown/code; P=/opt/mm/root/envs/place/bin
while ! grep -q QUEUE6DONE /opt/work/queue6.log; do sleep 20; done
[ -s /opt/work/rna_8k/backbone.fa ] || $P/python $C/make_subset_dataset.py /opt/work/rna50k 8000 /opt/work/rna_8k 1
$C/run_whole.sh /opt/work/rna_8k "stock fix" frag
$C/run_bscampp.sh /opt/work/rna_8k stock 2000 frag
echo QUEUE7DONE
