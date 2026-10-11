#!/bin/bash
# Orchestrator follow-ups: whole-backbone EPA-ng (10K-leaf random nt78 sub-backbone; the 77K one
# needs ~100 GB), BSCAMPP on the same 10K backbone, and BSCAMPP(p) (pplacer) on nt78.
C=/home/user/scratch/cs581/epangdown/code; P=/opt/mm/root/envs/place/bin
[ -s /opt/work/nt78_10k/backbone.fa ] || $P/python $C/make_subset_dataset.py /opt/work/nt78 10000 /opt/work/nt78_10k 1
$C/run_whole.sh /opt/work/nt78_10k "stock fix" frag
$C/run_bscampp.sh /opt/work/nt78_10k stock 2000 frag
$C/run_bscampp.sh /opt/work/nt78_10k fix 5000 frag
$C/run_bscampp.sh /opt/work/nt78 pplacer 2000 frag
$C/run_bscampp.sh /opt/work/nt78 pplacer 5000 frag
echo QUEUE6DONE
