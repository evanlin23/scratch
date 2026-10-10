#!/bin/bash
# After queue2: determinism/noise checks on 16S, RNASim 50K runs, standalone speed test.
C=/home/user/scratch/cs581/epangdown/code
while ! grep -q QUEUE2DONE /opt/work/queue2.log; do sleep 20; done
TAG=_rep2 $C/run_bscampp.sh /opt/work/16S stock 2000 frag          # unseeded replicate: noise floor
export PYTHONHASHSEED=0
TAG=_s0 $C/run_bscampp.sh /opt/work/16S "stock fix" "2000 5000" frag
TAG=_s0 $C/run_bscampp.sh /opt/work/rna50k "stock fix" "2000 5000" frag
TAG=_s0 $C/run_bscampp.sh /opt/work/nt78 "stock" 2000 frag
TAG=_s0b $C/run_bscampp.sh /opt/work/nt78 "stock" 2000 frag   # same seed twice: is it deterministic?
unset PYTHONHASHSEED
$C/speedtest.sh /opt/work/nt78 5000 200 1
$C/speedtest.sh /opt/work/nt78 10000 200 1
echo QUEUE3DONE
