#!/bin/bash
# b=10000 runs alone (EPA-ng needs ~12 GB at 10,000 tips x 1,286 sites; did not fit beside FastTree).
C=/home/user/scratch/cs581/epangdown/code
$C/run_bscampp.sh /opt/work/nt78 "fix stock" 10000 frag
echo QUEUE2DONE
