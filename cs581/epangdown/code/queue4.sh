#!/bin/bash
# Standalone speed test, then PICRUSt2 on a reduced reference (each alone on the machine).
C=/home/user/scratch/cs581/epangdown/code
$C/speedtest.sh /opt/work/nt78 5000 200 1
$C/speedtest.sh /opt/work/nt78 10000 200 1
$C/run_picrust2_subref.sh /opt/work/pic 4
echo QUEUE4DONE
