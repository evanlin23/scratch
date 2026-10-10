#!/bin/bash
# Create MAGUS inputs (subsets + backbones) for BBA0117, as bbtool_bench does, then run the 3Di merges.
cd /home/user/scratch/cs581/code
W=/opt/p3d/rv100runs/BBA0117; M=$W/magus_e2e
rm -rf $M $W/magus_e2e.fasta
python3 -m gcmx.run_magus --gcmx-fastgraph false -np 4 -d $M -i /opt/p3d/rv100/BBA0117.fa -o $W/magus_e2e.fasta \
  --maxsubsetsize 0 --maxnumsubsets 25 --decompstrategy pastastyle --decompskeletonsize 300 --graphbuildmethod mafft \
  --graphclustermethod mcl --graphtracemethod minclusters --graphtraceoptimize false -r 10 -m 200 -f 4 > $W/magus_e2e.log 2>&1
mkdir -p $W/inputs/backbones && cp -r $M/subalignments $W/inputs/ && cp $M/graph/backbone_*_mafft.txt $W/inputs/backbones/
python3 -m gcmx.score /home/user/scratch/cs581/data/balibase_clean/RV100_BBA0117.fasta $W/magus_e2e.fasta > $W/magus_e2e.score.json
rm -rf $M
python3 ../prost3di/code/rv100.py BBA0117 > /opt/p3d/rv100_BBA0117b.log 2>&1
