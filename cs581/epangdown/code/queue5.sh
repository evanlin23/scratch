#!/bin/bash
# PICRUSt2 controls: stock rerun on the same subset (noise), and a second 12K-tip subset (seed 2).
C=/home/user/scratch/cs581/epangdown/code; P=/opt/mm/root/envs/place/bin
mkdir -p /opt/work/pic_rerun /opt/work/pic_s2
cd /opt/work/pic_rerun && ln -sf /opt/work/pic/sub_ref sub_ref && cp /opt/work/pic/asv.fna /opt/work/pic/table.biom .
sed 's|for v in stock fix; do|for v in stock; do|' $C/run_picrust2_subref.sh > /tmp/claude-0/sp/run_stock_only.sh && bash /tmp/claude-0/sp/run_stock_only.sh /opt/work/pic_rerun 4
cd /opt/work/pic_s2 && cp /opt/work/pic/asv.fna /opt/work/pic/table.biom . && $P/python $C/picrust2_subref.py 12000 /opt/work/pic_s2/sub_ref 2 > subref.log 2>&1
$C/run_picrust2_subref.sh /opt/work/pic_s2 4
echo QUEUE5DONE
