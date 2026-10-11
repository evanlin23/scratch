#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# One SIMHIGH replicate: vote.py merges (magus, es4, hard-bb), FastSP check of magus vs the bank's
# recorded SP error, then IQ-TREE on true/magus/es4/hard-bb, 2 runs at a time. Restartable.
#   bash run_rep.sh SIMHIGH_R5      (rep unpacked in /opt/work/iqc/reps/<REP>)
set -euo pipefail
R=/home/user/scratch; W=/opt/work/iqc; name=$1; rep=$W/reps/$name
OUT=$R/cs581/iqtree/results_c/iqtree.jsonl; PY=/opt/mm/root/envs/pasta183/bin/python
for v in magus es4 hard-bb; do
  [ -s $rep/vote/$v/run.json ] || python3 $R/cs581/gcmvote/code/vote.py $rep $v --threads 4
done
if [ ! -s $rep/vote/magus/fastsp.txt ]; then
  java -Xmx6g -jar /opt/tools/FastSP/FastSP.jar -r $rep/true.fasta -e $rep/vote/magus/out.fasta > $rep/vote/magus/fastsp.txt 2>/dev/null
fi
echo "magus rebuilt:"; grep -E "^SP-Score|^SPFN|^SPFP" $rep/vote/magus/fastsp.txt
echo "magus recorded:"; python3 -c "import json;d=json.load(open('$rep/magus.json'))['magus'];print('SPFN',d['SPFN'],'SPFP',d['SPFP'])"
run() { $PY $R/cs581/iqtree/code_c/iqtree_one.py $name $1 $2 $rep/true_tree.nwk $W/iq/$name $OUT; }
run true $rep/true.fasta & run magus $rep/vote/magus/out.fasta & wait
run es4 $rep/vote/es4/out.fasta & run hard-bb $rep/vote/hard-bb/out.fasta & wait
echo "REP DONE $name"
