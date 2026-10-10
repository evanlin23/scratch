#!/bin/bash
# FastTree -lg -gamma (single-threaded, niced) on true / MAGUS / recipe / es3 / hard-filter alignments of every
# replicate whose variants are all done; nRF vs the true tree via protbench/code/trees.py. Restartable.
#   bash trees.sh [PARALLEL]     one pass;   LOOP=1 bash trees.sh  keeps scanning until $W/STOP_TREES exists
W=/opt/work/gcmtrees; G=/home/user/scratch/cs581/gcmtrees
one() {
  rep=$1; n=$(basename $rep); T=$W/trees/${n}_d0
  for v in linsi wsoft0.03_c_linsi_i_fftns2_es_4 linsi_es_3 linsi_i_fftns2-op3; do
    [ -s $rep/variants/$v/out.fasta ] || return 0; done
  mkdir -p $T
  ln -sf $rep/true.fasta $T/true.fasta
  ln -sf $rep/variants/linsi/out.fasta $T/magus.fasta
  ln -sf $rep/variants/wsoft0.03_c_linsi_i_fftns2_es_4/out.fasta $T/recipe.fasta
  ln -sf $rep/variants/linsi_es_3/out.fasta $T/es3.fasta
  ln -sf $rep/variants/linsi_i_fftns2-op3/out.fasta $T/hard.fasta
  b=${n%_d[0-9]*}; lvl=${b%_R*}; r=${b#*_R}
  nice -n 10 /opt/mm/root/envs/pasta183/bin/python /home/user/scratch/cs581/protbench/code/trees.py $T \
    /opt/data/sim/$lvl/R$r/tree.nwk $rep/trees.jsonl true magus recipe es3 hard >/dev/null
}
export -f one; export W
while true; do
  ls -d $W/reps/* 2>/dev/null | xargs -r -P ${1:-2} -I{} bash -c 'one {}'
  bash $G/code/collect.sh
  [ -z "${LOOP:-}" ] && break
  [ -f $W/STOP_TREES ] && break
  sleep 120
done
echo TREES_PASS_DONE
