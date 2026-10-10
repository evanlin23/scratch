#!/bin/bash
# FastTree -lg -gamma on the true alignment, MAGUS (linsi control merge) and the filtered alignments of each
# simulated dataset; nRF to the true tree (protbench/code/trees.py). Runs up to 4 datasets in parallel
# (FastTree is single-threaded). Restartable.
P=/home/user/scratch/cs581/protcons
W=/opt/work/protcons
one() {
  n=$1; R=$W/reps/$n; T=$W/trees/${n}_d0
  [ -f $R/variants/linsi_m_cons0.7/out.fasta ] || return
  mkdir -p $T
  ln -sf $R/true.fasta $T/true.fasta
  ln -sf $R/variants/linsi/out.fasta $T/magus.fasta
  ln -sf $R/variants/linsi_m_cons0.7/out.fasta $T/cons.fasta
  ln -sf $R/variants/linsi_i_fftns2-op3/out.fasta $T/fftint.fasta
  ln -sf $R/variants/linsi_p_clustalo_m_cons0.7/out.fasta $T/lccons.fasta
  r=${n#SIM*_}; lvl=${n%_R*}
  /opt/mm/root/envs/pasta183/bin/python $P/../protbench/code/trees.py $T /opt/data/sim/$lvl/$r/tree.nwk $R/trees.jsonl \
    true magus cons fftint lccons
}
export -f one; export P W
ls $W/reps | grep '^SIM' | xargs -P 4 -I{} bash -c 'one {}'
