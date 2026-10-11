#!/bin/bash
# pplacer on every nested subtree (FastTree GTR+Gamma20 numeric parameters via taxtastic, as in
# BSCAMPP(p)/SCAMPP(p)). Usage: run_pplacer.sh <datadir> "<qtypes>" [parallel=2]
D=$1; QT=$2; PAR=${3:-2}
P=/opt/mm/root/envs/place/bin
run_pp() {
  kd=$1; q=$2; out=$kd/pplacer_$q
  [ -s $out/pplacer.jplace ] && return 0
  mkdir -p $out
  [ -s $kd/ref_T.fa ] || tr 'U' 'T' < $kd/ref.fa > $kd/ref_T.fa
  [ -s $kd/${q}_T.fa ] || tr 'U' 'T' < $kd/$q.fa > $kd/${q}_T.fa
  if [ ! -d $kd/pp.refpkg ]; then
    $P/taxit create -l pp -P $kd/pp.refpkg --aln-fasta $kd/ref_T.fa --tree-file $kd/tree_ft.nwk --tree-stats $D/ft.log > /dev/null 2>&1
  fi
  s=$(date +%s.%N)
  $P/pplacer -m GTR -c $kd/pp.refpkg -o $out/pplacer.jplace -j 1 $kd/${q}_T.fa > $out/log.txt 2>&1 || rm -f $out/pplacer.jplace
  echo -e "$kd\tpplacer\t$q\t$(echo "$(date +%s.%N) - $s" | bc)\t0" >> $D/nested/times.tsv
}
export -f run_pp; export D P
for q in $QT; do for kd in $D/nested/c*/k*; do echo "$kd $q"; done; done | xargs -P $PAR -n 2 bash -c 'run_pp "$0" "$1"'
