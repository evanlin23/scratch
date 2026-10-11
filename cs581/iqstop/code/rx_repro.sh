#!/bin/bash
# Paired RAxML-NG 2.0.3 runs, sequential (same load for both): --fast (KH-mult) vs classic SPR search, 1 parsimony start.
RX=/opt/mm/root/envs/bio/bin/raxml-ng; W=/opt/work/runs
for f in /opt/work/emp/emp16S*.fa $(ls /opt/work/sim/*.phy | head -8); do
  b=$(basename $f); b=${b%.*}; d=$W/$b
  for tag in rxfastB rxclassic1; do
    [ -f $d/$tag.done ] && continue
    if [ $tag = rxfastB ]; then opts="--fast"; else opts="--search --tree pars{1} --opt-topology classic"; fi
    s=$(date +%s.%N); $RX $opts --msa $f --model GTR+FO+IO+G4 --threads 1 --seed 1 --prefix $d/$tag --redo --log PROGRESS > $d/$tag.stdout 2>&1; rc=$?
    e=$(date +%s.%N); echo "$b $tag $rc $(echo "$e - $s" | bc)" > $d/$tag.done
  done
done
