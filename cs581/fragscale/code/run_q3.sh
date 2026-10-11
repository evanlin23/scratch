#!/bin/bash
# Q3 tree arms on the WITCH alignments; waits for each replicate's alignment. Restartable.
#   bash run_q3.sh "M1HF 0" "M1HF 1" ...
cd "$(dirname "$0")"; export MLDATA=/opt/data/fscache
for x in "$@"; do
  set -- $x
  until [ -s $MLDATA/$1/R$2/witch.fasta ]; do sleep 20; done
  /opt/mm/root/envs/fml/bin/python pipe.py $1 $2 base_raxmlng place_ft_0.5_fix_rxfast base_fasttree constr_ft_0.5 --aln witch
done
