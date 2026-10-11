#!/bin/bash
# sequential: WITCH alignments for Q3 (MAFFT backbone)
cd "$(dirname "$0")"; export MLDATA=/opt/data/fscache
for x in "M1HF 0" "M1HF 1" "M1HF 2" "M1HF 3" "M1HF 4" "RNASimHF 0"; do
  /opt/mm/root/envs/fml/bin/python align_witch.py $x
done
