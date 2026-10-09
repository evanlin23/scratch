#!/bin/bash
# DISCO data (101 species) batch, restartable. Usage: bash batch_disco.sh [P]
# k = number of gene-family trees used (first k of 1000); the DISCO paper's default is 50.
D=/opt/data/disco/d/trees; R=/opt/runs/disco; C=$(cd "$(dirname "$0")" && pwd); P=${1:-1}
jobs=()
for spec in default:50 ils_2e8:50 gdl_1e-9_1:50 ils_1e4:50 gdl_1e-10_1:50; do
 c=${spec%%:*}; k=${spec##*:}
 for r in 01 02 03 04 05 06 07 08 09 10; do
  jobs+=("$D/$c/$r $R/$c-k$k/$r $k --iter")
 done; done
for r in 01 02 03 04 05; do jobs+=("$D/default/$r $R/default-k1000/$r 1000"); done
printf '%s\n' "${jobs[@]}" | xargs -P $P -L 1 bash -c 'mkdir -p "$(dirname "$1")"; [ -s "$1/result.json" ] || timeout 7200 python3 '"$C"'/run_rep.py "$0" g_100.trees "$1" --maxgenes "$2" --nsample 100 --tags $3 > /dev/null 2> "$1.err"'
