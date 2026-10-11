#!/bin/bash
# QR data with fewer genes (first K gene-family trees, 50 bp gene trees): harder regime. Usage: bash batch_qr_k.sh K [P]
D=/opt/data/disco/qr/trees; K=$1; R=/opt/runs/qr_k$K; C=$(cd "$(dirname "$0")" && pwd); P=${2:-4}
jobs=()
for c in 20_gdl_5e-10_1 20_gdl_1e-9_1 20_gdl_1e-10_1 20_gdl_1e-10_0 20_gdl_5e-10_0; do  # 20_gdl_1e-9_0 (mean 1800 leaves) only reps run manually
 for ils in _hILS ""; do
  for r in 01 02 03 04 05 06 07 08 09 10; do
    jobs+=("$D/$c$ils/$r $R/$c$ils/$r/g_50")
  done; done; done
printf '%s\n' "${jobs[@]}" | xargs -P $P -L 1 bash -c 'mkdir -p "$(dirname "$1")"; [ -s "$1/result.json" ] || timeout 7200 python3 '"$C"'/run_rep.py "$0" g_50.trees "$1" --maxgenes '"$K"' --tags --iter > /dev/null 2> "$1.err"'
