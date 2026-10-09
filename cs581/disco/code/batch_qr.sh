#!/bin/bash
# QR (DISCO+QR data, 21 species) batch: restartable; skips finished jobs. Usage: bash batch_qr.sh
D=/opt/data/disco/qr/trees; R=/opt/runs/qr; C=$(cd "$(dirname "$0")" && pwd)
jobs=()
for c in 20_gdl_1e-10_0 20_gdl_1e-10_1 20_gdl_5e-10_0 20_gdl_5e-10_1 20_gdl_1e-9_0 20_gdl_1e-9_1 20_gdl_1e-11_1; do
 for ils in "" _hILS; do
  for r in 01 02 03 04 05 06 07 08 09 10; do
   for g in g_100 g_50 g_true; do
    jobs+=("$D/$c$ils/$r $g.trees $R/$c$ils/$r/$g")
   done; done; done; done
printf '%s\n' "${jobs[@]}" | xargs -P 4 -L 1 bash -c 'mkdir -p "$(dirname "$2")"; [ -s "$2/result.json" ] || timeout 7200 python3 '"$C"'/run_rep.py "$0" "$1" "$2" --tags --iter > /dev/null 2> "$2.err"'
