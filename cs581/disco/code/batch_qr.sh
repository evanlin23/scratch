#!/bin/bash
# QR (DISCO+QR data, 21 species) batch: restartable; skips finished jobs. Usage: bash batch_qr.sh [P]
# Estimated gene trees (50, 100 bp): reps 01-10; true gene trees (tagging/root accuracy): reps 01-05.
D=/opt/data/disco/qr/trees; R=/opt/runs/qr; C=$(cd "$(dirname "$0")" && pwd); P=${1:-3}
jobs=()
for c in 20_gdl_5e-10_1 20_gdl_1e-9_1 20_gdl_1e-10_1 20_gdl_1e-10_0 20_gdl_5e-10_0; do  # 20_gdl_1e-9_0 (mean 1800 leaves) only reps run manually
 for ils in _hILS ""; do
  for r in 01 02 03 04 05 06 07 08 09 10; do
   for g in g_50 g_100 g_true; do
    [ $g = g_true ] && [ $r \> 05 ] && continue
    jobs+=("$D/$c$ils/$r $g.trees $R/$c$ils/$r/$g")
   done; done; done; done
printf '%s\n' "${jobs[@]}" | xargs -P $P -L 1 bash -c 'mkdir -p "$(dirname "$2")"; [ -s "$2/result.json" ] || timeout 7200 python3 '"$C"'/run_rep.py "$0" "$1" "$2" --tags --iter > /dev/null 2> "$2.err"'
