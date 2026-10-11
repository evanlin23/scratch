#!/bin/bash
# all published cases, 2 at a time; restartable (skips cases with result.json)
cd "$(dirname "$0")"
jobs=()
for c in 1000M1-HF Cox1-HET RNASim1000; do
  case $c in 1000M1-HF) R="R0 R1 R2 R3 R4";; Cox1-HET) R="R01 R02 R03 R04 R05 R06 R07 R08 R09 R10";; RNASim1000) R="1 2 3 4 5";; esac
  for r in $R; do for g in FT IQ; do
    [ -f /opt/gtms/pub/$c/$g/$r/result.json ] || jobs+=("$c $r $g")
  done; done
done
printf '%s\n' "${jobs[@]}" | xargs -P ${P:-2} -L1 sh -c 'python3 pub.py /opt/gtms/pub $0 $1 $2 >> /opt/gtms/pub.log 2>&1'
