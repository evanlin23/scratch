#!/bin/bash
# Pilot on all published replicates: (1) constrained parsimony SPR from GTM (radius 6),
# (2) parsimony constrained insertion + SPR polish (radius 6). 4 jobs in parallel.
cd "$(dirname "$0")"
OUT=${1:-/opt/gtmdata/pilot}
jobs=()
for c in 1000M1-HF Cox1-HET RNASim1000; do
  case $c in RNASim1000) R="1 2 3 4 5";; Cox1-HET) R="R01 R02 R03 R04 R05 R06 R07 R08 R09 R10";; *) R="R0 R1 R2 R3 R4";; esac
  for r in $R; do for g in FT IQ; do
    jobs+=("python3 run_blend.py $c $r $g $OUT 6 > /dev/null 2>&1; python3 run_insert.py $c $r $g pars $OUT 6 > /dev/null 2>&1")
  done; done
done
printf '%s\n' "${jobs[@]}" | xargs -P 4 -I{} bash -c "{}"
echo ALLDONE
