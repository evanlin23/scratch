#!/bin/bash
# Run WITCH's query-alignment + merge stage with a precomputed weights.txt and an existing eHMM dir.
# usage: witch_stage.sh INST HMMDIR WEIGHTS OUTDIR [K]
set -e
INST=$1; H=$2; WT=$3; OUT=$4; K=${5:-10}
rm -rf "$OUT"; mkdir -p "$OUT"
cp "$WT" "$OUT/weights.txt"
s=$(date +%s.%N)
witch.py -p "$H" -b "$INST/backbone.fasta" -q "$INST/queries.fasta" -d "$OUT" -o aln.fasta -t 4 \
  --molecule dna -k "$K" > "$OUT/stage.log" 2>&1
e=$(date +%s.%N)
echo "$e - $s" | bc > "$OUT/stage_wall.txt"
