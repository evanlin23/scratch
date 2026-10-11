#!/bin/bash
# PICRUSt2 2.6.3 full pipeline (bacterial reference, 26,868 tips -> EPA-ng rate scalers on)
# with stock vs fixed EPA-ng (each build put first on PATH as `epa-ng`).
# Usage: run_picrust2.sh <seqs.fna> <table.biom> <outroot> [threads=4]
S=$1; B=$2; O=$3; T=${4:-4}
E=/opt/mm/root/envs/picrust2/bin
for v in stock fix; do
  [ -s $O/$v/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz ] && continue
  rm -rf $O/$v; mkdir -p $O
  PATH=/opt/tools/epa/path_$v:$E:$PATH /usr/bin/time -v -o $O/$v.time \
    $E/picrust2_pipeline.py -s $S -i $B -o $O/$v -p $T --verbose > $O/$v.log 2>&1
done
