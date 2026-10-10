#!/bin/bash
# PICRUSt2 2.6.3 steps (place_seqs with EPA-ng -> hsp 16S+NSTI, hsp KO/EC -> metagenome_pipeline)
# on a reduced bacterial reference (picrust2_subref.py), stock vs fixed EPA-ng on PATH.
# Trait tables are the default bacterial ones subset to the reference tips (hsp's mp fails otherwise).
# Usage: run_picrust2_subref.sh <workdir with asv.fna, table.biom, sub_ref/> [threads=4]
W=$1; T=${2:-4}; E=/opt/mm/root/envs/picrust2/bin
cd $W
for v in stock fix; do
  o=$W/$v; [ -s $o/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz ] && continue
  rm -rf $o; mkdir -p $o
  export PATH=/opt/tools/epa/path_$v:$E:/usr/bin:/bin
  s=$(date +%s)
  /usr/bin/time -v -o $o/place.time place_seqs.py -s asv.fna -o $o/placed.tre -p $T -r $W/sub_ref \
      --intermediate $o/intermediate/place_seqs --verbose > $o/place.log 2>&1
  hsp.py --observed_trait_table $W/sub_ref/16S.txt.gz -t $o/placed.tre -o $o/marker_predicted_and_nsti.tsv.gz -n -p $T > $o/hsp16.log 2>&1
  for tr in KO EC; do tl=$(echo $tr | tr A-Z a-z); hsp.py --observed_trait_table $W/sub_ref/$tl.txt.gz -t $o/placed.tre -o $o/${tr}_predicted.tsv.gz -p $T > $o/hsp$tr.log 2>&1
    metagenome_pipeline.py -i table.biom -m $o/marker_predicted_and_nsti.tsv.gz -f $o/${tr}_predicted.tsv.gz \
      -o $o/${tr}_metagenome_out > $o/mg$tr.log 2>&1; done
  echo -e "$v\t$(( $(date +%s) - s ))" >> $W/times.tsv
done
