#!/bin/bash
# after after_loo2.sh: (1) re-run hsp for s1 with R on PATH and score s1; (2) hmp placement-only
# comparison with --chunk-size 500 (the patched run OOMs at PICRUSt2's 5000 on this 15 GB machine)
while pgrep -f after_loo2.sh >/dev/null; do sleep 15; done
E=/opt/mm/root/envs/picrust2; C=/home/user/scratch/cs581/picrust/code; S=/opt/work/loo/s1
export PATH=$E/bin:$PATH
cd $S
for v in stock fix; do
  s=$(date +%s); hsp.py -t placed_$v.tre -o ko_$v.tsv.gz --observed_trait_table ko_reduced.txt.gz -n -p 4 > hsp_$v.log 2>&1
  echo -e "$S\t$v\thsp_rerun2=$(( $(date +%s)-s ))" >> /opt/work/loo/times.tsv
done
python $C/score_loo.py $S > $S/summary.json 2> $S/score.err
R=$E/lib/python3.12/site-packages/picrust2/default_files/bacteria/bac_ref
I=/opt/work/hmp/pipe_stock/intermediate/place_seqs_bac
for v in stock fix; do
  mkdir -p /opt/work/hmp/epa500_$v
  s=$(date +%s)
  /opt/epa/$v/epa-ng.bin --tree $R/bac_ref.tre --ref-msa $I/ref_seqs_hmmalign.fasta --query $I/study_seqs_hmmalign.fasta \
    --chunk-size 500 -T 4 -m $R/bac_ref.model -w /opt/work/hmp/epa500_$v --filter-acc-lwr 0.99 --filter-max 100 --redo > /opt/work/hmp/epa500_$v/stdout 2>&1
  echo -e "hmp_epa500\t$v\t$(( $(date +%s)-s ))\t$?" >> /opt/work/pipe_times.tsv
done
