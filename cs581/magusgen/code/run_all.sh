#!/bin/bash
# Restartable driver for the merge-only comparisons (rerun after an interruption; finished rows are skipped).
#   bash run_all.sh [VARIANTS...]       default variant list below
set -u
M=/home/user/scratch/cs581/magusgen
C=/home/user/scratch/cs581/code
W=/opt/work/mg
mkdir -p $W/reps
while read -r name src k; do
  [ -z "$name" ] && continue
  R=$W/reps/$name
  if [ ! -d $R/sets/s0 ]; then
    if [ "$k" = fresh ]; then  # no cached MAGUS run: fresh MAGUS draw with the paper's flags
      echo "$name $src 25" > $W/job_$name.txt
      (cd $C && python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus_fresh.jsonl $W/fresh --draws 0 --tools '' \
        --e2e '' --union '' --threads 4) || continue
      python3 $M/../protcons/code/pc.py rep $W/fresh/${name}_d0 $R || continue
    else
      python3 $M/code/mg.py rep $name $R $src || continue
    fi
  fi
  V=$(grep -v '^#' $M/code/variants.txt | tr '\n' ' ')  # re-read per job, so variants can be added mid-run
  (cd $C && python3 $M/code/mg.py run $R $V) >> $W/run.log 2>&1
  mkdir -p $M/results
  cp $R/results.jsonl $M/results/$name.results.jsonl
done < ${JOBS:-$M/code/jobs.txt}
echo ALLDONE
