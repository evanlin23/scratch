#!/bin/bash
# Restartable driver (rerun after an interruption; every step skips finished work).
#   bash run_all.sh            protein jobs (jobs.txt) then nucleotide controls (nuc.txt)
set -u
P=/home/user/scratch/cs581/protcons
C=/home/user/scratch/cs581/code
W=/opt/work/protcons
V='linsi linsi|cons0.7 linsi&fftns2-op3 linsi+clustalo|cons0.7'
mkdir -p $W/reps
cd $C
while read -r name src k unal; do
  [ -z "$name" ] && continue
  echo "$name $src $k $unal" > $W/job_$name.txt
  python3 -m gcmx.bbtool_bench $W/job_$name.txt $W/magus.jsonl $W --draws 0 --tools '' --e2e '' --union '' --threads 4 || continue
  python3 $P/code/pc.py rep $W/${name}_d0 $W/reps/$name || continue
  python3 $P/code/pc.py run $W/reps/$name $V
  python3 $P/code/pc.py gate $W/reps/$name
  bash $P/code/collect.sh
done < ${JOBS:-$P/jobs.txt}
while read -r name src; do
  [ -z "$name" ] && continue
  R=$W/reps/$name
  if [ ! -f $R/magus.json ]; then
    rm -rf $R; mkdir -p $R
    ref=$(git -C $P for-each-ref --format='%(refname:short)' 'refs/remotes/origin/claude/cs581-worker-*' | while read r; do
      git -C $P cat-file -e $r:cs581/experiments/runs/$name/inputs.tar.xz 2>/dev/null && echo $r && break; done)
    git -C $P show $ref:cs581/experiments/runs/$name/inputs.tar.xz | tar xJ -C $R
    python3 $P/../bbevidence/code/bbe.py prep $R $src
    git -C $P show $ref:cs581/experiments/runs/$name/prep.json > $R/magus.json
  fi
  python3 $P/code/pc.py run $R $V
  python3 $P/code/pc.py gate $R
  bash $P/code/collect.sh
done < ${NUC:-$P/nuc.txt}
echo ALLDONE
