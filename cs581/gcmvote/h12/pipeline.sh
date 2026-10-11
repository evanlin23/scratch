#!/bin/bash
# AI-assisted (Claude), exploration code for CS581 project
# After the MAGUS draw: the first merge (gg.py linsi) alone, to measure its peak memory; then the other
# baselines and every pre-declared vote variant (B = 10; hard+mask skipped: its scored alignment is `hard`)
# 2-3 at a time (per the orchestrator), as memory allows. Each job runs under memwatch (wall, peak RSS, 4 h).
# Concurrent jobs share the 4 cores, so their walls are not single-job walls (recorded as "lanes").
# Vote rows go to a second rep dir (symlinked inputs) because run.py's results.jsonl rows carry "B".
# Restartable: finished jobs are skipped.
W=/opt/work/h12; H=/home/user/scratch/cs581/gcmvote/h12; REP=$W/reps/16S.T_R0; V=$W/reps/16S.T_R0v
VCODE=$W/vcode/cs581/gcmvote/code
until grep -q MAGUS_DONE $W/magus.out; do grep -q MAGUS_FAILED $W/magus.out && exit 1; sleep 30; done
mkdir -p $V $W/mem
for f in inputs true.fasta unaligned.fasta magus.json sets; do [ -e $V/$f ] || ln -s $REP/$f $V/$f; done
[ -f $W/inputs.tar.gz ] || tar czf $W/inputs.tar.gz -C $REP inputs/subalignments inputs/backbones
find $REP/variants -name "*.lock" -delete 2>/dev/null
cd /home/user/scratch/cs581/code
job() {  # job KIND VARIANT LANES
  local s=$(echo "$1_$2" | sed 's/[^A-Za-z0-9._-]/_/g')
  if [ "$1" = gg ]; then
    [ -f $REP/results.jsonl ] && grep -qF "\"variant\": \"$2\"," $REP/results.jsonl && return
    python3 $H/memwatch.py $W/mem/$s.json 14400 python3 /home/user/scratch/cs581/gcmgen/code/gg.py run $REP "$2"
  else
    [ -f $V/results.jsonl ] && grep -qF "\"variant\": \"$2\", \"B\": 10" $V/results.jsonl && return
    python3 $H/memwatch.py $W/mem/$s.json 14400 python3 $VCODE/run.py $V "$2"
  fi
  python3 -c "import json,sys; p=sys.argv[1]; r=json.load(open(p)); r['lanes']=int(sys.argv[2]); json.dump(r,open(p,'w'))" $W/mem/$s.json $3
}
export -f job; export W H REP V VCODE
job gg linsi 1
echo FIRST_MERGE_DONE
peak=$(python3 -c "import json; print(json.load(open('$W/mem/gg_linsi.json'))['peak_rss_gb'])")
P=$(python3 -c "print(max(1, min(3, int(12 // max($peak, 0.5)))))")
echo "LANES $P (first merge peak ${peak} GB)"
printf "%s\n" "gg linsi#es3" "gg linsi#es4" "gg linsi#es5" "vote magus" "vote es4" "vote hard" "vote hard-bb" \
  "vote soft" "vote soft-bb" "vote soft2" "vote soft4" "vote frac0.2" "vote frac0.3" "vote frac0.4" "vote frac0.5" \
  | xargs -P $P -L 1 bash -c 'job "$0" "$1" '$P
echo PIPELINE_DONE
