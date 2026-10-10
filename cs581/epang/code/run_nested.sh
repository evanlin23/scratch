#!/bin/bash
# Run EPA-ng for every (center, size, query type, mode) of the nested experiment; restartable.
# (ONLY_ALL=1: only the full-backbone runs, nested/all; run those alone, they need ~10 GB)
# Usage: run_nested.sh <datadir> "<modes>" "<qtypes>" [parallel=4]
#   modes: auto rs_on rs_off noheur baseball noheur_rsoff reest ...
D=$1; MODES=$2; QT=$3; PAR=${4:-4}
P=/opt/mm/root/envs/place/bin
EPA=${EPA_BIN:-$P/epa-ng}
flags() {
  case $1 in
    auto) echo "";; rs_on) echo "--rate-scalers on";; rs_off) echo "--rate-scalers off";;
    noheur) echo "--no-heur";; baseball) echo "--baseball-heur";; nopremask) echo "--no-pre-mask";;
    noheur_rsoff) echo "--no-heur --rate-scalers off";; reest) echo "";;
    fix_rson) echo "--rate-scalers on";; fix*) echo "";;
  esac
}
jobs=()
for m in $MODES; do for q in $QT; do
  if [ "${ONLY_ALL:-0}" = 1 ]; then jobs+=("$D/nested/all|$m|$q"); continue; fi
  for kd in "$D"/nested/c*/k*; do jobs+=("$kd|$m|$q"); done
done; done
run_one() {
  IFS='|' read kd m q <<< "$1"
  out=$kd/${m}_$q
  [ -s $out/epa_result.jplace ] && return 0
  mkdir -p $out
  if [ "$kd" = "$D/nested/all" ]; then tree=$D/rx.raxml.bestTree; ref=$D/backbone.fa; model=$D/rx.raxml.bestModel
  else tree=$kd/tree.nwk; ref=$kd/ref.fa; model=$D/rx.raxml.bestModel; fi
  if [ $m = reest ]; then
    if [ ! -s $kd/reest.raxml.bestModel ]; then
      $P/raxml-ng --evaluate --msa $ref --tree $tree --model GTR+G --prefix $kd/reest --threads 1 --redo --log ERROR > /dev/null 2>&1
    fi
    tree=$kd/reest.raxml.bestTree; model=$kd/reest.raxml.bestModel
  fi
  bin=$EPA; [ "${m#fix}" != "$m" ] && bin=${EPA_FIX}
  s=$(date +%s.%N)
  $bin -t $tree -s $ref -q $kd/$q.fa -m $model -w $out -T ${THREADS:-1} --redo $(flags $m) > $out/log.txt 2>&1
  rc=$?
  [ $rc != 0 ] && rm -f $out/epa_result.jplace
  e=$(date +%s.%N)
  echo -e "$kd\t$m\t$q\t$(echo "$e - $s" | bc)\t$rc" >> $D/nested/times.tsv
}
export -f run_one flags; export D P EPA EPA_FIX
printf '%s\n' "${jobs[@]}" | xargs -P $PAR -I{} bash -c 'run_one "{}"'
