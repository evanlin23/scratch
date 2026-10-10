#!/bin/bash
# Restartable job runner: one line per (dataset, method); 4 concurrent single-threaded jobs.
# Usage: bash run_jobs.sh <phase>   phase = main (IQ-TREE default seed1 + RAxML-NG --fast) | extra (seed2, --fast, nstop sanity)
set -u
IQ=/opt/work/bin/iqtree3        # IQ-TREE 3.1.4 + per-iteration trace (patch in code/iqtree_trace.patch)
RX=/opt/mm/root/envs/bio/bin/raxml-ng
W=/opt/work/runs; mkdir -p $W
PHASE=$1
jobs=(); jobs2=()
for f in /opt/work/sim/*.phy /opt/work/emp/*.fa; do
  b=$(basename $f); b=${b%.*}
  case $b in empAA*) M=LG+G4; MR=LG+G4;; *) M=GTR+F+I+G4; MR=GTR+FO+IO+G4;; esac
  if [ $PHASE = main ]; then
    jobs+=("$b iqdef1 $IQ -s $f -m $M -T 1 -seed 1 -pre $W/$b/iqdef1 -redo --quiet")
    jobs+=("$b rxfast $RX --fast --msa $f --model $MR --threads 1 --seed 1 --prefix $W/$b/rxfast --redo --log PROGRESS")
  elif [ $PHASE = repro ]; then
    jobs+=("$b rxclassic1 $RX --search --tree pars{1} --opt-topology classic --msa $f --model $MR --threads 1 --seed 1 --prefix $W/$b/rxclassic1 --redo --log PROGRESS")
    case $b in rg42163_n94|rg51878_n319|rg6754_n258|rg21213_n117|emp16S_16S3_1_n150|empAA_RV100_BBA0067_n150)
      jobs+=("$b iqnstop20 $IQ -s $f -m $M -T 1 -seed 1 -nstop 20 -pre $W/$b/iqnstop20 -redo --quiet");; esac
  else
    jobs+=("$b iqfast $IQ -s $f -m $M -T 1 -seed 1 -pre $W/$b/iqfast -redo --quiet --fast")
    wl=$(cut -d' ' -f4 $W/$b/iqdef1.done 2>/dev/null | cut -d. -f1); [ -n "$wl" ] && [ "$wl" -lt 700 ] && \
    jobs2+=("$b iqdef2 $IQ -s $f -m $M -T 1 -seed 2 -pre $W/$b/iqdef2 -redo --quiet")
  fi
done
jobs+=("${jobs2[@]}")
run1() {
  b=$1; tag=$2; shift 2; d=/opt/work/runs/$b; mkdir -p $d
  [ -f $d/$tag.done ] && exit 0
  s=$(date +%s.%N); "$@" > $d/$tag.stdout 2>&1; rc=$?; e=$(date +%s.%N)
  echo "$b $tag $rc $(echo "$e - $s" | bc)" > $d/$tag.done
}
export -f run1
printf '%s\n' "${jobs[@]}" | xargs -P 4 -I{} bash -c 'run1 {}'
