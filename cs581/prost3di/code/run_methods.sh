#!/bin/bash
# Run baseline and predicted-3Di aligners on one BAliBASE set (restartable: skips existing outputs).
# Usage: run_methods.sh SETID OUTROOT [methods...]
#   SETID like BB11001; needs /opt/p3d/bb3/SETID.{aa,3di}.fa from predict3di.sh
# Methods (all single-threaded):
#   linsi      MAFFT L-INS-i on amino acids
#   clustalo   Clustal Omega (--threads 1)
#   mafft      MAFFT default (FFT-NS-2)
#   fm         FoldMason structuremsa on the ProstT5 3Di+AA database (no coordinates)
#   fm_aa      same FoldMason engine with the 3Di score switched off (--bitfactor-3di 0.01; 0 is rejected): AA-only control
#   linsi3di   MAFFT L-INS-i on predicted 3Di strings with Foldseek's 3Di matrix, mapped back to AA
#   mc3di      M-Coffee (T-Coffee consistency library) combining linsi + linsi3di + fm
#   fm_true, linsi3di_true   oracle controls: same pipelines with 3Di from the experimental PDB structure
#              (true3di.py) instead of ProstT5
#   mcaa       control: M-Coffee combining three AA-only alignments, linsi + clustalo + mafft
set -u
ID=$1; OUT=$2; shift 2
METHODS=${@:-linsi clustalo mafft fm fm_aa linsi3di mc3di mcaa fm_true linsi3di_true}
FS=/opt/mm/root/envs/fs/bin
P=/opt/p3d/bb3
CODE=$(cd "$(dirname "$0")" && pwd)
D=$OUT/$ID; mkdir -p "$D"
AA=$P/$ID.aa.fa; TDI=$P/$ID.3di.fa
TRUE=/opt/p3d/bb3true/$ID.3di.fa
case "$METHODS" in *_true*) [ -s $TRUE ] || python3 "$CODE/true3di.py" /opt/p3d/pdb/sdb_aa.fa /opt/p3d/pdb/sdb_3di.fa /opt/p3d/bb3true $AA > /dev/null ;; esac
T() { /usr/bin/time -f "%e %U %S" -o "$D/$1.time" "${@:2}"; }
for m in $METHODS; do
  [ -s "$D/$m.fa" ] && continue
  case $m in
    linsi)    T $m mafft --quiet --localpair --maxiterate 1000 --thread 1 "$AA" > "$D/$m.fa" ;;
    mafft)    T $m mafft --quiet --thread 1 "$AA" > "$D/$m.fa" ;;
    clustalo) T $m clustalo -i "$AA" --threads 1 --force -o "$D/$m.fa" ;;
    fm|fm_aa|fm_true)
      W=$(mktemp -d); cp $P/$ID.db/db $P/$ID.db/db.* $P/$ID.db/db_h* $P/$ID.db/db_ss $P/$ID.db/db_ss.* "$W/"
      [ $m = fm_true ] && python3 "$CODE/write_ss.py" "$W/db" $TRUE
      extra=""; [ $m = fm_aa ] && extra="--bitfactor-3di 0.01"
      T $m $FS/foldmason structuremsa "$W/db" "$W/out" --threads 1 $extra > "$W/log" 2>&1
      cp "$W/out_aa.fa" "$D/$m.fa"; rm -rf "$W" ;;
    linsi3di|linsi3di_true)
      IN=$TDI; [ $m = linsi3di_true ] && IN=$TRUE
      T $m mafft --quiet --localpair --maxiterate 1000 --thread 1 --aamatrix /opt/p3d/mat/mat3di.mafft "$IN" > "$D/$m.3di.fa"
      python3 "$CODE/map3di.py" "$D/$m.3di.fa" "$AA" > "$D/$m.fa" ;;
    mc3di|mcaa)
      [ $m = mc3di ] && ins="linsi linsi3di fm" || ins="linsi clustalo mafft"
      W=$(mktemp -d); args=""; for x in $ins; do args="$args A$D/$x.fa"; done
      (cd "$W" && T $m /opt/mm/root/envs/tc/bin/t_coffee -in $args -output fasta_aln -outfile "$W/out.fa" -n_core 1 -quiet > /dev/null 2>&1)
      [ -s "$W/out.fa" ] && cp "$W/out.fa" "$D/$m.fa"; rm -rf "$W" ;;
  esac
done
