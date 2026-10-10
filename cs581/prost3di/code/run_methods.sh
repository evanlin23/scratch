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
set -u
ID=$1; OUT=$2; shift 2
METHODS=${@:-linsi clustalo mafft fm fm_aa linsi3di}
FS=/opt/mm/root/envs/fs/bin
P=/opt/p3d/bb3
CODE=$(cd "$(dirname "$0")" && pwd)
D=$OUT/$ID; mkdir -p "$D"
AA=$P/$ID.aa.fa; TDI=$P/$ID.3di.fa
T() { /usr/bin/time -f "%e %U %S" -o "$D/$1.time" "${@:2}"; }
for m in $METHODS; do
  [ -s "$D/$m.fa" ] && continue
  case $m in
    linsi)    T $m mafft --quiet --localpair --maxiterate 1000 --thread 1 "$AA" > "$D/$m.fa" ;;
    mafft)    T $m mafft --quiet --thread 1 "$AA" > "$D/$m.fa" ;;
    clustalo) T $m clustalo -i "$AA" --threads 1 --force -o "$D/$m.fa" ;;
    fm|fm_aa)
      W=$(mktemp -d); cp $P/$ID.db/db $P/$ID.db/db.* $P/$ID.db/db_h* $P/$ID.db/db_ss $P/$ID.db/db_ss.* "$W/"
      extra=""; [ $m = fm_aa ] && extra="--bitfactor-3di 0.01"
      T $m $FS/foldmason structuremsa "$W/db" "$W/out" --threads 1 $extra > "$W/log" 2>&1
      cp "$W/out_aa.fa" "$D/$m.fa"; rm -rf "$W" ;;
    linsi3di)
      T $m mafft --quiet --localpair --maxiterate 1000 --thread 1 --aamatrix /opt/p3d/mat/mat3di.mafft "$TDI" > "$D/$m.3di.fa"
      python3 "$CODE/map3di.py" "$D/$m.3di.fa" "$AA" > "$D/$m.fa" ;;
  esac
done
