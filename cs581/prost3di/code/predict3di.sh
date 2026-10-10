#!/bin/bash
# Predict 3Di strings from amino-acid FASTA with ProstT5 (Foldseek's CPU gguf path).
# Usage: predict3di.sh OUTDIR in1.tfa [in2.tfa ...]   (restartable: skips sets already done)
# Writes OUTDIR/<set>.db/ (Foldseek db: AA + 3Di, no coordinates), OUTDIR/<set>.3di.fa, OUTDIR/<set>.aa.fa and OUTDIR/<set>.time (wall/user/sys seconds).
set -u
F=/opt/mm/root/envs/fs/bin/foldseek
W=/opt/prostt5/weights
OUT=$1; shift
mkdir -p "$OUT"
for f in "$@"; do
  b=$(basename "$f"); b=${b%.*}
  [ -s "$OUT/$b.3di.fa" ] && continue
  T=$(mktemp -d)
  # one sequence per line, header = first word
  awk '/^>/{if(s)print s; print $1; s=""; next}{gsub(/[ \t\r-]/,""); s=s toupper($0)}END{if(s)print s}' "$f" > "$T/in.fa"
  /usr/bin/time -f "%e %U %S" -o "$OUT/$b.time" $F createdb "$T/in.fa" "$T/db" --prostt5-model $W --threads 4 > "$T/log" 2>&1
  $F lndb "$T/db_h" "$T/db_ss_h" > /dev/null 2>&1
  $F convert2fasta "$T/db_ss" "$T/3di.fa" > /dev/null 2>&1
  $F convert2fasta "$T/db" "$T/aa.fa" > /dev/null 2>&1
  if [ -s "$T/3di.fa" ]; then cp "$T/aa.fa" "$OUT/$b.aa.fa"; mkdir -p "$OUT/$b.db"; cp -P "$T"/db* "$OUT/$b.db/"; mv "$T/3di.fa" "$OUT/$b.3di.fa"; echo "$b done $(cat $OUT/$b.time)"; else echo "$b FAILED"; tail -3 "$T/log"; fi
  rm -rf "$T"
done
