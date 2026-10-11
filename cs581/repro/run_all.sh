#!/bin/bash
# AI-assisted (Claude), code for CS581 project
#
# Recompute every core number from the rep bank (after fetch_bank.sh). Restartable: finished pieces are skipped.
#   bash cs581/repro/run_all.sh [aln] [bb] [trees] [tables]      (no argument = all four, in this order)
# Environment: JOBS = parallel jobs (default 4), DATA = bank root (default cs581/repro/data),
#              OUT = output root (default cs581/repro/out), VARIANTS = merge variants (default "magus es4 vote-hard-bb").
#
#   aln     every bank replicate at B = 10: filter_gcm.py with each variant, scored with FastSP
#           -> OUT/rows/SRC/REP/VARIANT@B10.json, alignment in OUT/aln/SRC/REP/VARIANT@B10.fasta
#   bb      backbone-count sensitivity on the replicates that have 5/10/20 backbones (h6, h6b, h9, h9b):
#           magus, frac0.2-0.5, vote-hard-bb at B = 5, 10, 20 -> OUT/rows/SRC/REP/VARIANT@B{5,10,20}.json
#   trees   FastTree -lg -gamma on the true alignment and the magus / es4 / vote-hard-bb alignments of every
#           replicate with a true tree (SIMHIGH_*, 1000M1_*)            -> OUT/trees/SRC/REP/METHOD.nwk
#   tables  python3 tables.py: tables + comparison with the logged rows -> OUT/TABLES.md, OUT/MISMATCHES.md
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DATA=${DATA:-$HERE/data}
OUT=${OUT:-$HERE/out}
JOBS=${JOBS:-4}
VARIANTS=${VARIANTS:-"magus es4 vote-hard-bb"}
BB_VARIANTS="magus frac0.2 frac0.3 frac0.4 frac0.5 vote-hard-bb"
stages=${*:-aln bb trees tables}
mkdir -p "$OUT"

# one merge job: SRC REP_DIR VARIANT B BACKBONE_DIR  ->  row JSON (skipped when it exists)
merge_job() {
  local src=$1 rep=$2 v=$3 B=$4 bb=$5
  local name; name=$(basename "$rep")
  local row="$OUT/rows/$src/$name/$v@B$B.json" aln="$OUT/aln/$src/$name/$v@B$B.fasta"
  [ -s "$row" ] && return 0
  mkdir -p "$(dirname "$row")" "$(dirname "$aln")"
  if python3 "$HERE/filter_gcm.py" "$rep" "$v" --backbones "$bb" -np 1 --out "$aln" --ref "$rep/true.fasta" \
       > "$row.tmp" 2> "$row.err"; then
    python3 -c "import json,sys; r=json.load(open(sys.argv[1])); r.update(src=sys.argv[2], rep=sys.argv[3], B=int(sys.argv[4])); \
json.dump(r, open(sys.argv[1], 'w'))" "$row.tmp" "$src" "$name" "$B"
    mv "$row.tmp" "$row"; rm -f "$row.err"
  else
    echo "FAILED $src/$name $v B=$B (see $row.err)"; rm -f "$row.tmp"
  fi
}

# one tree job: SRC REP_DIR METHOD ALIGNMENT  ->  FastTree newick (skipped when it exists)
tree_job() {
  local src=$1 rep=$2 m=$3 aln=$4
  local t="$OUT/trees/$src/$(basename "$rep")/$m.nwk"
  [ -s "$t" ] && return 0
  [ -s "$aln" ] || { echo "missing alignment $aln (run the aln stage first)"; return 0; }
  mkdir -p "$(dirname "$t")"
  if FastTree -lg -gamma -quiet "$aln" > "$t.tmp" 2> "$t.log"; then mv "$t.tmp" "$t"; else echo "FAILED tree $t"; fi
}
export -f merge_job tree_job
export HERE OUT

reps() { find "$DATA/bank" -mindepth 3 -maxdepth 3 -name inputs -path "*/bank/*/*/inputs" | sed 's|/inputs$||' | sort; }

for stage in $stages; do
  echo "== stage $stage ($(date -u +%H:%M:%S) UTC)"
  jobs_file="$OUT/jobs.$stage.txt"; : > "$jobs_file"
  case $stage in
    aln)
      for rep in $(reps); do
        src=$(basename "$(dirname "$rep")"); name=$(basename "$rep")
        case $name in *_B5|*_B20) continue;; esac          # h9/h9b B=5/20 copies: bb stage only
        bb=$rep/inputs/backbones; [ -d "$rep/bb10" ] && bb=$rep/bb10  # h6/h6b: 20-backbone draw, use its first 10
        for v in $VARIANTS; do echo "$src $rep $v 10 $bb" >> "$jobs_file"; done
      done ;;
    bb)
      for rep in $(reps); do
        src=$(basename "$(dirname "$rep")"); name=$(basename "$rep")
        if [ -d "$rep/bb5" ]; then                               # h6, h6b: bb5 / bb10 / bb20 folders
          for B in 5 10 20; do for v in $BB_VARIANTS; do echo "$src $rep $v $B $rep/bb$B" >> "$jobs_file"; done; done
        elif [[ $name =~ _B(5|10|20)$ ]]; then                   # h9, h9b: one replicate copy per B
          B=${BASH_REMATCH[1]}
          for v in $BB_VARIANTS; do echo "$src $rep $v $B $rep/inputs/backbones" >> "$jobs_file"; done
        fi
      done ;;
    trees)
      for rep in $(reps); do
        src=$(basename "$(dirname "$rep")"); name=$(basename "$rep")
        [ -s "$rep/true_tree.nwk" ] || continue
        case $name in SIMHIGH_*|1000M1_*) ;; *) continue;; esac
        [[ $src == gcmvote-* ]] || continue                   # gcmtrees reps: no vote rows were logged for them
        echo "$src $rep true $rep/true.fasta" >> "$jobs_file"
        for v in magus es4 vote-hard-bb; do echo "$src $rep $v $OUT/aln/$src/$name/$v@B10.fasta" >> "$jobs_file"; done
      done ;;
    tables)
      python3 "$HERE/tables.py" --data "$DATA" --out "$OUT"; continue ;;
    *) echo "unknown stage $stage"; exit 1 ;;
  esac
  fn=merge_job; [ "$stage" = trees ] && fn=tree_job
  echo "$(wc -l < "$jobs_file") jobs, $JOBS at a time"
  # bigger replicates first (16S.T / 16S.3 take ~10 min per merge), so the tail of the run stays parallel
  awk '{n=($2 ~ /16S\.(T|3)/) ? 0 : 1; print n "\t" $0}' "$jobs_file" | sort -s -k1,1n | cut -f2- \
    | xargs -P "$JOBS" -L 1 bash -c "$fn \"\$@\"" _
done
echo "== done ($(date -u +%H:%M:%S) UTC)"
