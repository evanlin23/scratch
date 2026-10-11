#!/bin/bash
# AI-assisted (Claude), code for CS581 project
#
# Collect the rep bank (MAGUS merge inputs) and the logged result rows from the experiment branches.
#   bash cs581/repro/fetch_bank.sh [DEST]          (DEST defaults to cs581/repro/data)
#
# For every remote branch claude/cs581-gcmvote-* and claude/cs581-gcmtrees*  (SRC = gcmvote-h5, gcmtrees-h1, ...):
#   DEST/bank/tar/SRC/REP.tar.gz   the tarball, checked against the sha256 in that branch's MANIFEST.tsv
#   DEST/bank/SRC/REP/             the extracted replicate (inputs/subalignments, inputs/backbones, true.fasta, ...)
#   DEST/logged/SRC/...            the result rows that branch logged (*.jsonl under results_*), for tables.py
# plus DEST/logged/gcmvote/ (main gcmvote branch: results/raw/*.results.jsonl, results/tables.md).
#
# Path fixes: subsets.json (MAGUS's subset order) is rewritten to paths relative to the replicate (several were
# saved as absolute paths of the machine that ran them) and created from the tarball's file order where it is
# missing (step 4); any other absolute path in the small JSON files is reported. Re-running is safe.
set -euo pipefail
REPO=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
DEST=$(realpath -m "${1:-$REPO/cs581/repro/data}")
mkdir -p "$DEST/bank/tar" "$DEST/logged"
cd "$REPO"

# 1. fetch the branches (read only; nothing is checked out)
git fetch -q origin 'refs/heads/claude/cs581-gcmvote*:refs/remotes/origin/claude/cs581-gcmvote*' \
                    'refs/heads/claude/cs581-gcmtrees*:refs/remotes/origin/claude/cs581-gcmtrees*'
branches=$(git for-each-ref --format='%(refname:short)' 'refs/remotes/origin/claude/cs581-gcmvote*' \
                                                       'refs/remotes/origin/claude/cs581-gcmtrees*')

ntar=0; nbad=0
for br in $branches; do
  src=${br#origin/claude/cs581-}
  # 2. logged rows of this branch: only its own results_<suffix> folder (helper branches also carry copies of
  #    other branches' folders, which are skipped here so each row is taken once, from its source)
  suffix=${src#gcmvote}; suffix=${suffix#gcmtrees}; suffix=${suffix#-}
  top=gcmvote; [[ $src == gcmtrees* ]] && top=gcmtrees
  if [ -z "$suffix" ]; then
    pattern="^cs581/$top/(results/raw/.*\.results\.jsonl|results/tables\.md|results_h[0-9]+/.*\.jsonl)$"
  else
    pattern="^cs581/$top/results_${suffix}/.*\.jsonl$"
  fi
  for f in $(git ls-tree -r --name-only "$br" | grep -E "$pattern" || true); do
    out="$DEST/logged/$src/${f#cs581/$top/}"
    mkdir -p "$(dirname "$out")"
    git show "$br:$f" > "$out"
  done

  # 3. bank tarballs of this branch
  manifest="cs581/$top/bank/MANIFEST.tsv"
  git cat-file -e "$br:$manifest" 2>/dev/null || continue
  git show "$br:$manifest" > "$DEST/bank/tar/$src.MANIFEST.tsv"
  for f in $(git ls-tree -r --name-only "$br" "cs581/$top/bank/" | grep '\.tar\.gz$'); do
    rep=$(basename "$f" .tar.gz)
    tgz="$DEST/bank/tar/$src/$rep.tar.gz"
    mkdir -p "$DEST/bank/tar/$src" "$DEST/bank/$src"
    [ -s "$tgz" ] || git show "$br:$f" > "$tgz"
    # sha256 check against the MANIFEST row of this replicate (the 64-hex-digit field)
    want=$(awk -F'\t' -v r="$rep" '$1==r{for(i=2;i<=NF;i++) if($i ~ /^[0-9a-f]{64}$/) print $i}' \
           "$DEST/bank/tar/$src.MANIFEST.tsv" | head -1)
    have=$(sha256sum "$tgz" | cut -d' ' -f1)
    if [ -n "$want" ] && [ "$want" != "$have" ]; then
      echo "SHA256 MISMATCH $src/$rep (manifest $want, file $have)"; nbad=$((nbad + 1))
    fi
    if [ ! -d "$DEST/bank/$src/$rep/inputs" ]; then
      tar xzf "$tgz" -C "$DEST/bank/$src"
    fi
    ntar=$((ntar + 1))
  done
done

# 4. subset order and path fixes. MAGUS numbers the subsets in the order os.listdir returned them on the machine
#    that ran it, and the result depends (slightly) on that order. subsets.json holds it where the run saved it;
#    otherwise it is taken from the order of the files in the tarball (tar also lists a directory with readdir;
#    where both exist they agree on all but one replicate, SIMHIGH_R12). subsets.source records which was used.
python3 - "$DEST/bank" <<'EOF'
import glob, json, os, sys, tarfile
bank = sys.argv[1]
for rep in sorted(glob.glob(os.path.join(bank, "*", "*", "inputs"))):
    rep = os.path.dirname(rep)
    src, name = rep.split(os.sep)[-2:]
    order_file = os.path.join(rep, "subsets.json")
    if os.path.exists(os.path.join(rep, "subsets.source")):
        continue  # already fixed by an earlier run
    if os.path.exists(order_file):
        names, how = [os.path.basename(p) for p in json.load(open(order_file))], "subsets.json saved by the run"
    else:
        with tarfile.open(os.path.join(bank, "tar", src, name + ".tar.gz")) as t:
            names = [os.path.basename(m) for m in t.getnames() if "/inputs/subalignments/" in m and m.endswith(".txt")]
        how = "order of the files in the tarball"
    have = set(os.listdir(os.path.join(rep, "inputs", "subalignments")))
    assert set(names) == have, "subset files differ from the recorded order in " + rep
    json.dump(["inputs/subalignments/" + n for n in names], open(order_file, "w"))
    open(os.path.join(rep, "subsets.source"), "w").write(how + "\n")
for js in glob.glob(os.path.join(bank, "*", "*", "*.json")):
    if os.path.basename(js) != "subsets.json" and os.path.getsize(js) < 1e6:
        txt = open(js).read()
        if '"/opt/' in txt or '"/home/' in txt or '"/tmp/' in txt:
            print("note: absolute path left in", os.path.relpath(js, bank), "(not used by filter_gcm.py)")
EOF

echo "bank: $ntar replicates in $DEST/bank ($nbad sha256 mismatches); logged rows in $DEST/logged"
