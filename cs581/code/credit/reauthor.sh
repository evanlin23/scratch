#!/bin/bash
# AI-assisted (Claude), helper script for the CS581 project repository.
# Re-authors the Claude-authored commits of the finished cs581 exploration branches to Evan Lin, on COPIES
# named credit/<branch> (the original branches are left untouched). Every re-authored commit keeps a
# "Co-authored-by: Claude <noreply@anthropic.com>" line, as the course's AI-disclosure rule requires.
# Commits already on the default branch are not touched.
#
# Run it yourself from a fresh clone:
#   git clone https://github.com/evanlin23/scratch && cd scratch && bash cs581/code/credit/reauthor.sh
# It ends by pushing the credit/* branches; they are then merged into the default branch with PRs.
set -euo pipefail
DEFAULT=origin/claude/charming-pasteur-yl6k2v
BRANCHES="alncrit astralpro2 camus consensus-1 consensus-2 consensus-3 consensus-4 decodiphy2 disco dldtm
e2e-1 e2e-2 e2e-3 e2e-4 e2e-5 e2e-6 e2e-7 fastmagus forest gdl gdlcons gtm lenhet ling ml models pairmerge
quartets rogue samplecx sota2026 supertree timing-a timing-b worker-1 worker-2 worker-3 worker-4 worker-5
worker-6 worker-7 worker-8 worker-9 worker-10 worker-11 worker-12"

git fetch origin
for b in $BRANCHES; do git branch -f "credit/cs581-$b" "origin/claude/cs581-$b"; done
refs=$(git for-each-ref --format='%(refname)' refs/heads/credit)
echo "Claude-authored commits to re-author: $(git log --format=%ae $refs --not $DEFAULT | grep -c anthropic || true)"

FILTER_BRANCH_SQUELCH_WARNING=1 git filter-branch -f --env-filter '
if [ "$GIT_AUTHOR_EMAIL" = "noreply@anthropic.com" ]; then
  export WAS_CLAUDE=1
  export GIT_AUTHOR_NAME="Evan Lin" GIT_AUTHOR_EMAIL="113861384+evanlin23@users.noreply.github.com"
fi
if [ "$GIT_COMMITTER_EMAIL" = "noreply@anthropic.com" ]; then
  export GIT_COMMITTER_NAME="Evan Lin" GIT_COMMITTER_EMAIL="113861384+evanlin23@users.noreply.github.com"
fi' --msg-filter '
m=$(cat)
printf "%s\n" "$m"
if [ "${WAS_CLAUDE:-}" = "1" ] && ! printf "%s" "$m" | grep -qi "co-authored-by: claude"; then
  printf "\nCo-authored-by: Claude <noreply@anthropic.com>\n"
fi' -- $refs --not $DEFAULT

# Checks: no Claude authors left, and every copy has exactly the same files as its original.
left=$(git log --format=%ae $refs --not $DEFAULT | grep -c anthropic || true)
echo "Claude-authored commits left: $left"
for b in $BRANCHES; do
  [ "$(git rev-parse "credit/cs581-$b^{tree}")" = "$(git rev-parse "origin/claude/cs581-$b^{tree}")" ] || { echo "TREE MISMATCH $b"; exit 1; }
done
echo "all $(echo $BRANCHES | wc -w) copies match their originals file-for-file"
git push origin 'refs/heads/credit/*:refs/heads/credit/*'
