#!/bin/bash
# Render proposal.md to proposal.pdf (pandoc -> HTML -> headless Chromium), US Letter.
set -euo pipefail
cd "$(dirname "$0")"
pandoc proposal.md -s --metadata pagetitle="CS581 proposal" -c style.css --embed-resources -o proposal.html
CHROME=$(ls -d /opt/pw-browsers/chromium-*/chrome-linux/chrome 2>/dev/null | head -1)
"$CHROME" --headless --no-sandbox --disable-gpu --no-pdf-header-footer \
  --print-to-pdf=proposal.pdf "file://$PWD/proposal.html" 2>/dev/null
pdfinfo proposal.pdf | grep Pages
