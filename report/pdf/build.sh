#!/usr/bin/env bash
# Build report/report.pdf and report/supplement_A3.pdf from their Markdown sources.
#   bash report/pdf/build.sh            (from the repository root, or from anywhere)
# Needs pandoc (3.1 or later) and Chromium or Google Chrome. Markdown -> HTML with pandoc
# (filter.lua, style.css, images embedded); HTML -> PDF with headless Chromium (print.js through
# Playwright if Node and Playwright are available, otherwise Chromium's own --print-to-pdf).
# Set CHROME=/path/to/chrome to choose the browser.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
report="$(dirname "$here")"
command -v pandoc >/dev/null || { echo "pandoc not found" >&2; exit 1; }

find_chrome() {
  if [ -n "${CHROME:-}" ]; then echo "$CHROME"; return; fi
  for c in chromium chromium-browser google-chrome google-chrome-stable \
           /opt/pw-browsers/chromium-*/chrome-linux*/chrome \
           "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"; do
    if command -v "$c" >/dev/null 2>&1; then command -v "$c"; return; fi
    if [ -x "$c" ]; then echo "$c"; return; fi
  done
}

build() {
  local name="$1" title="$2"
  local html="$here/$name.html" pdf="$report/$name.pdf"
  pandoc "$report/$name.md" \
    --from markdown-implicit_figures+pipe_tables+lists_without_preceding_blankline \
    --to html5 --standalone --embed-resources --mathml \
    --resource-path "$report" \
    --lua-filter "$here/filter.lua" \
    --css "$here/style.css" \
    --metadata pagetitle="$title" \
    --metadata lang=en-US \
    --output "$html"
  # pandoc 3.1 marks every parenthesis in MathML as stretchy; inline ones then grow around hats
  # and subscripts. Nothing in the report needs stretching, so switch it off.
  sed -i 's/stretchy="true"/stretchy="false"/g' "$html"
  if command -v node >/dev/null && node -e "require('playwright')" >/dev/null 2>&1; then
    node "$here/print.js" "$html" "$pdf" "$title"
  else
    local chrome; chrome="$(find_chrome)"
    [ -n "$chrome" ] || { echo "Chromium/Chrome not found; set CHROME=" >&2; exit 1; }
    "$chrome" --headless --no-sandbox --disable-gpu --no-pdf-header-footer \
      --print-to-pdf="$pdf" "file://$html" 2>/dev/null
  fi
  echo "Wrote $pdf"
}

build report "A wider window, the same wall: Two decades of item-level DIF tests for binary items"
build supplement_A3 "Supplement A3: retrospective ADEMP-PreReg form"
