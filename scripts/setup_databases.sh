#!/usr/bin/env bash
# Simple helper to download reference DB files (SILVA, PR2, optional BOLD)
# Usage: ./scripts/setup_databases.sh [DEST_DIR]

set -euo pipefail

DEST_DIR=${1:-${DB_ROOT:-./data/dbs}}
SILVA_URL=${SILVA_URL:-https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz}
PR2_URL=${PR2_URL:-https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_mothur.fasta.gz}
# Optional: set BOLD_URL to download a prepared BOLD fasta (project-specific)
BOLD_URL=${BOLD_URL:-}

mkdir -p "${DEST_DIR}/silva" "${DEST_DIR}/pr2" "${DEST_DIR}/bold_db" "${DEST_DIR}/blast_db"

download_if_missing() {
  local url=$1
  local outpath=$2
  if [ -f "$outpath" ]; then
    echo "Skipping existing file: $outpath"
    return 0
  fi
  echo "Downloading $url -> $outpath"
  curl -fSL "$url" -o "$outpath"
}

echo "Using destination root: ${DEST_DIR}"

SILVA_OUT="${DEST_DIR}/silva/$(basename "$SILVA_URL")"
download_if_missing "$SILVA_URL" "$SILVA_OUT"

PR2_OUT="${DEST_DIR}/pr2/$(basename "$PR2_URL")"
download_if_missing "$PR2_URL" "$PR2_OUT"

if [ -n "$BOLD_URL" ]; then
  BOLD_OUT="${DEST_DIR}/bold_db/$(basename "$BOLD_URL")"
  download_if_missing "$BOLD_URL" "$BOLD_OUT"
else
  echo "BOLD_URL not set; skipping BOLD download. To download BOLD set BOLD_URL and re-run."
fi

echo "Database setup complete. Files placed under: ${DEST_DIR}"
echo "- SILVA: ${SILVA_OUT}"
echo "- PR2:  ${PR2_OUT}"
[ -n "${BOLD_OUT:-}" ] && echo "- BOLD: ${BOLD_OUT}"

echo "You can mount these with DB_ROOT pointing to ${DEST_DIR} in your .env or compose environment."
