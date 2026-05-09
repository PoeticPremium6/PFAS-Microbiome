#!/usr/bin/env bash
set -euo pipefail

CONFIG=${1:-config/config.sh}
source "$CONFIG"

mkdir -p "$CHECKM_DIR/all_bins" "$CHECKM_DIR/checkm_output" "$HIGH_QUALITY_MAG_DIR"

# Flatten bins from all sample folders into one CheckM input directory.
find "$BINNING_DIR" -type f \( -name "*.fa" -o -name "*.fasta" \) -exec cp {} "$CHECKM_DIR/all_bins/" \;

checkm lineage_wf \
  -x fa \
  --reduced_tree \
  -t "${THREADS:-8}" \
  "$CHECKM_DIR/all_bins" \
  "$CHECKM_DIR/checkm_output"

checkm qa \
  "$CHECKM_DIR/checkm_output/lineage.ms" \
  "$CHECKM_DIR/checkm_output" \
  -f "$CHECKM_DIR/qa.tsv" \
  --tab_table \
  -o 2

awk -F '\t' -v minc="${MIN_MAG_COMPLETENESS:-80}" -v maxc="${MAX_MAG_CONTAMINATION:-10}" \
  'NR > 1 && $5 >= minc && $6 <= maxc { print $1".fa" }' "$CHECKM_DIR/qa.tsv" | while read -r binfa; do
    if [ -f "$CHECKM_DIR/all_bins/$binfa" ]; then
      cp "$CHECKM_DIR/all_bins/$binfa" "$HIGH_QUALITY_MAG_DIR/"
    fi
  done

echo "High-quality MAG filtering complete: $HIGH_QUALITY_MAG_DIR"
