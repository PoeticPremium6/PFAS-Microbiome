#!/usr/bin/env bash
set -euo pipefail

CONFIG=${1:-config/config.sh}
source "$CONFIG"

mkdir -p "$TAXONOMY_DIR"

if [ ! -f "$KRAKEN2_DB/database${BRACKEN_READ_LENGTH:-150}mers.kmer_distrib" ]; then
  bracken-build \
    -d "$KRAKEN2_DB" \
    -t "${THREADS:-8}" \
    -k 35 \
    -l "${BRACKEN_READ_LENGTH:-150}"
fi

for fasta in "$HIGH_QUALITY_MAG_DIR"/*.fa "$HIGH_QUALITY_MAG_DIR"/*.fasta; do
  [ -f "$fasta" ] || continue
  base=$(basename "$fasta")
  base=${base%.*}

  echo "Classifying $base"
  kraken2 \
    --db "$KRAKEN2_DB" \
    --threads "${THREADS:-8}" \
    --report "$TAXONOMY_DIR/${base}_kraken_report.tsv" \
    --output "$TAXONOMY_DIR/${base}_kraken.out" \
    "$fasta"

  bracken \
    -d "$KRAKEN2_DB" \
    -i "$TAXONOMY_DIR/${base}_kraken_report.tsv" \
    -o "$TAXONOMY_DIR/${base}_bracken_species.tsv" \
    -r "${BRACKEN_READ_LENGTH:-150}" \
    -l S
done

echo "Taxonomic classification complete: $TAXONOMY_DIR"
