#!/usr/bin/env bash
set -euo pipefail

CONFIG=${1:-config/config.sh}
source "$CONFIG"

mkdir -p "$ASSEMBLY_DIR" "$BINNING_DIR"

awk -F',' 'NR>1 {print $1}' "$SAMPLE_MANIFEST" | while read -r sample_id; do
  r1="$CLEAN_READ_DIR/${sample_id}_R1.clean.fastq.gz"
  r2="$CLEAN_READ_DIR/${sample_id}_R2.clean.fastq.gz"
  assembly_out="$ASSEMBLY_DIR/${sample_id}"
  binning_out="$BINNING_DIR/${sample_id}"

  echo "Assembling $sample_id"
  megahit \
    -1 "$r1" \
    -2 "$r2" \
    -o "$assembly_out" \
    -t "${THREADS:-8}"

  echo "Binning $sample_id"
  mkdir -p "$binning_out"
  metabat2 \
    -i "$assembly_out/final.contigs.fa" \
    -o "$binning_out/${sample_id}_bin" \
    -t "${THREADS:-8}" \
    -m 1500
done

echo "Assembly and binning complete."
