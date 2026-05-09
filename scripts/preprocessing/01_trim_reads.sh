#!/usr/bin/env bash
set -euo pipefail

CONFIG=${1:-config/config.sh}
source "$CONFIG"

mkdir -p "$TRIMMED_DIR"

# sample_manifest.csv columns required: sample_id,read1,read2
awk -F',' 'NR>1 {print $1","$2","$3}' "$SAMPLE_MANIFEST" | while IFS=',' read -r sample_id read1 read2; do
  echo "Trimming $sample_id"
  cutadapt \
    --cores "${THREADS:-8}" \
    --quality-cutoff 10 \
    --minimum-length 30 \
    -o "$TRIMMED_DIR/${sample_id}_R1.trimmed.fastq.gz" \
    -p "$TRIMMED_DIR/${sample_id}_R2.trimmed.fastq.gz" \
    "$read1" "$read2" \
    > "$TRIMMED_DIR/${sample_id}.cutadapt.log"
done

echo "Trimming complete: $TRIMMED_DIR"
