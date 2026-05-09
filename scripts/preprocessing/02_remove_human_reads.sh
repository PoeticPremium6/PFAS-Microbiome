#!/usr/bin/env bash
set -euo pipefail

CONFIG=${1:-config/config.sh}
source "$CONFIG"

mkdir -p "$CLEAN_READ_DIR"

# Build the Bowtie2 index once if it does not already exist.
if [ ! -f "${HUMAN_BOWTIE2_INDEX}.1.bt2" ] && [ ! -f "${HUMAN_BOWTIE2_INDEX}.1.bt2l" ]; then
  bowtie2-build "$HUMAN_REFERENCE_FASTA" "$HUMAN_BOWTIE2_INDEX"
fi

awk -F',' 'NR>1 {print $1}' "$SAMPLE_MANIFEST" | while read -r sample_id; do
  r1="$TRIMMED_DIR/${sample_id}_R1.trimmed.fastq.gz"
  r2="$TRIMMED_DIR/${sample_id}_R2.trimmed.fastq.gz"

  echo "Removing human reads from $sample_id"
  bowtie2 \
    -x "$HUMAN_BOWTIE2_INDEX" \
    -1 "$r1" \
    -2 "$r2" \
    --very-sensitive \
    --threads "${THREADS:-8}" \
    --un-conc-gz "$CLEAN_READ_DIR/${sample_id}.host_removed_R%.fastq.gz" \
    -S /dev/null \
    2> "$CLEAN_READ_DIR/${sample_id}.bowtie2_host_removal.log"

  # Bowtie2 writes R1/R2 as *_R1.fastq.gz and *_R2.fastq.gz using the % placeholder.
  mv "$CLEAN_READ_DIR/${sample_id}.host_removed_R1.fastq.gz" "$CLEAN_READ_DIR/${sample_id}_R1.clean.fastq.gz" 2>/dev/null || true
  mv "$CLEAN_READ_DIR/${sample_id}.host_removed_R2.fastq.gz" "$CLEAN_READ_DIR/${sample_id}_R2.clean.fastq.gz" 2>/dev/null || true
done

echo "Human read removal complete: $CLEAN_READ_DIR"
