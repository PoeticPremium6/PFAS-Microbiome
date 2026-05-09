#!/usr/bin/env bash
set -euo pipefail

CONFIG=${1:-config/config.sh}
source "$CONFIG"

mkdir -p "$QC_DIR/fastqc" "$QC_DIR/multiqc"

# Run FastQC on all raw FASTQ files.
find "$RAW_READ_DIR" -type f \( -name "*.fastq.gz" -o -name "*.fq.gz" -o -name "*.fastq" -o -name "*.fq" \) \
  -print0 | xargs -0 -n 1 -P "${THREADS:-8}" fastqc -o "$QC_DIR/fastqc"

# Summarize QC results.
multiqc "$QC_DIR/fastqc" -o "$QC_DIR/multiqc"

echo "QC complete: $QC_DIR"
