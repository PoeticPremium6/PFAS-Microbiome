#!/usr/bin/env bash
# Example configuration file. Copy to config/config.sh and edit.

# Project directories
export PROJECT_DIR="/path/to/project"
export RAW_READ_DIR="$PROJECT_DIR/data/raw"
export QC_DIR="$PROJECT_DIR/results/00_qc"
export TRIMMED_DIR="$PROJECT_DIR/results/01_trimmed"
export CLEAN_READ_DIR="$PROJECT_DIR/results/02_host_removed"
export ASSEMBLY_DIR="$PROJECT_DIR/results/03_assemblies"
export BINNING_DIR="$PROJECT_DIR/results/04_bins"
export CHECKM_DIR="$PROJECT_DIR/results/05_checkm"
export HIGH_QUALITY_MAG_DIR="$PROJECT_DIR/results/06_high_quality_mags"
export TAXONOMY_DIR="$PROJECT_DIR/results/07_taxonomy"
export PHYLOSEQ_DIR="$PROJECT_DIR/results/08_phyloseq_tables"

# Input metadata
export SAMPLE_MANIFEST="$PROJECT_DIR/metadata/sample_manifest.csv"

# Databases/references
export HUMAN_REFERENCE_FASTA="/path/to/GRCh38.fasta"
export HUMAN_BOWTIE2_INDEX="/path/to/bowtie2_indices/GRCh38"
export KRAKEN2_DB="/path/to/kraken2_database"

# Resources
export THREADS=16
export BRACKEN_READ_LENGTH=150

# Quality thresholds
export MIN_MAG_COMPLETENESS=80
export MAX_MAG_CONTAMINATION=10
