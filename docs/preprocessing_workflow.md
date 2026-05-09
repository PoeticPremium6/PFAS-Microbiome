# Preprocessing workflow notes

This document mirrors the numbered scripts in `scripts/preprocessing/`.

## 00. QC
Raw reads are inspected using FastQC and summarized using MultiQC.

## 01. Trimming
Adapters and low-quality sequence are removed. The example script uses `cutadapt`, but equivalent trimming software can be substituted if parameters are reported.

## 02. Human read removal
Trimmed reads are aligned to a human reference genome using Bowtie2. Only unmapped read pairs are retained for downstream metagenomic analysis.

## 03. Assembly and binning
Cleaned reads are assembled per sample and metagenomic bins are generated from the resulting contigs. The example script uses MEGAHIT and MetaBAT2.

## 04. MAG quality filtering
MAG quality is assessed with CheckM. Bins passing the completeness and contamination thresholds in the config file are copied to the high-quality MAG directory.

## 05. Taxonomic classification
MAGs or read-derived assemblies are classified with Kraken2. Bracken is then used to estimate species-level abundances.

## 06. Phyloseq tables
Bracken outputs are converted into an OTU/count table and taxonomy table for downstream R/phyloseq analysis.

## 07. Community metabolic reconstruction
Detected species are matched to APOLLO genome-scale reconstructions. Sample-specific community models are reconstructed with mgPipe and simulated to export reaction and metabolite flux summaries.
