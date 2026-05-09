# PFAS gut microbiome manuscript analysis repository

This repository contains the minimal scripts needed to document the preprocessing and microbiome metabolic-model reconstruction workflow used for the manuscript. It is intended for peer reviewers and readers who want to follow the analytical steps, not as a full raw-data archive.

## Workflow overview

1. **Read QC**: inspect raw sequencing quality with FastQC/MultiQC.
2. **Adapter/quality trimming**: trim low-quality bases and adapters.
3. **Human read removal**: map reads to GRCh38 and retain unmapped reads.
4. **Assembly and binning**: assemble metagenomes, bin contigs, and assess MAG quality.
5. **Taxonomic profiling**: classify high-quality MAGs or sample reads with Kraken2/Bracken.
6. **Phyloseq-ready tables**: convert Bracken reports into OTU and taxonomy tables.
7. **Community metabolic reconstruction**: map detected species to APOLLO genome-scale models, reconstruct sample-specific communities with mgPipe, and export reaction/metabolite flux summaries.

## Repository layout

```text
config/                  Example configuration file
metadata/                Example sample manifest
scripts/preprocessing/   Sequencing preprocessing scripts
scripts/modeling/        Community metabolic reconstruction scripts
docs/                    Workflow notes for manuscript methods
data/                    Placeholder only; raw data are not stored here
results/                 Placeholder only; generated outputs are ignored
```

## Quick start

Copy the example config and edit it for your system:

```bash
cp config/config.example.sh config/config.sh
# edit config/config.sh
```

Then run scripts in numerical order, for example:

```bash
bash scripts/preprocessing/00_fastqc_multiqc.sh config/config.sh
bash scripts/preprocessing/01_trim_reads.sh config/config.sh
bash scripts/preprocessing/02_remove_human_reads.sh config/config.sh
bash scripts/preprocessing/03_assemble_and_bin.sh config/config.sh
bash scripts/preprocessing/04_filter_high_quality_mags.sh config/config.sh
bash scripts/preprocessing/05_taxonomy_kraken_bracken.sh config/config.sh
bash scripts/preprocessing/06_bracken_to_phyloseq_tables.sh config/config.sh
```

Community model reconstruction is run in MATLAB:

```matlab
run('scripts/modeling/07_community_reconstruction.m')
```

Before running the MATLAB script, edit the configuration block at the top of `07_community_reconstruction.m` or pass equivalent paths interactively.

## Data availability note

