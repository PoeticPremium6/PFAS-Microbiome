# Software requirements

This repository uses a mixed workflow involving Bash command-line preprocessing, R-based analysis and figure generation, and MATLAB-based community metabolic reconstruction.

The exact versions used may depend on the local HPC/software environment. The versions below are recommended minimums for reproducibility.

---

## 1. System requirements

Recommended operating system:

- Linux or macOS for preprocessing scripts
- Windows, macOS, or Linux for R figure scripts
- MATLAB-compatible operating system for community metabolic reconstruction

Recommended command-line tools:

```bash
bash >= 4.0
coreutils
awk
sed
grep
gzip
tar
wget or curl
```

---

## 2. Preprocessing software

The preprocessing scripts may require the following tools, depending on which steps are run:

```text
FastQC >= 0.11.9
MultiQC >= 1.13
Trim Galore >= 0.6.7
Cutadapt >= 3.0
Bowtie2 >= 2.4
Samtools >= 1.12
MEGAHIT >= 1.2.9
MetaBAT2 >= 2.15
CheckM >= 1.2
Kraken2 >= 2.1
Bracken >= 2.7
```

Optional, depending on local preprocessing choices:

```text
SPAdes/metaSPAdes
BBMap
BWA
QUAST/metaQUAST
GTDB-Tk
```

External databases/reference resources are not included in this repository and must be downloaded separately:

```text
GRCh38/human genome reference for host-read removal
Kraken2/Bracken taxonomic database
GTDB database, if GTDB-Tk is used
CheckM marker database
```

---

## 3. R requirements

Recommended R version:

```text
R >= 4.2.0
```

Core CRAN packages:

```r
install.packages(c(
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "readr",
  "stringr",
  "forcats",
  "ggplot2",
  "scales",
  "patchwork",
  "cowplot",
  "viridis",
  "vegan",
  "caret",
  "ranger",
  "igraph",
  "tidygraph",
  "ggraph",
  "pheatmap",
  "RColorBrewer",
  "yaml"
))
```

Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "phyloseq",
  "DESeq2",
  "Biostrings",
  "microbiome"
))
```

Optional R packages that may be useful depending on local script edits:

```r
install.packages(c(
  "here",
  "janitor",
  "ggrepel",
  "ggpubr",
  "ggtext",
  "ComplexHeatmap",
  "circlize"
))
```

To save the exact R environment after running the analyses:

```r
writeLines(capture.output(sessionInfo()), "docs/sessionInfo_R.txt")
```

---

## 4. MATLAB requirements

Recommended MATLAB version:

```text
MATLAB R2021b or newer
```

Required MATLAB resources:

```text
COBRA Toolbox
mgPipe
APOLLO gut microbiome genome-scale reconstruction catalogue
Gurobi, IBM CPLEX, or another COBRA-compatible linear programming solver
```

Recommended MATLAB setup commands:

```matlab
initCobraToolbox(false)
changeCobraSolver('gurobi', 'LP')
```

If Gurobi is unavailable, use another solver supported by the COBRA Toolbox and document the solver used.

---

## 5. Suggested conda environment for preprocessing

A minimal conda environment for command-line preprocessing can be created as:

```yaml
name: pfas_microbiome_preprocessing
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - fastqc
  - multiqc
  - cutadapt
  - trim-galore
  - bowtie2
  - samtools
  - megahit
  - metabat2
  - checkm-genome
  - kraken2
  - bracken
  - seqkit
```

Save this as:

```text
environment/preprocessing_environment.yml
```

Then create the environment with:

```bash
conda env create -f environment/preprocessing_environment.yml
conda activate pfas_microbiome_preprocessing
```

---

## 6. Suggested R package setup file

A simple `requirements.R` file can contain:

```r
cran_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "tibble",
  "readr",
  "stringr",
  "forcats",
  "ggplot2",
  "scales",
  "patchwork",
  "cowplot",
  "viridis",
  "vegan",
  "caret",
  "ranger",
  "igraph",
  "tidygraph",
  "ggraph",
  "pheatmap",
  "RColorBrewer",
  "yaml",
  "here",
  "janitor",
  "ggrepel",
  "ggpubr",
  "ggtext"
)

install.packages(setdiff(cran_packages, rownames(installed.packages())))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "phyloseq",
  "DESeq2",
  "Biostrings",
  "microbiome"
)

BiocManager::install(setdiff(bioc_packages, rownames(installed.packages())))
```

---
