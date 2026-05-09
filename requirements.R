# R package requirements for PFAS microbiome manuscript analyses

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

missing_cran <- setdiff(cran_packages, rownames(installed.packages()))
if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "phyloseq",
  "DESeq2",
  "Biostrings",
  "microbiome"
)

missing_bioc <- setdiff(bioc_packages, rownames(installed.packages()))
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc)
}

writeLines(capture.output(sessionInfo()), "docs/sessionInfo_R.txt")
