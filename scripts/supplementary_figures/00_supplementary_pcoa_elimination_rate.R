############################################################
# 00_supplementary_pcoa_elimination_rate.R
#
# Supplementary ordination panel used during manuscript development.
#
############################################################

# ==========================================================
# Build Phyloseq object and perform PCoA with PFAS metadata
# ==========================================================

library(phyloseq)
library(tidyverse)

# -------------------
# File paths
# -------------------
otu_file  <- file.path("data", "processed", "otu_table.tsv")
tax_file  <- file.path("data", "processed", "taxonomy_table.tsv")
meta_file <- file.path("data", "metadata", "PFAS_Metadata.csv")

# -------------------
# 1. Import OTU table
# -------------------
otu_df <- read.delim(otu_file, header = TRUE, check.names = FALSE)
otu_df$Taxa <- trimws(otu_df$Taxa)
rownames(otu_df) <- make.unique(otu_df$Taxa)
otu_df$Taxa <- NULL
otu <- otu_table(as.matrix(otu_df), taxa_are_rows = TRUE)

# -------------------
# 2. Import taxonomy
# -------------------
tax_df <- read.delim(tax_file, header = TRUE, check.names = FALSE)
tax_df$Taxa <- trimws(tax_df$Taxa)
rownames(tax_df) <- make.unique(tax_df$Taxa)
tax_df$Taxa <- NULL
tax <- tax_table(as.matrix(tax_df))

# -------------------
# 3. Import metadata
# -------------------
meta_df <- read.csv(meta_file, header = TRUE, check.names = FALSE)

# Align rownames with sample IDs
meta_df$Sample_Name <- trimws(meta_df$Sample_Name)
rownames(meta_df) <- meta_df$Sample_Name
meta_df$Sample_Name <- NULL

samp <- sample_data(meta_df)

# -------------------
# 4. Build phyloseq object
# -------------------
physeq <- phyloseq(otu, tax, samp)

# Quick check
sample_variables(physeq)[1:10]
head(as.data.frame(sample_data(physeq))[,1:6])

# -------------------
# 5. Subset to PFAS_Exposed + Reference
# -------------------
physeq_subset <- subset_samples(physeq, Condition %in% c("PFAS_Exposed", "Reference"))

# -------------------
# 6. Compute average elimination rate constant
# -------------------
pfas_columns <- c("k_PFOA", "k_PFPeS", "k_PFHxS", "k_PFHpS",
                  "k_LPFOS", "k_PFOS_MP1", "k_PFOS_MP345", "k_PFOS_MP26")

available_cols <- intersect(pfas_columns, colnames(sample_data(physeq_subset)))
if (length(available_cols) == 0) {
  stop("No PFAS elimination rate constant columns found in metadata.")
} else if (length(available_cols) < length(pfas_columns)) {
  warning("Some PFAS columns missing. Using only: ", paste(available_cols, collapse = ", "))
}

sample_data(physeq_subset)$average_k <- rowMeans(
  as.data.frame(sample_data(physeq_subset))[, available_cols, drop = FALSE],
  na.rm = TRUE
)

# -------------------
# 7. Ordination (Bray-Curtis PCoA)
# -------------------
ordination_bray <- ordinate(physeq_subset, method = "PCoA", distance = "bray")
axis_labels <- ordination_bray$values$Relative_eig * 100

# -------------------
# 8. Plot
# -------------------
p <- plot_ordination(physeq_subset, ordination_bray) +
  geom_point(
    aes(fill = average_k, color = Condition),
    shape = 21, size = 6
  ) +
  scale_fill_gradient(low = "gray", high = "purple", name = "Avg Elimination Rate") +
  scale_color_manual(values = c("PFAS_Exposed" = "purple", "Reference" = "black"),
                     name = "Condition") +
  guides(
    fill = guide_colorbar(title = "Elimination Rate", barwidth = 2, barheight = 4, title.position = "top"),
    color = guide_legend(title = "Condition", override.aes = list(shape = 21, size = 6))
  ) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(axis_labels[1], 2), "%)"),
    y = paste0("PC2 (", round(axis_labels[2], 2), "%)")
  ) +
  theme(
    plot.title   = element_text(size = 16, face = "bold"),
    axis.title   = element_text(size = 14, face = "bold"),
    axis.text    = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12)
  )

# -------------------
# 9. Save Plot
# -------------------
ggsave(
  filename = "Figure2B_Bray_Curtis_PCoA_with_Elimination_Rate.png",
  plot = p,
  path = file.path("results", "supplementary_figures"),
  width = 8, height = 6, dpi = 300, bg = "white"
)

# -------------------
# 10. Print Plot
# -------------------
print(p)

# =========================================
