# ==========================================================
# Build Phyloseq object and prepare PFAS metadata
# ==========================================================

library(phyloseq)
library(tidyverse)

# -------------------
# File paths
# -------------------
otu_file  <- file.path("data", "processed", "otu_table.tsv")
tax_file  <- file.path("data", "processed", "taxonomy_table.tsv")
meta_file <- file.path("data", "metadata", "PFAS_Metadata.csv")
sample_metadata <- file.path("data", "metadata", "PFAS_Metadata.csv")

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
head(as.data.frame(sample_data(physeq))[, 1:6])

# -------------------
# 5. Subset to PFAS_Exposed + Reference
# -------------------
physeq_subset <- subset_samples(
  physeq,
  Condition %in% c("PFAS_Exposed", "Reference")
)

# -------------------
# 6. Compute average elimination rate constant
# -------------------
pfas_columns <- c(
  "k_PFOA",
  "k_PFPeS",
  "k_PFHxS",
  "k_PFHpS",
  "k_LPFOS",
  "k_PFOS_MP1",
  "k_PFOS_MP345",
  "k_PFOS_MP26"
)

available_cols <- intersect(
  pfas_columns,
  colnames(sample_data(physeq_subset))
)

if (length(available_cols) == 0) {
  stop("No PFAS elimination rate constant columns found in metadata.")
} else if (length(available_cols) < length(pfas_columns)) {
  warning(
    "Some PFAS columns missing. Using only: ",
    paste(available_cols, collapse = ", ")
  )
}

sample_data(physeq_subset)$average_k <- rowMeans(
  as.data.frame(sample_data(physeq_subset))[, available_cols, drop = FALSE],
  na.rm = TRUE
)

#################################
#################################
# Figure 1
# A/B Order-level relative abundance
# C Beta-dispersion distance to centroid
#
# Final version:
# - Removed Bray-Curtis PCoA panel because it is redundant
# - A/B have biological titles: PFAS-exposed / Reference
# - C has no title
# - All panels have visible standalone labels: A, B, C
# - Saves individual panels and combined portrait figure
#################################
#################################

# -------------------
# 0. Load Libraries
# -------------------
library(phyloseq)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(vegan)
library(viridis)
library(patchwork)
library(tibble)
library(cowplot)

# -------------------
# 1. File paths
# -------------------
otu_file <- file.path("data", "processed", "otu_table.tsv")
tax_file <- file.path("data", "processed", "taxonomy_table.tsv")
meta_file <- file.path("data", "metadata", "PFAS_Metadata.csv")

fig_dir <- file.path("results", "figures")
csv_dir <- file.path("results", "figures")

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
if (!dir.exists(csv_dir)) dir.create(csv_dir, recursive = TRUE)

# -------------------
# 2. Load OTU Table
# -------------------
otu_df <- read.delim(
  otu_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

otu_mat <- as.matrix(otu_df)
otu_table_ps <- otu_table(otu_mat, taxa_are_rows = TRUE)
taxa_names(otu_table_ps) <- trimws(taxa_names(otu_table_ps))

# -------------------
# 3. Load Taxonomy Table
# -------------------
tax_df <- read.delim(
  tax_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dup_taxa <- tax_df[, 1][duplicated(tax_df[, 1])]

if (length(dup_taxa) > 0) {
  cat("Warning: Found duplicate taxa names. Making them unique.\n")
  tax_df[, 1] <- make.unique(tax_df[, 1])
}

rownames(tax_df) <- trimws(tax_df[, 1])
tax_df <- tax_df[, -1]
tax_mat <- as.matrix(tax_df)
tax_table_ps <- tax_table(tax_mat)

# -------------------
# 4. Load Sample Metadata
# -------------------
sample_metadata <- read.csv(
  meta_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

sample_data_ps <- sample_data(sample_metadata)

# -------------------
# 5. Keep only matching taxa
# -------------------
common_taxa <- intersect(
  taxa_names(otu_table_ps),
  rownames(tax_table_ps)
)

otu_table_ps <- prune_taxa(common_taxa, otu_table_ps)
tax_table_ps <- tax_table_ps[common_taxa, ]

# -------------------
# 6. Build Phyloseq Object
# -------------------
physeq <- phyloseq(
  otu_table_ps,
  tax_table_ps,
  sample_data_ps
)

# -------------------
# 7. Filter taxa and samples
# -------------------
physeq_filtered <- filter_taxa(
  physeq,
  function(x) sum(x > 0) >= 3,
  TRUE
)

physeq_filtered <- prune_samples(
  sample_sums(physeq_filtered) > 1000,
  physeq_filtered
)

covariates <- c("Age", "Sex", "weight", "BMI")
meta_df <- as.data.frame(sample_data(physeq_filtered))

available_covariates <- intersect(covariates, colnames(meta_df))

if (length(available_covariates) > 0) {
  keep_samples <- rownames(meta_df)[
    complete.cases(meta_df[, available_covariates, drop = FALSE])
  ]
  physeq_filtered <- prune_samples(keep_samples, physeq_filtered)
}

meta_df <- as.data.frame(sample_data(physeq_filtered))

sample_data(physeq_filtered)$Condition <- factor(
  sample_data(physeq_filtered)$Condition,
  levels = c("PFAS_Exposed", "Reference")
)

# -------------------
# 8. DESeq2 VST with covariates
# Kept from original workflow
# -------------------
covariates <- intersect(covariates, colnames(meta_df))

design_formula <- if (length(covariates) > 0) {
  as.formula(
    paste("~", paste(covariates, collapse = " + "), "+ Condition")
  )
} else {
  ~ Condition
}

cat("DESeq2 design formula:", deparse(design_formula), "\n")

dds <- phyloseq_to_deseq2(physeq_filtered, design_formula)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds)

vst <- varianceStabilizingTransformation(dds, blind = FALSE)

vst_physeq <- phyloseq(
  otu_table(assay(vst), taxa_are_rows = TRUE),
  tax_table(physeq_filtered),
  sample_data(physeq_filtered)
)

# -------------------
# 9. Relative Abundance and Order Aggregation
# -------------------
physeq_rel_abund <- transform_sample_counts(
  physeq_filtered,
  function(x) x / sum(x)
)

physeq_rel_abund <- prune_taxa(
  taxa_sums(physeq_rel_abund) > 0,
  physeq_rel_abund
)

physeq_Order_rel_abund <- tax_glom(
  physeq_rel_abund,
  taxrank = "Order"
)

physeq_case <- subset_samples(
  physeq_Order_rel_abund,
  Condition == "PFAS_Exposed"
)

physeq_control <- subset_samples(
  physeq_Order_rel_abund,
  Condition == "Reference"
)

# -------------------
# 10. Shared themes
# -------------------
base_theme <- theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    strip.text = element_text(size = 14, face = "bold")
  )

# -------------------
# 11. Helper function to add panel labels
# -------------------
add_panel_label <- function(plot_obj, label_text) {
  cowplot::ggdraw(plot_obj) +
    cowplot::draw_label(
      label_text,
      x = 0.015,
      y = 0.985,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = 24
    )
}

# -------------------
# 12. Plot Relative Abundance Function
# A/B keep biological titles
# -------------------
plot_rel_abund <- function(physeq_obj, title_text) {
  
  plot_bar(physeq_obj, fill = "Order") +
    geom_bar(
      stat = "identity",
      position = "stack",
      width = 0.85
    ) +
    scale_fill_viridis(
      discrete = TRUE,
      option = "magma",
      direction = -1
    ) +
    base_theme +
    labs(
      title = title_text,
      x = "Sample",
      y = "Relative abundance",
      fill = "Order"
    ) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text.x = element_text(
        size = 10,
        face = "bold",
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.text.y = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11, face = "bold"),
      plot.margin = margin(14, 10, 10, 18)
    )
}

p_case_raw <- plot_rel_abund(
  physeq_case,
  "PFAS-exposed"
)

p_control_raw <- plot_rel_abund(
  physeq_control,
  "Reference"
)

p_case <- add_panel_label(p_case_raw, "A")
p_control <- add_panel_label(p_control_raw, "B")

print(p_case)
print(p_control)

ggsave(
  file.path(fig_dir, "Figure1A_taxa_rel_abund_case.png"),
  plot = p_case,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "Figure1B_taxa_rel_abund_control.png"),
  plot = p_control,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

# -------------------
# 13. Export Mean Relative Abundance per Order
# -------------------
calc_mean_rel_abund <- function(physeq_obj) {
  
  otu_df <- as.data.frame(otu_table(physeq_obj))
  order_tax <- as.character(tax_table(physeq_obj)[, "Order"])
  mean_abund <- rowMeans(otu_df)
  
  data.frame(
    Order = order_tax,
    Mean_Relative_Abundance = mean_abund
  )
}

write.csv(
  calc_mean_rel_abund(physeq_case),
  file.path(csv_dir, "Figure1_Order_relative_abundance_case.csv"),
  row.names = FALSE
)

write.csv(
  calc_mean_rel_abund(physeq_control),
  file.path(csv_dir, "Figure1_Order_relative_abundance_control.csv"),
  row.names = FALSE
)

# ==========================================================
# Figure 1C: beta-dispersion only
# Bray-Curtis PCoA removed from Figure 1 because it is redundant
# ==========================================================

# -------------------
# 14. Prepare beta-diversity object
# -------------------
physeq_beta <- subset_samples(
  physeq_filtered,
  Condition %in% c("PFAS_Exposed", "Reference")
)

physeq_beta <- prune_taxa(
  taxa_sums(physeq_beta) > 0,
  physeq_beta
)

physeq_beta <- prune_samples(
  sample_sums(physeq_beta) > 0,
  physeq_beta
)

sample_data(physeq_beta)$Condition <- factor(
  sample_data(physeq_beta)$Condition,
  levels = c("PFAS_Exposed", "Reference")
)

physeq_beta_rel <- transform_sample_counts(
  physeq_beta,
  function(x) x / sum(x)
)

# -------------------
# 15. Bray-Curtis distance
# -------------------
bray_dist <- phyloseq::distance(
  physeq_beta_rel,
  method = "bray"
)

meta_beta <- data.frame(
  sample_data(physeq_beta_rel),
  check.names = FALSE
)

meta_beta$SampleID <- rownames(meta_beta)

meta_beta <- meta_beta[labels(bray_dist), , drop = FALSE]
meta_beta <- as.data.frame(meta_beta)

stopifnot(all(rownames(meta_beta) == labels(bray_dist)))

# -------------------
# 16. PERMANOVA
# Kept for reporting, but PCoA panel is no longer plotted in Figure 1
# -------------------
set.seed(123)

permanova_res <- vegan::adonis2(
  bray_dist ~ Condition,
  data = meta_beta,
  permutations = 999
)

print(permanova_res)

permanova_R2 <- permanova_res$R2[1]
permanova_p <- permanova_res$`Pr(>F)`[1]

permanova_label <- paste0(
  "PERMANOVA: R² = ",
  signif(permanova_R2, 2),
  ", p = ",
  signif(permanova_p, 3)
)

# -------------------
# 17. Beta-dispersion
# -------------------
dispersion <- vegan::betadisper(
  bray_dist,
  group = meta_beta$Condition,
  type = "centroid"
)

set.seed(123)

dispersion_perm <- vegan::permutest(
  dispersion,
  permutations = 999
)

print(anova(dispersion))
print(dispersion_perm)

dispersion_p <- dispersion_perm$tab$`Pr(>F)`[1]

dispersion_label <- paste0(
  "Betadisper: p = ",
  signif(dispersion_p, 3)
)

disp_df <- data.frame(
  SampleID = names(dispersion$distances),
  Distance_to_centroid = dispersion$distances,
  Condition = meta_beta[names(dispersion$distances), "Condition"]
)

disp_summary <- disp_df %>%
  group_by(Condition) %>%
  summarise(
    n = n(),
    mean_distance = mean(Distance_to_centroid, na.rm = TRUE),
    sd_distance = sd(Distance_to_centroid, na.rm = TRUE),
    median_distance = median(Distance_to_centroid, na.rm = TRUE),
    .groups = "drop"
  )

print(disp_summary)

write.csv(
  disp_df,
  file.path(csv_dir, "Figure1C_BetaDispersion_DistanceToCentroid_Data.csv"),
  row.names = FALSE
)

write.csv(
  disp_summary,
  file.path(csv_dir, "Figure1C_BetaDispersion_Summary.csv"),
  row.names = FALSE
)

# -------------------
# 18. Distinct colours for Figure 1C
# PFAS = purple; Reference = teal
# -------------------
condition_colors <- c(
  "PFAS_Exposed" = "#6A1B9A",
  "Reference" = "#00796B"
)

condition_fills <- c(
  "PFAS_Exposed" = "#CE93D8",
  "Reference" = "#80CBC4"
)

# -------------------
# 19. Figure 1C: beta-dispersion distance to centroid
# No plot title; panel label added separately
# -------------------
p_1C_raw <- ggplot(
  disp_df,
  aes(x = Condition, y = Distance_to_centroid, fill = Condition)
) +
  geom_boxplot(
    width = 0.55,
    alpha = 0.9,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.6
  ) +
  geom_jitter(
    aes(color = Condition),
    width = 0.12,
    size = 3.4,
    alpha = 0.9
  ) +
  scale_fill_manual(values = condition_fills) +
  scale_color_manual(values = condition_colors) +
  annotate(
    "label",
    x = 1.5,
    y = max(disp_df$Distance_to_centroid, na.rm = TRUE) * 1.08,
    label = dispersion_label,
    size = 4.1,
    fontface = "bold",
    label.size = 0.25,
    fill = "white"
  ) +
  base_theme +
  labs(
    title = NULL,
    x = NULL,
    y = "Bray-Curtis distance to group centroid"
  ) +
  theme(
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 13, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold"),
    legend.position = "none",
    plot.margin = margin(14, 10, 10, 18)
  )

p_1C <- add_panel_label(p_1C_raw, "C")

print(p_1C)

ggsave(
  file.path(fig_dir, "Figure1C_BetaDispersion_DistanceToCentroid.png"),
  plot = p_1C,
  width = 6.5,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "Figure1C_BetaDispersion_DistanceToCentroid.pdf"),
  plot = p_1C,
  width = 6.5,
  height = 6,
  bg = "white"
)

# ==========================================================
# 20. Combined Figure 1
# Portrait layout: A, B, C stacked vertically
# A/B have titles; C has no title
# A/B/C labels are added manually and preserved
# ==========================================================

combined_figure1 <- 
  p_case /
  p_control /
  p_1C +
  plot_layout(
    heights = c(1.15, 1.00, 0.90)
  )

print(combined_figure1)

ggsave(
  file.path(fig_dir, "Figure1_Combined_ABC_Portrait.png"),
  plot = combined_figure1,
  width = 9,
  height = 18,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "Figure1_Combined_ABC_Portrait.pdf"),
  plot = combined_figure1,
  width = 9,
  height = 18,
  bg = "white"
)

cat("Figure 1A-C plots printed and saved.\n")
cat("PERMANOVA calculated but PCoA panel removed from Figure 1: ", permanova_label, "\n")
cat("Betadisper label shown in Figure 1C: ", dispersion_label, "\n")


#######################
#######################
# ==========================================================
