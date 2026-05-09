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
# Supplementary Figure S1
# Species-level DESeq2 log2 fold-changes
# PFAS-exposed vs Reference
#
# Landscape-page optimized version:
# - Designed for one full-width landscape supplementary page
# - No confidence intervals
# - Adds zero reference line
# - Adds nominal p-value label for each bar
# - Saves individual subplots
# - Saves combined A/B figure optimized for landscape layout
# - Exports plotted taxa with p-value and adjusted p-value
# =========================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(patchwork)
library(tibble)
library(stringr)

# -------------------
# Output directory
# -------------------
fig_dir <- file.path("results", "supplementary_figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# -------------------
# 1. Build DESeq2 object with covariates
# -------------------
covariates <- c("Age", "Sex", "weight", "BMI")

meta_df <- as.data.frame(sample_data(physeq_filtered))
covariates <- intersect(covariates, colnames(meta_df))

design_formula <- if (length(covariates) > 0) {
  as.formula(paste("~", paste(covariates, collapse = " + "), "+ Condition"))
} else {
  ~ Condition
}

cat("DESeq2 design formula:", deparse(design_formula), "\n")

dds <- phyloseq_to_deseq2(physeq_filtered, design_formula)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds)

# -------------------
# 2. Extract results: PFAS_Exposed vs Reference
# -------------------
res <- results(
  dds,
  contrast = c("Condition", "PFAS_Exposed", "Reference")
)

res_df <- as.data.frame(res) %>%
  rownames_to_column("OTU")

# -------------------
# 3. Add p-value labels
# -------------------
format_p <- function(p) {
  case_when(
    is.na(p) ~ "p=NA",
    p < 0.001 ~ "p<0.001",
    p < 0.01 ~ paste0("p=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("p=", formatC(p, format = "f", digits = 2))
  )
}

format_padj <- function(p) {
  case_when(
    is.na(p) ~ "padj=NA",
    p < 0.001 ~ "padj<0.001",
    p < 0.01 ~ paste0("padj=", formatC(p, format = "f", digits = 3)),
    TRUE ~ paste0("padj=", formatC(p, format = "f", digits = 2))
  )
}

res_df <- res_df %>%
  mutate(
    p_label = format_p(pvalue),
    padj_label = format_padj(padj),
    Significance = case_when(
      !is.na(padj) & padj < 0.05 ~ "FDR < 0.05",
      !is.na(pvalue) & pvalue < 0.05 ~ "Nominal p < 0.05",
      TRUE ~ "Not significant"
    )
  )

cat("\nDESeq2 significance summary:\n")
print(table(res_df$Significance, useNA = "ifany"))

cat("\nNumber FDR significant taxa:\n")
print(sum(!is.na(res_df$padj) & res_df$padj < 0.05))

cat("\nNumber nominal p < 0.05 taxa:\n")
print(sum(!is.na(res_df$pvalue) & res_df$pvalue < 0.05))

# -------------------
# 4. Merge taxonomy
# -------------------
tax_df <- as.data.frame(tax_table(physeq_filtered))
tax_df$OTU <- rownames(tax_df)

res_df <- res_df %>%
  left_join(tax_df, by = "OTU")

res_df$Species <- ifelse(
  is.na(res_df$Species) | res_df$Species == "",
  res_df$OTU,
  res_df$Species
)

res_df <- res_df %>%
  filter(!is.na(log2FoldChange))

# -------------------
# 5. Identify higher- vs lower-abundance species
# Based on mean normalized counts
# -------------------
norm_counts <- counts(dds, normalized = TRUE)
mean_counts <- rowMeans(norm_counts)

rare_threshold <- quantile(mean_counts, 0.25, na.rm = TRUE)

lower_abundance_otus <- names(mean_counts)[mean_counts <= rare_threshold]
higher_abundance_otus <- names(mean_counts)[mean_counts > rare_threshold]

res_lower <- res_df %>%
  filter(OTU %in% lower_abundance_otus)

res_higher <- res_df %>%
  filter(OTU %in% higher_abundance_otus)

# -------------------
# 6. Select top species by absolute log2 fold-change
# -------------------
top_n <- 10

top_higher <- res_higher %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = top_n) %>%
  arrange(log2FoldChange) %>%
  mutate(
    Species_clean = str_replace_all(Species, "_", " "),
    Species_clean = str_wrap(Species_clean, width = 32),
    Species_clean = factor(Species_clean, levels = Species_clean)
  )

top_lower <- res_lower %>%
  arrange(desc(abs(log2FoldChange))) %>%
  slice_head(n = top_n) %>%
  arrange(log2FoldChange) %>%
  mutate(
    Species_clean = str_replace_all(Species, "_", " "),
    Species_clean = str_wrap(Species_clean, width = 32),
    Species_clean = factor(Species_clean, levels = Species_clean)
  )

# -------------------
# 7. Export plotted taxa with statistics
# -------------------
write.csv(
  top_higher,
  file.path(fig_dir, "SuppFigS1A_higher_abundance_species_DESeq2.csv"),
  row.names = FALSE
)

write.csv(
  top_lower,
  file.path(fig_dir, "SuppFigS1B_lower_abundance_species_DESeq2.csv"),
  row.names = FALSE
)

# -------------------
# 8. Colours
# -------------------
bar_fill <- "#6A1B9A"
bar_outline <- "black"
zero_line <- "grey35"

# -------------------
# 9. Landscape-optimized shared plotting theme
# -------------------
supp_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 10.5, face = "bold", color = "black", lineheight = 0.9),
    axis.title.x = element_text(size = 13, face = "bold", color = "black", margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, face = "bold", color = "black", margin = margin(r = 8)),
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey88", linewidth = 0.35),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey25"),
    plot.margin = margin(8, 18, 8, 8)
  )

# -------------------
# 10. Shared x-axis limits across both panels
# This helps the landscape figure look balanced
# -------------------
all_top <- bind_rows(top_higher, top_lower)

x_min_global <- floor(min(all_top$log2FoldChange, na.rm = TRUE) - 1.4)
x_max_global <- ceiling(max(all_top$log2FoldChange, na.rm = TRUE) + 1.4)

# Keep axis visually balanced around zero
x_abs <- max(abs(c(x_min_global, x_max_global)), na.rm = TRUE)
x_limits <- c(-x_abs, x_abs)

# -------------------
# 11. Plot function
# -------------------
plot_deseq_lfc <- function(df, y_lab) {
  
  label_offset <- 0.35
  
  df <- df %>%
    mutate(
      label_x = ifelse(
        log2FoldChange >= 0,
        log2FoldChange + label_offset,
        log2FoldChange - label_offset
      ),
      label_hjust = ifelse(log2FoldChange >= 0, 0, 1)
    )
  
  ggplot(
    df,
    aes(
      y = Species_clean,
      x = log2FoldChange
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = zero_line,
      linewidth = 0.65
    ) +
    geom_col(
      width = 0.62,
      fill = bar_fill,
      color = bar_outline,
      linewidth = 0.35,
      alpha = 0.92
    ) +
    geom_text(
      aes(
        x = label_x,
        label = p_label,
        hjust = label_hjust
      ),
      size = 3.1,
      fontface = "bold",
      color = "black"
    ) +
    scale_x_continuous(
      limits = x_limits,
      breaks = pretty(x_limits, n = 5),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    supp_theme +
    labs(
      x = "DESeq2 log2 fold-change: PFAS-exposed vs Reference",
      y = y_lab,
      subtitle = ""
    ) +
    coord_cartesian(clip = "off")
}

# -------------------
# 12. Higher-abundance species plot
# -------------------
p_higher <- plot_deseq_lfc(
  df = top_higher,
  y_lab = "Higher-abundance species\n(top 75% by mean normalized count)"
)

print(p_higher)

ggsave(
  file.path(fig_dir, "SuppFigS1A_HigherAbundance_Species_log2FC_pvalues_landscape.png"),
  plot = p_higher,
  width = 9.5,
  height = 5.4,
  dpi = 300,
  bg = "white"
)

# -------------------
# 13. Lower-abundance species plot
# -------------------
p_lower <- plot_deseq_lfc(
  df = top_lower,
  y_lab = "Lower-abundance species\n(bottom 25% by mean normalized count)"
)

print(p_lower)

ggsave(
  file.path(fig_dir, "SuppFigS1B_LowerAbundance_Species_log2FC_pvalues_landscape.png"),
  plot = p_lower,
  width = 9.5,
  height = 5.4,
  dpi = 300,
  bg = "white"
)

# -------------------
# 14. Combined plot with labels A and B
# Landscape-page optimized
# -------------------
combined_plot <- (p_higher | p_lower) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 20, face = "bold"),
      plot.margin = margin(6, 6, 6, 6)
    )
  )

print(combined_plot)

# Standard landscape page: close to A4 landscape aspect
ggsave(
  file.path(fig_dir, "SuppFigS1_AB_HigherLowerAbundance_Species_log2FC_pvalues_LANDSCAPE.png"),
  plot = combined_plot,
  width = 15.5,
  height = 7.2,
  dpi = 300,
  bg = "white"
)

# Extra-wide version, useful if species names or p-values are clipped
ggsave(
  file.path(fig_dir, "SuppFigS1_AB_HigherLowerAbundance_Species_log2FC_pvalues_EXTRA_WIDE.png"),
  plot = combined_plot,
  width = 17,
  height = 7.2,
  dpi = 300,
  bg = "white"
)

cat("\nSupplementary Figure S1 complete.\n")
cat("Individual panels and combined landscape A/B figure saved to:\n", fig_dir, "\n")
cat("Shared x-axis limits used: ", paste(x_limits, collapse = " to "), "\n")

library(phyloseq)
library(ggplot2)
library(dplyr)

# Define the output directory for saving plots
output_dir <- file.path("results", "supplementary_figures", "family_comparison")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Subset data for Case and Control groups
physeq_case <- subset_samples(physeq_rel_abund, Condition == "PFAS_Exposed")
physeq_control <- subset_samples(physeq_rel_abund, Condition == "Reference")

# Extract OTU and Taxonomy tables
otu_case <- as.data.frame(otu_table(physeq_case))
otu_control <- as.data.frame(otu_table(physeq_control))
taxonomy_case <- as.data.frame(tax_table(physeq_case))
taxonomy_control <- as.data.frame(tax_table(physeq_control))

# Identify unique phyla in both datasets
unique_phyla <- unique(c(taxonomy_case$Family, taxonomy_control$Family))
unique_phyla <- unique_phyla[!is.na(unique_phyla)]  # Remove any NAs

# Loop through each Family and generate a plot
for (Family in unique_phyla) {
  # Get indices of OTUs belonging to the current Family
  case_idx <- which(taxonomy_case$Family == Family)
  control_idx <- which(taxonomy_control$Family == Family)
  
  # Extract OTU abundance for the Family from both groups
  Family_case <- otu_case[case_idx, , drop = FALSE]
  Family_control <- otu_control[control_idx, , drop = FALSE]
  
  # Sum the abundance for the Family in each sample
  abundance_case <- rowSums(Family_case)
  abundance_control <- rowSums(Family_control)
  
  # Create a data frame for plotting
  Family_data <- data.frame(
    Abundance = c(abundance_case, abundance_control),
    Group = rep(c("PFAS-Exposed", "Reference"), c(length(abundance_case), length(abundance_control)))
  )
  
  # Skip if no data for this Family
  if (nrow(Family_data) == 0) next
  
  # Perform statistical test (Wilcoxon test)
  wilcox_test <- wilcox.test(Abundance ~ Group, data = Family_data)
  
  # Generate the plot
  p <- ggplot(Family_data, aes(x = Group, y = Abundance, fill = Group)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste("Relative Abundance of", Family), 
         x = "Group", y = "Abundance") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10, face = "bold")
    ) +
    annotate("text", x = 1.5, y = max(Family_data$Abundance, na.rm = TRUE), 
             label = paste("p-value =", round(wilcox_test$p.value, 4)), 
             size = 5, hjust = 0.5)
  
  # Save the plot
  filename <- paste0(output_dir, "/", gsub(" ", "_", Family), "_abundance_comparison.png")
  ggsave(filename = filename, plot = p, width = 12, height = 8, dpi = 300)
  
  # Print message for tracking progress
  message("Saved plot for: ", Family)
}
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

sample_metadata_path <- file.path("data", "metadata", "Sample_Metadata_New.tsv")

# Load and Preprocess Data
sample_metadata <- read.delim(sample_metadata_path, check.names = FALSE, sep = "\t", row.names = 1)
names(sample_metadata) <- make.names(names(sample_metadata), unique = TRUE)

# Filter PFAS_Exposed group
pfas_exposed_data <- sample_metadata %>% 
  filter(Condition == "PFAS_Exposed")

# Select degradation rate columns
pfas_exposed_data_k <- pfas_exposed_data %>%
  as_tibble() %>%
  dplyr::select(starts_with("K_"))

# Reshape to long format
pfas_exposed_data_long <- pfas_exposed_data_k %>%
  pivot_longer(cols = everything(), 
               names_to = "PFAS_Type", 
               values_to = "Elimination_Rate_Value")

# Define dark purple-based palette
dark_palette <- c(
  "#3E1F47", "#5D3A6B", "#7A4D89", "#9A67A1",
  "#BC8CC3", "#6E2C7B", "#421C52", "#8E5BA3"
)

# Plot
plot <- ggplot(pfas_exposed_data_long, aes(x = PFAS_Type, y = Elimination_Rate_Value, fill = PFAS_Type)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "white", outlier.color = "black", width = 0.6) +
  scale_fill_manual(values = dark_palette) +
  theme_minimal(base_size = 14) +
  labs(
    x = "PFAS Type", 
    y = "Elimination Rate Value", 
    title = ""
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "#2D0033"),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12, color = "#2D0033"),
    axis.text.y = element_text(face = "bold", size = 12, color = "#2D0033"),
    axis.title = element_text(face = "bold", size = 14, color = "#2D0033"),
    legend.position = "none", # remove redundant legend
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 10, 10, 10) # tighter margins, less whitespace
  )

# Save
ggsave(
  file.path("results", "supplementary_figures", "Supp2_pfas_elimination_boxplot.png"), 
  plot = plot, 
  width = 9, height = 5, dpi = 300, bg = 'white'
)


#########################################
#########################################
# Figure 3 + Supplementary Figure
# PFAS–taxa correlations
#
# FINAL VERSION WITH HIERARCHICAL TAXON ORDERING + CLUSTER LINES
#
# MAIN FIGURE 3:
#   A. ΣPFAS: PFAS-exposed
#   B. ΣPFAS: Reference
#   C. Average PFAS elimination rate: PFAS-exposed only
#
# SUPPLEMENTARY FIGURE:
#   A. Individual serum PFAS values: PFAS-exposed
#   B. Individual serum PFAS values: Reference
#   C. Individual PFAS elimination-rate constants: PFAS-exposed only
#
# Key points:
# - Main figure uses only three simplified endpoints.
# - Supplementary figure shows compound-specific PFAS and k values.
# - Taxa are ordered by hierarchical clustering of the plotted correlation matrix.
# - Cluster boundaries are shown with horizontal separator lines.
# - Reduced horizontal white space.
# - One shared Spearman rho legend and one shared -log10(p-value) legend.
# - Black outline around circles.
# - Bold italic species labels.
# - Exports full and plotted correlation tables.
#########################################
#########################################

# -----------------------------
# 0. Load packages
# -----------------------------
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(tibble)
library(scales)
library(cowplot)

# -----------------------------
# 1. Output directory
# -----------------------------
output_dir <- file.path("results", "supplementary_figures")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# -----------------------------
# 2. User inputs
# -----------------------------

sum_pfas_columns <- c(
  "PFOS",
  "PFHxS",
  "PFOA",
  "PFNA"
)

serum_pfas_columns <- c(
  "PFOA",
  "PFHxS",
  "PFOS",
  "PFNA",
  "PFOS_branched_MP11",
  "PFOS_branchedMP3_4_51",
  "PFOS_branchedMP2_61"
)

k_columns <- c(
  "k_PFOA",
  "k_PFPeS",
  "k_PFHxS",
  "k_PFHpS",
  "k_LPFOS",
  "k_PFOS_MP1",
  "k_PFOS_MP345",
  "k_PFOS_MP26"
)

candidate_average_k_columns <- c(
  "average_k",
  "Average_k",
  "k_average",
  "K_average",
  "mean_k",
  "Mean_k",
  "k_sum",
  "K_sum",
  "sum_k",
  "Sum_k",
  "k_sumPFAS",
  "k_SumPFAS"
)

top_n_per_main_endpoint <- 18
max_taxa_main <- 24

plot_all_taxa_in_supplement <- TRUE
top_supp_taxa_total <- 100

min_n_cor <- 5
max_minus_log10_p <- 3

# Recommended: hclust
# Other options: "alphabetical", "effect_size"
taxon_ordering_method <- "hclust"

# Number of visible clusters used only for drawing separator lines
main_k_clusters <- 4
supp_k_clusters <- 6

# -----------------------------
# 3. Helper functions
# -----------------------------

abbreviate_species <- function(name) {
  if (is.na(name) || name == "") return(NA_character_)
  name <- gsub("_", " ", name)
  parts <- strsplit(name, " ")[[1]]
  if (length(parts) == 1) return(parts[1])
  paste0(substr(parts[1], 1, 1), ". ", paste(parts[-1], collapse = " "))
}

to_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "nd", "ND", "Na", "na", "n/a", "N/A",
             "<LOD", "< LOQ", "<LOQ", "LOD", "LOQ")] <- NA
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- gsub("^<\\s*", "", x)
  suppressWarnings(as.numeric(x))
}

safe_spearman <- function(data, x_col, y_col, min_n = 5) {
  
  dat <- data %>%
    dplyr::select(all_of(c(x_col, y_col))) %>%
    mutate(
      x_tmp = to_numeric_safe(.data[[x_col]]),
      y_tmp = to_numeric_safe(.data[[y_col]])
    ) %>%
    filter(
      !is.na(x_tmp),
      !is.na(y_tmp),
      is.finite(x_tmp),
      is.finite(y_tmp)
    )
  
  if (nrow(dat) < min_n) {
    return(tibble(n = nrow(dat), Correlation = NA_real_, P.Value = NA_real_))
  }
  
  if (length(unique(dat$x_tmp)) < 2 || length(unique(dat$y_tmp)) < 2) {
    return(tibble(n = nrow(dat), Correlation = NA_real_, P.Value = NA_real_))
  }
  
  test <- suppressWarnings(
    cor.test(dat$x_tmp, dat$y_tmp, method = "spearman", exact = FALSE)
  )
  
  tibble(
    n = nrow(dat),
    Correlation = unname(test$estimate),
    P.Value = test$p.value
  )
}

pretty_endpoint <- function(x) {
  out <- x
  
  out[out == "sum_PFAS_exposed"] <- "\u03A3PFAS\nPFAS-exposed"
  out[out == "sum_PFAS_reference"] <- "\u03A3PFAS\nReference"
  out[out == "average_k_main"] <- "Average PFAS\nelimination rate"
  
  out[out == "PFOA"] <- "PFOA"
  out[out == "PFHxS"] <- "PFHxS"
  out[out == "PFOS"] <- "PFOS"
  out[out == "PFNA"] <- "PFNA"
  out[out == "PFOS_branched_MP11"] <- "PFOS MP1"
  out[out == "PFOS_branchedMP3_4_51"] <- "PFOS MP3/4/5"
  out[out == "PFOS_branchedMP2_61"] <- "PFOS MP2/6"
  
  out[out == "k_PFOA"] <- "k PFOA"
  out[out == "k_PFPeS"] <- "k PFPeS"
  out[out == "k_PFHxS"] <- "k PFHxS"
  out[out == "k_PFHpS"] <- "k PFHpS"
  out[out == "k_LPFOS"] <- "k L-PFOS"
  out[out == "k_PFOS_MP1"] <- "k PFOS MP1"
  out[out == "k_PFOS_MP345"] <- "k PFOS MP3/4/5"
  out[out == "k_PFOS_MP26"] <- "k PFOS MP2/6"
  
  out
}

# Hierarchical ordering helper.
# Returns:
#   order         = taxon order
#   cluster_lines = horizontal boundary positions between clusters
#   clusters      = cluster assignment in plotted order
order_taxa_for_plot <- function(df, method = "hclust", k_clusters = 4) {
  
  required_cols <- c("Taxon_Readable", "Endpoint_Label", "Correlation")
  if (!all(required_cols %in% colnames(df))) {
    stop("order_taxa_for_plot requires Taxon_Readable, Endpoint_Label, and Correlation columns.")
  }
  
  df_clean <- df %>%
    filter(
      !is.na(Taxon_Readable),
      !is.na(Endpoint_Label),
      !is.na(Correlation),
      is.finite(Correlation)
    ) %>%
    mutate(
      Taxon_Readable = as.character(Taxon_Readable),
      Endpoint_Label = as.character(Endpoint_Label)
    )
  
  if (nrow(df_clean) == 0) {
    stop("No valid rows available for taxon ordering.")
  }
  
  if (method == "alphabetical") {
    taxon_order <- sort(unique(df_clean$Taxon_Readable), decreasing = TRUE)
    return(list(
      order = taxon_order,
      cluster_lines = numeric(0),
      clusters = NULL
    ))
  }
  
  if (method == "effect_size") {
    taxon_order <- df_clean %>%
      group_by(Taxon_Readable) %>%
      summarise(
        max_abs_rho = max(abs(Correlation), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(max_abs_rho) %>%
      pull(Taxon_Readable)
    
    return(list(
      order = taxon_order,
      cluster_lines = numeric(0),
      clusters = NULL
    ))
  }
  
  if (method == "hclust") {
    
    # Collapse duplicates before pivoting.
    wide_df <- df_clean %>%
      group_by(Taxon_Readable, Endpoint_Label) %>%
      summarise(
        Correlation = mean(Correlation, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      pivot_wider(
        names_from = Endpoint_Label,
        values_from = Correlation,
        values_fill = 0
      ) %>%
      as.data.frame()
    
    rownames(wide_df) <- wide_df$Taxon_Readable
    wide_df$Taxon_Readable <- NULL
    
    wide_mat <- as.matrix(wide_df)
    storage.mode(wide_mat) <- "numeric"
    wide_mat[is.na(wide_mat)] <- 0
    
    if (nrow(wide_mat) < 2 || ncol(wide_mat) < 2) {
      taxon_order <- df_clean %>%
        group_by(Taxon_Readable) %>%
        summarise(
          max_abs_rho = max(abs(Correlation), na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(max_abs_rho) %>%
        pull(Taxon_Readable)
      
      return(list(
        order = taxon_order,
        cluster_lines = numeric(0),
        clusters = NULL
      ))
    }
    
    hc <- hclust(
      dist(wide_mat, method = "euclidean"),
      method = "ward.D2"
    )
    
    taxon_order <- rownames(wide_mat)[hc$order]
    
    k_use <- min(k_clusters, length(taxon_order))
    clusters <- cutree(hc, k = k_use)
    clusters_ordered <- clusters[taxon_order]
    
    cluster_lines <- which(
      clusters_ordered[-1] != clusters_ordered[-length(clusters_ordered)]
    ) + 0.5
    
    return(list(
      order = taxon_order,
      cluster_lines = cluster_lines,
      clusters = clusters_ordered
    ))
  }
  
  stop("Unknown taxon ordering method. Use 'hclust', 'alphabetical', or 'effect_size'.")
}

# -----------------------------
# 4. Extract abundance and metadata
# -----------------------------

taxa_abundances <- as.data.frame(otu_table(vst_physeq))

if (taxa_are_rows(vst_physeq)) {
  taxa_abundances <- as.data.frame(t(taxa_abundances))
} else {
  taxa_abundances <- as.data.frame(taxa_abundances)
}

metadata <- data.frame(sample_data(vst_physeq), check.names = FALSE)
metadata$SampleID <- rownames(metadata)

common_samples <- intersect(metadata$SampleID, rownames(taxa_abundances))

metadata <- metadata %>%
  filter(SampleID %in% common_samples)

taxa_abundances <- taxa_abundances[metadata$SampleID, , drop = FALSE]

metadata <- metadata %>%
  filter(
    !is.na(Condition),
    Condition %in% c("PFAS_Exposed", "Reference")
  ) %>%
  mutate(
    Condition = factor(Condition, levels = c("PFAS_Exposed", "Reference"))
  )

taxa_abundances <- taxa_abundances[metadata$SampleID, , drop = FALSE]

taxa_cols <- colnames(taxa_abundances)

keep_samples <- !apply(
  taxa_abundances[, taxa_cols, drop = FALSE],
  1,
  function(x) all(is.na(x)) || all(x == 0, na.rm = TRUE)
)

metadata <- metadata[keep_samples, , drop = FALSE]
taxa_abundances <- taxa_abundances[keep_samples, , drop = FALSE]

merged_data <- bind_cols(metadata, as.data.frame(taxa_abundances))

cat("\nSamples included by condition:\n")
print(table(merged_data$Condition, useNA = "ifany"))

# -----------------------------
# 5. Convert endpoint columns to numeric
# -----------------------------

available_sum_pfas_cols <- intersect(sum_pfas_columns, colnames(merged_data))

if (length(available_sum_pfas_cols) == 0) {
  stop("No PFAS columns found for ΣPFAS. Check column names in metadata.")
}

if (length(available_sum_pfas_cols) < 3) {
  warning(
    "Fewer than three PFAS columns available for ΣPFAS: ",
    paste(available_sum_pfas_cols, collapse = ", ")
  )
}

available_serum_pfas_cols <- intersect(serum_pfas_columns, colnames(merged_data))

if (length(available_serum_pfas_cols) == 0) {
  warning("No individual serum PFAS columns found for supplementary serum panels.")
}

available_k_cols <- intersect(k_columns, colnames(merged_data))

if (length(available_k_cols) == 0) {
  stop("No k columns found. Check elimination-rate column names in metadata.")
}

available_average_k_col <- intersect(candidate_average_k_columns, colnames(merged_data))

endpoint_cols_to_convert <- unique(c(
  available_sum_pfas_cols,
  available_serum_pfas_cols,
  available_k_cols,
  available_average_k_col
))

for (col in endpoint_cols_to_convert) {
  merged_data[[col]] <- to_numeric_safe(merged_data[[col]])
}

# -----------------------------
# 6. Build ΣPFAS and average PFAS elimination rate
# -----------------------------

merged_data <- merged_data %>%
  mutate(
    n_sum_pfas_available = rowSums(!is.na(across(all_of(available_sum_pfas_cols)))),
    sum_PFAS = rowSums(across(all_of(available_sum_pfas_cols)), na.rm = TRUE),
    sum_PFAS = ifelse(n_sum_pfas_available == 0, NA, sum_PFAS)
  )

cat("\nΣPFAS calculated using:\n")
print(available_sum_pfas_cols)

if (length(available_average_k_col) > 0) {
  
  average_k_col <- available_average_k_col[1]
  merged_data$average_k_main <- to_numeric_safe(merged_data[[average_k_col]])
  
  cat("\nUsing existing average/composite k column:\n")
  print(average_k_col)
  
} else {
  
  merged_data <- merged_data %>%
    mutate(
      n_k_available = rowSums(!is.na(across(all_of(available_k_cols)))),
      average_k_main = rowMeans(across(all_of(available_k_cols)), na.rm = TRUE),
      average_k_main = ifelse(n_k_available == 0, NA, average_k_main)
    )
  
  cat("\nCalculated average_k_main from:\n")
  print(available_k_cols)
}

# -----------------------------
# 7. Taxonomy mapping
# -----------------------------

tax_table_df <- as.data.frame(tax_table(vst_physeq))
tax_table_df$Taxon <- rownames(tax_table_df)

if (!"Species" %in% colnames(tax_table_df)) {
  tax_table_df$Species <- tax_table_df$Taxon
}

tax_name_map <- tax_table_df %>%
  mutate(
    Taxon_Readable = ifelse(is.na(Species) | Species == "", Taxon, Species),
    Taxon_Readable = sapply(Taxon_Readable, abbreviate_species)
  ) %>%
  dplyr::select(Taxon, Taxon_Readable)

# -----------------------------
# 8. Run MAIN correlations
# -----------------------------

main_results <- list()

# Main A: ΣPFAS, PFAS-exposed
for (taxon in taxa_cols) {
  
  dat_i <- merged_data %>%
    filter(Condition == "PFAS_Exposed")
  
  res_i <- safe_spearman(
    data = dat_i,
    x_col = "sum_PFAS",
    y_col = taxon,
    min_n = min_n_cor
  )
  
  main_results[[length(main_results) + 1]] <- res_i %>%
    mutate(
      Taxon = taxon,
      Endpoint = "sum_PFAS_exposed",
      Panel = "A",
      Panel_Title = "\u03A3PFAS: PFAS-exposed",
      Condition = "PFAS_Exposed",
      Figure_Set = "Main"
    )
}

# Main B: ΣPFAS, Reference
for (taxon in taxa_cols) {
  
  dat_i <- merged_data %>%
    filter(Condition == "Reference")
  
  res_i <- safe_spearman(
    data = dat_i,
    x_col = "sum_PFAS",
    y_col = taxon,
    min_n = min_n_cor
  )
  
  main_results[[length(main_results) + 1]] <- res_i %>%
    mutate(
      Taxon = taxon,
      Endpoint = "sum_PFAS_reference",
      Panel = "B",
      Panel_Title = "\u03A3PFAS: Reference",
      Condition = "Reference",
      Figure_Set = "Main"
    )
}

# Main C: average k, PFAS-exposed only
for (taxon in taxa_cols) {
  
  dat_i <- merged_data %>%
    filter(Condition == "PFAS_Exposed")
  
  res_i <- safe_spearman(
    data = dat_i,
    x_col = "average_k_main",
    y_col = taxon,
    min_n = min_n_cor
  )
  
  main_results[[length(main_results) + 1]] <- res_i %>%
    mutate(
      Taxon = taxon,
      Endpoint = "average_k_main",
      Panel = "C",
      Panel_Title = "Average PFAS elimination rate:\nPFAS-exposed only",
      Condition = "PFAS_Exposed",
      Figure_Set = "Main"
    )
}

main_results <- bind_rows(main_results) %>%
  left_join(tax_name_map, by = "Taxon") %>%
  group_by(Panel) %>%
  mutate(
    Adj.P.Value = p.adjust(P.Value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    Endpoint_Label = pretty_endpoint(Endpoint),
    minus_log10_p_raw = -log10(P.Value),
    minus_log10_p = pmin(minus_log10_p_raw, max_minus_log10_p)
  )

write.csv(
  main_results,
  file.path(output_dir, "Figure3_MAIN_All_TaxaCorrelation_Data.csv"),
  row.names = FALSE
)

# -----------------------------
# 9. Select taxa for MAIN figure
# -----------------------------

main_taxa_to_plot <- main_results %>%
  filter(!is.na(Correlation), !is.na(P.Value)) %>%
  group_by(Panel) %>%
  slice_max(order_by = abs(Correlation), n = top_n_per_main_endpoint, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(Taxon_Readable) %>%
  summarise(
    max_abs_rho = max(abs(Correlation), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_rho)) %>%
  slice_head(n = max_taxa_main) %>%
  pull(Taxon_Readable)

main_plot_df <- main_results %>%
  filter(Taxon_Readable %in% main_taxa_to_plot) %>%
  filter(!is.na(Correlation), !is.na(P.Value))

main_taxon_info <- order_taxa_for_plot(
  main_plot_df,
  method = taxon_ordering_method,
  k_clusters = main_k_clusters
)

main_taxon_order <- main_taxon_info$order
main_cluster_lines <- main_taxon_info$cluster_lines

main_plot_df <- main_plot_df %>%
  mutate(
    Taxon_Readable = factor(Taxon_Readable, levels = main_taxon_order)
  )

endpoint_order_main_A <- pretty_endpoint("sum_PFAS_exposed")
endpoint_order_main_B <- pretty_endpoint("sum_PFAS_reference")
endpoint_order_main_C <- pretty_endpoint("average_k_main")

# -----------------------------
# 10. MAIN plotting scales
# -----------------------------

main_rho_limit <- max(abs(main_plot_df$Correlation), na.rm = TRUE)
main_rho_limit <- max(main_rho_limit, 0.1)

main_size_breaks <- c(0.5, 1.0, 1.5, 2.0)
main_size_breaks <- main_size_breaks[
  main_size_breaks <= max(main_plot_df$minus_log10_p, na.rm = TRUE) + 0.05
]
if (length(main_size_breaks) == 0) main_size_breaks <- c(0.5)

# -----------------------------
# 11. Plotting function for main panels
# -----------------------------

make_panel <- function(
    df,
    panel_title,
    endpoint_levels,
    rho_limit,
    size_breaks,
    show_y = TRUE,
    cluster_lines = NULL
) {
  
  df <- df %>%
    mutate(
      Endpoint_Label = factor(Endpoint_Label, levels = endpoint_levels)
    )
  
  cluster_df <- data.frame(yintercept = cluster_lines)
  
  p <- ggplot(
    df,
    aes(
      x = Endpoint_Label,
      y = Taxon_Readable,
      fill = Correlation,
      size = minus_log10_p
    )
  ) +
    {
      if (!is.null(cluster_lines) && length(cluster_lines) > 0) {
        geom_hline(
          data = cluster_df,
          aes(yintercept = yintercept),
          inherit.aes = FALSE,
          linewidth = 0.45,
          color = "grey45",
          linetype = "solid"
        )
      }
    } +
    geom_point(
      shape = 21,
      color = "black",
      stroke = 0.35,
      alpha = 0.96
    ) +
    scale_fill_gradient2(
      low = "#B2182B",
      mid = "white",
      high = "#6A1B9A",
      midpoint = 0,
      limits = c(-rho_limit, rho_limit),
      oob = squish,
      name = "Spearman\nrho"
    ) +
    scale_size_continuous(
      range = c(1.5, 5.4),
      breaks = size_breaks,
      name = "-log10\np-value"
    ) +
    scale_y_discrete(
      drop = FALSE,
      expand = expansion(add = 0.38)
    ) +
    scale_x_discrete(
      drop = FALSE,
      expand = expansion(add = 0.22)
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = if (show_y) "Species" else NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 9.5,
        face = "bold",
        color = "grey20"
      ),
      axis.text.y = element_text(
        size = 9.4,
        face = "bold.italic",
        color = "grey20"
      ),
      axis.title.y = element_text(size = 13, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "grey84", linewidth = 0.28),
      panel.grid.major.y = element_line(color = "grey94", linewidth = 0.20),
      panel.border = element_rect(color = "grey45", fill = NA, linewidth = 0.45),
      legend.position = "none",
      plot.margin = margin(4, 3, 4, 3)
    )
  
  if (!show_y) {
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()
      )
  }
  
  p
}

# -----------------------------
# 12. Build MAIN panels
# -----------------------------

p_main_A <- make_panel(
  main_plot_df %>% filter(Panel == "A"),
  "\u03A3PFAS: PFAS-exposed",
  endpoint_order_main_A,
  main_rho_limit,
  main_size_breaks,
  show_y = TRUE,
  cluster_lines = main_cluster_lines
)

p_main_B <- make_panel(
  main_plot_df %>% filter(Panel == "B"),
  "\u03A3PFAS: Reference",
  endpoint_order_main_B,
  main_rho_limit,
  main_size_breaks,
  show_y = FALSE,
  cluster_lines = main_cluster_lines
)

p_main_C <- make_panel(
  main_plot_df %>% filter(Panel == "C"),
  "Average PFAS elimination rate:\nPFAS-exposed only",
  endpoint_order_main_C,
  main_rho_limit,
  main_size_breaks,
  show_y = FALSE,
  cluster_lines = main_cluster_lines
)

# -----------------------------
# 13. Shared legend for MAIN
# -----------------------------

legend_main <- ggplot(
  main_plot_df,
  aes(
    x = Endpoint_Label,
    y = Taxon_Readable,
    fill = Correlation,
    size = minus_log10_p
  )
) +
  geom_point(shape = 21, color = "black", stroke = 0.35, alpha = 0.96) +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#6A1B9A",
    midpoint = 0,
    limits = c(-main_rho_limit, main_rho_limit),
    oob = squish,
    name = "Spearman\nrho"
  ) +
  scale_size_continuous(
    range = c(1.5, 5.4),
    breaks = main_size_breaks,
    name = "-log10\np-value"
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barwidth = 5.5,
      barheight = 0.55
    ),
    size = guide_legend(
      title.position = "top",
      override.aes = list(fill = "grey70", color = "black")
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

shared_legend_main <- cowplot::get_legend(legend_main)

main_panel_row <- (p_main_A | p_main_B | p_main_C) +
  plot_layout(
    ncol = 3,
    widths = c(0.72, 0.72, 0.82)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 20, face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )
  )

main_combined <- cowplot::plot_grid(
  main_panel_row,
  shared_legend_main,
  ncol = 1,
  rel_heights = c(1, 0.13)
)

print(main_combined)

# -----------------------------
# 14. Save MAIN
# -----------------------------

n_taxa_main <- length(levels(main_plot_df$Taxon_Readable))

main_height <- max(7.0, min(10.2, 2.4 + 0.22 * n_taxa_main))
main_width <- 8.8

ggsave(
  file.path(output_dir, "Figure3_MAIN_SumPFAS_AverageK_TaxaCorrelations_ABC.png"),
  plot = main_combined,
  width = main_width,
  height = main_height,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Figure3_MAIN_SumPFAS_AverageK_TaxaCorrelations_ABC.pdf"),
  plot = main_combined,
  width = main_width,
  height = main_height,
  bg = "white"
)

write.csv(
  main_plot_df,
  file.path(output_dir, "Figure3_MAIN_PLOTTED_SumPFAS_AverageK_TaxaCorrelation_Data.csv"),
  row.names = FALSE
)

# ==========================================================
# SUPPLEMENTARY FIGURE: individual PFAS and individual k values
# ==========================================================

supp_results <- list()

# Supplementary A: individual serum PFAS, PFAS-exposed
for (taxon in taxa_cols) {
  for (endpoint in available_serum_pfas_cols) {
    
    dat_i <- merged_data %>%
      filter(Condition == "PFAS_Exposed")
    
    res_i <- safe_spearman(
      data = dat_i,
      x_col = endpoint,
      y_col = taxon,
      min_n = min_n_cor
    )
    
    supp_results[[length(supp_results) + 1]] <- res_i %>%
      mutate(
        Taxon = taxon,
        Endpoint = endpoint,
        Panel = "A",
        Panel_Title = "Individual serum PFAS:\nPFAS-exposed",
        Condition = "PFAS_Exposed",
        Figure_Set = "Supplementary"
      )
  }
}

# Supplementary B: individual serum PFAS, Reference
for (taxon in taxa_cols) {
  for (endpoint in available_serum_pfas_cols) {
    
    dat_i <- merged_data %>%
      filter(Condition == "Reference")
    
    res_i <- safe_spearman(
      data = dat_i,
      x_col = endpoint,
      y_col = taxon,
      min_n = min_n_cor
    )
    
    supp_results[[length(supp_results) + 1]] <- res_i %>%
      mutate(
        Taxon = taxon,
        Endpoint = endpoint,
        Panel = "B",
        Panel_Title = "Individual serum PFAS:\nReference",
        Condition = "Reference",
        Figure_Set = "Supplementary"
      )
  }
}

# Supplementary C: individual k, PFAS-exposed only
for (taxon in taxa_cols) {
  for (endpoint in available_k_cols) {
    
    dat_i <- merged_data %>%
      filter(Condition == "PFAS_Exposed")
    
    res_i <- safe_spearman(
      data = dat_i,
      x_col = endpoint,
      y_col = taxon,
      min_n = min_n_cor
    )
    
    supp_results[[length(supp_results) + 1]] <- res_i %>%
      mutate(
        Taxon = taxon,
        Endpoint = endpoint,
        Panel = "C",
        Panel_Title = "Individual elimination rates:\nPFAS-exposed only",
        Condition = "PFAS_Exposed",
        Figure_Set = "Supplementary"
      )
  }
}

supp_results <- bind_rows(supp_results) %>%
  left_join(tax_name_map, by = "Taxon") %>%
  group_by(Panel) %>%
  mutate(
    Adj.P.Value = p.adjust(P.Value, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    Endpoint_Label = pretty_endpoint(Endpoint),
    minus_log10_p_raw = -log10(P.Value),
    minus_log10_p = pmin(minus_log10_p_raw, max_minus_log10_p)
  )

write.csv(
  supp_results,
  file.path(output_dir, "Figure3_SUPPLEMENTARY_All_CompoundSpecific_TaxaCorrelation_Data.csv"),
  row.names = FALSE
)

# -----------------------------
# 15. Select data for supplementary plotting
# -----------------------------

supp_plot_df <- supp_results %>%
  filter(!is.na(Correlation), !is.na(P.Value))

if (!plot_all_taxa_in_supplement) {
  
  supp_taxa_to_plot <- supp_plot_df %>%
    group_by(Taxon_Readable) %>%
    summarise(
      max_abs_rho = max(abs(Correlation), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(max_abs_rho)) %>%
    slice_head(n = top_supp_taxa_total) %>%
    pull(Taxon_Readable)
  
  supp_plot_df <- supp_plot_df %>%
    filter(Taxon_Readable %in% supp_taxa_to_plot)
}

supp_taxon_info <- order_taxa_for_plot(
  supp_plot_df,
  method = taxon_ordering_method,
  k_clusters = supp_k_clusters
)

supp_taxon_order <- supp_taxon_info$order
supp_cluster_lines <- supp_taxon_info$cluster_lines

supp_plot_df <- supp_plot_df %>%
  mutate(
    Taxon_Readable = factor(Taxon_Readable, levels = supp_taxon_order),
    Panel_Title = factor(
      Panel_Title,
      levels = c(
        "Individual serum PFAS:\nPFAS-exposed",
        "Individual serum PFAS:\nReference",
        "Individual elimination rates:\nPFAS-exposed only"
      )
    )
  )

# -----------------------------
# 16. Supplementary plotting scales
# -----------------------------

supp_rho_limit <- max(abs(supp_plot_df$Correlation), na.rm = TRUE)
supp_rho_limit <- max(supp_rho_limit, 0.1)

supp_size_breaks <- c(0.5, 1.0, 1.5, 2.0)
supp_size_breaks <- supp_size_breaks[
  supp_size_breaks <= max(supp_plot_df$minus_log10_p, na.rm = TRUE) + 0.05
]
if (length(supp_size_breaks) == 0) supp_size_breaks <- c(0.5)

# -----------------------------
# 17. Supplementary figure
# -----------------------------

p_supp <- ggplot(
  supp_plot_df,
  aes(
    x = Endpoint_Label,
    y = Taxon_Readable,
    fill = Correlation,
    size = minus_log10_p
  )
) +
  geom_hline(
    data = data.frame(yintercept = supp_cluster_lines),
    aes(yintercept = yintercept),
    inherit.aes = FALSE,
    linewidth = 0.35,
    color = "grey45",
    linetype = "solid"
  ) +
  geom_point(
    shape = 21,
    color = "black",
    stroke = 0.24,
    alpha = 0.95
  ) +
  facet_grid(
    . ~ Panel_Title,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "white",
    high = "#6A1B9A",
    midpoint = 0,
    limits = c(-supp_rho_limit, supp_rho_limit),
    oob = squish,
    name = "Spearman\nrho"
  ) +
  scale_size_continuous(
    range = c(0.85, 4.5),
    breaks = supp_size_breaks,
    name = "-log10\np-value"
  ) +
  scale_y_discrete(
    drop = FALSE,
    expand = expansion(add = 0.38)
  ) +
  scale_x_discrete(
    drop = FALSE,
    expand = expansion(add = 0.18)
  ) +
  labs(
    x = NULL,
    y = "Species"
  ) +
  coord_cartesian(clip = "off") +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barwidth = 6,
      barheight = 0.55
    ),
    size = guide_legend(
      title.position = "top",
      override.aes = list(fill = "grey70", color = "black")
    )
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12.5, face = "bold", hjust = 0.5),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 8.4,
      face = "bold",
      color = "grey20"
    ),
    axis.text.y = element_text(
      size = 7.2,
      face = "bold.italic",
      color = "grey20"
    ),
    axis.title.y = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey84", linewidth = 0.25),
    panel.grid.major.y = element_line(color = "grey93", linewidth = 0.18),
    panel.border = element_rect(color = "grey45", fill = NA, linewidth = 0.45),
    panel.spacing.x = unit(0.35, "lines"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.title = element_text(size = 10.5, face = "bold"),
    legend.text = element_text(size = 9.5),
    plot.margin = margin(6, 6, 6, 6)
  )

print(p_supp)

n_taxa_supp <- length(levels(supp_plot_df$Taxon_Readable))

supp_height <- if (plot_all_taxa_in_supplement) {
  max(10, 2.8 + 0.115 * n_taxa_supp)
} else {
  max(9, 2.8 + 0.15 * n_taxa_supp)
}

supp_width <- 12.8

ggsave(
  file.path(output_dir, "Figure3_SUPPLEMENTARY_All_CompoundSpecific_TaxaCorrelations.png"),
  plot = p_supp,
  width = supp_width,
  height = supp_height,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(output_dir, "Figure3_SUPPLEMENTARY_All_CompoundSpecific_TaxaCorrelations.pdf"),
  plot = p_supp,
  width = supp_width,
  height = supp_height,
  bg = "white",
  limitsize = FALSE
)

write.csv(
  supp_plot_df,
  file.path(output_dir, "Figure3_SUPPLEMENTARY_PLOTTED_CompoundSpecific_TaxaCorrelation_Data.csv"),
  row.names = FALSE
)

# -----------------------------
# 18. Final console summary
# -----------------------------

cat("\nDone.\n")

cat("\nTaxon ordering method used:\n")
print(taxon_ordering_method)

cat("\nMain cluster separator lines:\n")
print(main_cluster_lines)

cat("\nSupplementary cluster separator lines:\n")
print(supp_cluster_lines)

cat("\nMain figure saved as:\n")
cat(file.path(output_dir, "Figure3_MAIN_SumPFAS_AverageK_TaxaCorrelations_ABC.png"), "\n")

cat("\nSupplementary figure saved as:\n")
cat(file.path(output_dir, "Figure3_SUPPLEMENTARY_All_CompoundSpecific_TaxaCorrelations.png"), "\n")

cat("\nΣPFAS calculated using:\n")
print(available_sum_pfas_cols)

cat("\nIndividual serum PFAS shown in supplementary figure:\n")
print(available_serum_pfas_cols)

cat("\nIndividual k values shown in supplementary figure:\n")
print(available_k_cols)

cat("\nNumber of taxa plotted in main figure:\n")
print(n_taxa_main)

cat("\nNumber of taxa plotted in supplementary figure:\n")
print(n_taxa_supp)

cat("\nFull main correlation data saved as:\n")
cat(file.path(output_dir, "Figure3_MAIN_All_TaxaCorrelation_Data.csv"), "\n")

cat("\nFull supplementary correlation data saved as:\n")
cat(file.path(output_dir, "Figure3_SUPPLEMENTARY_All_CompoundSpecific_TaxaCorrelation_Data.csv"), "\n")

############################################################
# Supplementary Figure S4 species relative-abundance dot plot
#
# Improved version:
# - Samples grouped by condition
# - Dot size = relative abundance (%)
# - Dot fill = log10(relative abundance %)
# - Avoids confusing log10 legend where the maximum is 0
# - Species ordered by prevalence and mean abundance
# - Samples ordered by condition and total abundance of plotted species
# - Produces readable top-species figure and optional full figure
############################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(readr)
  library(stringr)
  library(forcats)
  library(scales)
  library(grid)
})

# ============================================================
# 1. File paths
# ============================================================

data_file_path <- "C:/Users/jonat/OneDrive - University of Glasgow/PFAS_Microbiome/Manuscript/Inputs/NormCoverage_Fig3.csv"
condition_file_path <- file.path("data", "metadata", "Sample_Metadata_New.tsv")

output_dir <- file.path("results", "supplementary_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

out_top_png  <- file.path(output_dir, "Supp_Species_Relative_Abundance_DotPlot_TopSpecies.png")
out_top_pdf  <- file.path(output_dir, "Supp_Species_Relative_Abundance_DotPlot_TopSpecies.pdf")
out_full_png <- file.path(output_dir, "Supp_Species_Relative_Abundance_DotPlot_AllSpecies.png")
out_full_pdf <- file.path(output_dir, "Supp_Species_Relative_Abundance_DotPlot_AllSpecies.pdf")
out_csv      <- file.path(output_dir, "Supp_Species_Relative_Abundance_DotPlot_Data.csv")
out_summary  <- file.path(output_dir, "Supp_Species_Relative_Abundance_DotPlot_SpeciesSummary.csv")

# ============================================================
# 2. Settings
# ============================================================

top_n_species <- 55
min_abundance <- 0
save_full_plot <- TRUE

# If TRUE, exposed and reference samples are shown in separate facets.
# This usually improves readability.
facet_by_condition <- TRUE

condition_levels <- c("PFAS_Exposed", "Reference")

condition_labels <- c(
  "PFAS_Exposed" = "PFAS-exposed",
  "Reference" = "Reference"
)

# Dot-size breaks in relative abundance percentage.
# These are easier to read than raw proportions.
size_breaks_percent <- c(0.1, 1, 5, 10, 25, 50)

# ============================================================
# 3. Helper functions
# ============================================================

clean_numeric <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "nd", "ND", "n/a", "N/A", "<LOD", "<LOQ")] <- NA
  x <- gsub(",", ".", x, fixed = TRUE)
  suppressWarnings(as.numeric(x))
}

abbreviate_species <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  x <- str_squish(x)
  
  vapply(
    x,
    function(z) {
      if (is.na(z) || z == "") return(NA_character_)
      parts <- strsplit(z, " ")[[1]]
      if (length(parts) < 2) return(z)
      paste0(substr(parts[1], 1, 1), ". ", paste(parts[-1], collapse = " "))
    },
    character(1)
  )
}

percent_label_clean <- function(x) {
  paste0(format(x, trim = TRUE, scientific = FALSE), "%")
}

# ============================================================
# 4. Load data
# ============================================================

abund_data <- read.csv(
  data_file_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

condition_data <- read.delim(
  condition_file_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"Species" %in% colnames(abund_data)) {
  stop("Could not find a column named 'Species' in the abundance table.")
}

if (!"Sample_Name" %in% colnames(condition_data)) {
  stop("Could not find a column named 'Sample_Name' in the metadata table.")
}

if (!"Condition" %in% colnames(condition_data)) {
  stop("Could not find a column named 'Condition' in the metadata table.")
}

# ============================================================
# 5. Reshape and clean
# ============================================================

data_long_raw <- abund_data %>%
  pivot_longer(
    cols = -Species,
    names_to = "Sample_Name",
    values_to = "Relative_Abundance"
  ) %>%
  mutate(
    Species = as.character(Species),
    Species_Label = abbreviate_species(Species),
    Sample_Name = as.character(Sample_Name),
    Relative_Abundance = clean_numeric(Relative_Abundance)
  ) %>%
  left_join(
    condition_data %>%
      select(Sample_Name, Condition) %>%
      distinct(),
    by = "Sample_Name"
  ) %>%
  filter(
    !is.na(Condition),
    Condition %in% condition_levels,
    !is.na(Relative_Abundance),
    is.finite(Relative_Abundance),
    Relative_Abundance >= 0
  ) %>%
  mutate(
    Condition = factor(Condition, levels = condition_levels),
    Condition_Label = recode(as.character(Condition), !!!condition_labels),
    Condition_Label = factor(Condition_Label, levels = unname(condition_labels)),
    Relative_Abundance_Percent = Relative_Abundance * 100
  )

# Apply abundance filter after creating percent values.
data_long <- data_long_raw %>%
  filter(Relative_Abundance > min_abundance)

if (nrow(data_long) == 0) {
  stop("No abundance values remained after filtering.")
}

# Use a data-driven pseudocount for log10 percent abundance.
# This avoids using a fixed 1e-6 that may not match the data scale.
min_nonzero_percent <- min(
  data_long$Relative_Abundance_Percent[data_long$Relative_Abundance_Percent > 0],
  na.rm = TRUE
)

pseudocount_percent <- min_nonzero_percent / 2

data_long <- data_long %>%
  mutate(
    log10_abundance_percent = log10(Relative_Abundance_Percent + pseudocount_percent)
  )

cat("\nRelative abundance range, raw proportion:\n")
print(range(data_long$Relative_Abundance, na.rm = TRUE))

cat("\nRelative abundance range, percent:\n")
print(range(data_long$Relative_Abundance_Percent, na.rm = TRUE))

cat("\nlog10(relative abundance %) range:\n")
print(range(data_long$log10_abundance_percent, na.rm = TRUE))

cat("\nPseudocount used for log10 percent abundance:\n")
print(pseudocount_percent)

# ============================================================
# 6. Summarize species and select top species
# ============================================================

species_summary <- data_long %>%
  group_by(Species, Species_Label) %>%
  summarise(
    prevalence_n = n_distinct(Sample_Name[Relative_Abundance > 0]),
    prevalence_percent = 100 * prevalence_n / n_distinct(data_long$Sample_Name),
    mean_abundance = mean(Relative_Abundance, na.rm = TRUE),
    mean_abundance_percent = mean(Relative_Abundance_Percent, na.rm = TRUE),
    max_abundance = max(Relative_Abundance, na.rm = TRUE),
    max_abundance_percent = max(Relative_Abundance_Percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    desc(prevalence_n),
    desc(mean_abundance),
    Species
  )

write_csv(species_summary, out_summary)

top_species <- species_summary %>%
  slice_head(n = top_n_species) %>%
  pull(Species)

species_order_top <- species_summary %>%
  filter(Species %in% top_species) %>%
  arrange(prevalence_n, mean_abundance, Species) %>%
  pull(Species_Label)

species_order_full <- species_summary %>%
  arrange(prevalence_n, mean_abundance, Species) %>%
  pull(Species_Label)

# ============================================================
# 7. Order samples
# ============================================================

# For the top-species plot, order samples by condition and total abundance
# of the plotted species. This usually creates a cleaner visual structure.
sample_order_top <- data_long %>%
  filter(Species %in% top_species) %>%
  group_by(Sample_Name, Condition, Condition_Label) %>%
  summarise(
    total_top_species_percent = sum(Relative_Abundance_Percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    Condition,
    desc(total_top_species_percent),
    Sample_Name
  ) %>%
  pull(Sample_Name)

sample_order_full <- data_long %>%
  group_by(Sample_Name, Condition, Condition_Label) %>%
  summarise(
    total_abundance_percent = sum(Relative_Abundance_Percent, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(
    Condition,
    desc(total_abundance_percent),
    Sample_Name
  ) %>%
  pull(Sample_Name)

plot_data_top <- data_long %>%
  filter(Species %in% top_species) %>%
  mutate(
    Sample_Name = factor(Sample_Name, levels = sample_order_top),
    Species_Label = factor(Species_Label, levels = species_order_top)
  )

plot_data_full <- data_long %>%
  mutate(
    Sample_Name = factor(Sample_Name, levels = sample_order_full),
    Species_Label = factor(Species_Label, levels = species_order_full)
  )

write_csv(data_long, out_csv)

# ============================================================
# 8. Plot function
# ============================================================

make_species_dotplot <- function(
    plot_df,
    title_text = NULL,
    y_text_size = 7.5,
    facet_by_condition = TRUE
) {
  
  max_percent <- max(plot_df$Relative_Abundance_Percent, na.rm = TRUE)
  
  size_breaks_use <- size_breaks_percent[
    size_breaks_percent <= max_percent
  ]
  
  if (length(size_breaks_use) == 0) {
    size_breaks_use <- pretty_breaks(n = 4)(
      c(0, max_percent)
    )
    size_breaks_use <- size_breaks_use[size_breaks_use > 0]
  }
  
  fill_breaks <- pretty_breaks(n = 5)(
    range(plot_df$log10_abundance_percent, na.rm = TRUE)
  )
  
  p <- ggplot(
    plot_df,
    aes(
      x = Sample_Name,
      y = Species_Label,
      size = Relative_Abundance_Percent,
      fill = log10_abundance_percent
    )
  ) +
    geom_point(
      shape = 21,
      color = "black",
      stroke = 0.25,
      alpha = 0.95
    ) +
    scale_size_area(
      max_size = 8.5,
      breaks = size_breaks_use,
      limits = c(0, max_percent),
      labels = percent_label_clean,
      name = "Relative\nabundance"
    ) +
    scale_fill_viridis_c(
      option = "magma",
      direction = -1,
      breaks = fill_breaks,
      name = expression(log[10]~"(relative abundance, %)"),
      guide = guide_colorbar(
        order = 2,
        barheight = unit(4.2, "cm"),
        barwidth = unit(0.45, "cm")
      )
    ) +
    labs(
      title = title_text,
      x = "Sample",
      y = "Species"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      strip.background = element_rect(
        fill = "grey92",
        color = "black",
        linewidth = 0.4
      ),
      strip.text = element_text(
        face = "bold",
        size = 13
      ),
      axis.title.x = element_text(
        face = "bold",
        size = 14,
        margin = margin(t = 10)
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 14,
        margin = margin(r = 10)
      ),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 6.8,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(
        size = y_text_size,
        face = "bold.italic",
        color = "black",
        margin = margin(r = 4)
      ),
      panel.grid.major.x = element_line(
        color = "grey88",
        linewidth = 0.25
      ),
      panel.grid.major.y = element_line(
        color = "grey88",
        linewidth = 0.25
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.5
      ),
      legend.position = "right",
      legend.title = element_text(
        face = "bold",
        size = 11
      ),
      legend.text = element_text(size = 10),
      plot.margin = margin(t = 10, r = 20, b = 10, l = 10)
    ) +
    guides(
      size = guide_legend(
        order = 1,
        override.aes = list(
          fill = "grey70",
          color = "black",
          alpha = 1
        )
      )
    )
  
  if (facet_by_condition) {
    p <- p +
      facet_grid(
        . ~ Condition_Label,
        scales = "free_x",
        space = "free_x"
      )
  }
  
  p
}

# ============================================================
# 9. Generate top-species plot
# ============================================================

p_top <- make_species_dotplot(
  plot_data_top,
  title_text = NULL,
  y_text_size = 8.0,
  facet_by_condition = facet_by_condition
)

print(p_top)

ggsave(
  filename = out_top_png,
  plot = p_top,
  width = 15.5,
  height = 10.5,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = out_top_pdf,
  plot = p_top,
  width = 15.5,
  height = 10.5,
  bg = "white"
)

# ============================================================
# 10. Optional full all-species plot
# ============================================================

if (save_full_plot) {
  
  p_full <- make_species_dotplot(
    plot_data_full,
    title_text = NULL,
    y_text_size = 5.4,
    facet_by_condition = facet_by_condition
  )
  
  print(p_full)
  
  ggsave(
    filename = out_full_png,
    plot = p_full,
    width = 16.5,
    height = 18,
    dpi = 300,
    bg = "white",
    limitsize = FALSE
  )
  
  ggsave(
    filename = out_full_pdf,
    plot = p_full,
    width = 16.5,
    height = 18,
    bg = "white",
    limitsize = FALSE
  )
}

# ============================================================
# 11. Completion message
# ============================================================

cat("\nDone.\n")
cat("Top-species figure saved to:\n", out_top_png, "\n")
cat("Top-species PDF saved to:\n", out_top_pdf, "\n")

if (save_full_plot) {
  cat("Full all-species figure saved to:\n", out_full_png, "\n")
  cat("Full all-species PDF saved to:\n", out_full_pdf, "\n")
}

cat("Long-format plot data saved to:\n", out_csv, "\n")
cat("Species summary saved to:\n", out_summary, "\n")

cat("\nLegend note:\n")
cat("Dot size shows relative abundance as percent.\n")
cat("Dot fill shows log10(relative abundance percentage + pseudocount).\n")
cat("Using percent abundance avoids the confusing raw log10 scale where 0 is the maximum possible value for proportions.\n")

############################################################
############################################################
# FIGURE 6 MASTER SCRIPT — COMPLETE UPDATED VERSION

# Script also include Supplementary Figure S5 & S6
# Includes:
#   - Taxa features from vst_physeq
#   - Reaction features from ReactionAbundance.csv
#   - Subsystem features from SubsystemAbundance.csv
#   - Metabolite features from Metabolite_Clean_Summarized_Test.csv
#
# Main analyses:
#   - Spearman correlations with average PFAS elimination rate
#   - Feature network
#   - Repeated random forest prediction and stability analysis
#   - RF stable-feature correlation heatmap
#   - Supplementary RF QC and stability plots
#
# Important fixes:
#   - Handles mismatched sample columns across feature matrices.
#   - Robustly builds heatmap table without relying on Var1/Var2/Freq.
#   - Adds metabolite-label map for VMH-like metabolite IDs.
#   - Exports untranslated metabolite labels for manual curation.
#
# Requirements:
#   - vst_physeq already loaded in the R environment
############################################################
############################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggraph)
  library(igraph)
  library(tidygraph)
  library(scales)
  library(patchwork)
  library(readr)
  library(stringr)
  library(forcats)
  library(caret)
  library(ranger)
})

# ============================================================
# 1. PATHS
# ============================================================

reaction_path <- file.path("data", "model_outputs", "ReactionAbundance.csv")
subsystem_path <- file.path("data", "model_outputs", "SubsystemAbundance.csv")
metabolic_data_path <- file.path("data", "model_outputs", "Metabolite_Clean_Summarized_Test.csv")

out_dir <- file.path("results", "supplementary_figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_fig6_png <- file.path(out_dir, "Figure6_Final_A_to_D_WithMetabolites.png")
out_fig6_pdf <- file.path(out_dir, "Figure6_Final_A_to_D_WithMetabolites.pdf")

out_csv_all <- file.path(out_dir, "Figure6_PFAS_Elimination_Correlations_AllFeatures_WithMetabolites.csv")
out_csv_sig <- file.path(out_dir, "Figure6_PFAS_Elimination_Correlations_FilteredForPlot_WithMetabolites.csv")
out_csv_nodes <- file.path(out_dir, "Figure6_NetworkNodes_TopNPerType_WithMetabolites.csv")
out_csv_edges <- file.path(out_dir, "Figure6_NetworkEdges_TopNPerType_WithMetabolites.csv")

csv_perf <- file.path(out_dir, "Supp_RF_Performance_ByRun_WithMetabolites.csv")
csv_resid_stats <- file.path(out_dir, "Supp_RF_ResidualStats_ByRun_WithMetabolites.csv")
csv_importance_all <- file.path(out_dir, "Supp_RF_Importance_AllRuns_WithMetabolites.csv")
csv_selected_topk <- file.path(out_dir, "Supp_RF_Selected_Top25_ByRun_WithMetabolites.csv")
csv_stability <- file.path(out_dir, "Supp_RF_Stability_Summary_Top25Freq_WithMetabolites.csv")
csv_topstable_type <- file.path(out_dir, "Supp_RF_TopStableFeatures_PerType_Top25_WithMetabolites.csv")
csv_topdisplay <- file.path(out_dir, "Supp_RF_TopFeatures_DisplayedInStabilityFig_Top25_WithMetabolites.csv")
csv_stablecand <- file.path(out_dir, "Supp_RF_StableCandidates_Top25Freq_ge0.25_WithMetabolites.csv")
csv_stab_overall <- file.path(out_dir, "Supp_RF_Stability_OverallSummary_WithMetabolites.csv")
csv_stab_bytype <- file.path(out_dir, "Supp_RF_Stability_ByType_WithMetabolites.csv")
csv_rf_filtering <- file.path(out_dir, "Supp_RF_FeatureFiltering_Diagnostics_WithMetabolites.csv")
csv_rf_feature_heatmap <- file.path(out_dir, "Supp_RF_TopStableFeature_CorrelationHeatmapData_WithMetabolites.csv")

csv_metabolite_untranslated <- file.path(out_dir, "Figure6_MetaboliteLabels_StillNeedTranslation.csv")
txt_metabolite_template <- file.path(out_dir, "Figure6_MetaboliteLabels_ManualMap_Template.txt")

fig_qc_4panel <- file.path(out_dir, "SuppFig_RF_QC_100Runs_4Panel_WithMetabolites.png")
fig_qc_4panel_pdf <- file.path(out_dir, "SuppFig_RF_QC_100Runs_4Panel_WithMetabolites.pdf")

fig_stability_4panel <- file.path(out_dir, "SuppFig_RF_Stability_100Runs_4Panel_WithMetabolites.png")
fig_stability_4panel_pdf <- file.path(out_dir, "SuppFig_RF_Stability_100Runs_4Panel_WithMetabolites.pdf")

fig_rf_extra <- file.path(out_dir, "SuppFig_RF_TopStable_Features_Extra_WithMetabolites.png")
fig_rf_extra_pdf <- file.path(out_dir, "SuppFig_RF_TopStable_Features_Extra_WithMetabolites.pdf")

# ============================================================
# 2. PARAMETERS
# ============================================================

rho_threshold <- 0.12
pval_threshold <- 0.20
max_nodes_per_type <- 5
top_summary_per_type <- 10

min_shared_samples <- 4
min_pairs_cor <- min_shared_samples

seed_base <- 42
train_fraction <- 0.8
n_runs <- 100
top_k <- 25

num_trees <- 1000
min_node_size <- 5
cv_folds <- 5
cv_repeats <- 1

top_rf_features_per_type_main <- 5
top_rf_features_heatmap <- 18

debug_fast_rf <- FALSE
if (debug_fast_rf) {
  n_runs <- 5
  num_trees <- 250
}

show_plots_in_rstudio <- FALSE

feature_type_levels <- c("Taxa", "Reaction", "Subsystem", "Metabolite")
feature_type_levels_with_unknown <- c(feature_type_levels, "Unknown")

# ============================================================
# 3. MANUAL LABEL MAPS
# ============================================================

reaction_label_map <- c(
  "3OAACPR2" = "3-oxoacyl-ACP reductase",
  "3HACPR1" = "3-hydroxyacyl-ACP reductase",
  "4HTHRS" = "4-hydroxy-L-threonine synthase",
  "AAGAAH" = "N-acetylglucosamine aminohydrolase",
  "AMY2e" = "Extracellular amylase",
  "BACCL" = "Biotin carboxyl carrier protein ligase",
  "CEPA" = "Cephalosporin biosynthesis reaction",
  "EX core3(e)" = "Exchange of O-glycan core 3",
  "EX core4(e)" = "Exchange of O-glycan core 4",
  "EX core6(e)" = "Exchange of O-glycan core 6",
  "EX core7(e)" = "Exchange of O-glycan core 7",
  "EX glycogen2(e)" = "Exchange of glycogen structure 2",
  "GACPCD" = "Glutaryl-ACP decarboxylase",
  "GALASE OGLYCAN2epp" = "Periplasmic O-glycan galactosidase 2",
  "GALASE OGLYCAN3epp" = "Periplasmic O-glycan galactosidase 3",
  "GAMYe" = "Extracellular glucoamylase",
  "GLCNACASE OGLYCAN1epp" = "Periplasmic O-glycan N-acetylglucosaminidase 1",
  "GLYt4r" = "Reversible glycine transport",
  "Glycine transport" = "Glycine transport",
  "L-aspartate:nad[c]P+ oxidoreductase (deaminating)" = "L-aspartate NAD(P)+ oxidoreductase",
  "MALCOAMT" = "Malonyl-CoA methyltransferase",
  "RBK Dr" = "Ribokinase",
  "r1116" = "Uncharacterised reaction r1116",
  "SALCS2" = "Salicylate CoA synthase",
  "SIAASE OGLYCAN2epp" = "Periplasmic O-glycan sialidase 2",
  "sink s" = "Sink reaction",
  "SQLS" = "Squalene synthase",
  "S-Propane-1,2-diol facilitated transport" = "S-propane-1,2-diol facilitated transport",
  "T ANTIGENtex" = "T-antigen transport",
  "DM 4d 5" = "Demand reaction DM 4d 5",
  "DM 2HYM2EPH" = "Demand reaction DM 2HYM2EPH",
  "1,4-dihydroxy-2-naphthoate octaprenyltransferase" =
    "1,4-dihydroxy-2-naphthoate octaprenyltransferase",
  "1,4-dihydroxy-2-naphthoate Periplasm Transport" =
    "1,4-dihydroxy-2-naphthoate periplasmic transport",
  "Exchange of 1,4-Dihydroxy-2-naphthoate" =
    "Exchange of 1,4-dihydroxy-2-naphthoate",
  "10DMMCNFDOR" = "Dimethylmenaquinone:ferricytochrome oxidoreductase",
  "12DGR180t" = "1,2-diacylglycerol transport",
  "12PPD Stex" = "S-propane-1,2-diol extracellular transport",
  "12PPDtpp" = "S-propane-1,2-diol periplasmic transport",
  "13PPDH" = "Propanediol dehydrogenase",
  "13PPDtex" = "R-propane-1,3-diol extracellular transport",
  "13PPDtpp" = "R-propane-1,3-diol periplasmic transport",
  "15DAPtex" = "Diaminopentane extracellular transport",
  "1H2NPTH" = "1-hydroxy-2-naphthoate hydrolase",
  "1P4H2CBXLAH" = "1-pyrroline-4-hydroxy-2-carboxylate dehydrogenase",
  "22IDPOR" = "2,2-iminodipropanoate oxidoreductase",
  "23DAPAL" = "2,3-diaminopropionate ammonia-lyase",
  "24DCOAR" = "2,4-dienoyl-CoA reductase",
  "26DAPt2r" = "2,6-diaminopimelate reversible transport",
  "2DDPENTHL" = "2-dehydro-3-deoxy-D-pentonate hydrolase",
  "2DHPL" = "2-dehydropantoate aldolase",
  "2FBZOR" = "2-fluorobenzoate oxidoreductase",
  "2FBZTOR" = "2-fluorobenzoyl-CoA reductase",
  "2IMZS" = "2-isopropylmalate synthase",
  "2INSD" = "Inosine deaminase",
  "2IPDPIPT" = "Dimethylallyltranstransferase",
  "2OBUTFDXOR" = "2-oxobutyrate:ferredoxin oxidoreductase",
  "2OXOADOX" = "2-oxoadipate dehydrogenase",
  "33HPHACt2r" = "3-(3-hydroxyphenyl)propionate reversible transport",
  "34DHOXPEGOX" = "3,4-dihydroxyphenylglycol oxidase",
  "34DHPHAt2r" = "3,4-dihydroxyphenylacetate reversible transport",
  "35CGMPt2" = "3,5-cyclic GMP transport",
  "36DAHXI" = "3,6-dideoxy-D-arabino-hexose isomerase",
  "3DGUR" = "3-dehydroglucuronate reductase",
  "3FBZOR" = "3-fluorobenzoate oxidoreductase",
  "3FBZTOR" = "3-fluorobenzoyl-CoA reductase",
  "3FCHLOR" = "3-fluorocatechol chlorohydrolase",
  "3HACPR2" = "3-hydroxyacyl-ACP reductase 2",
  "3HAD100" = "3-hydroxyacyl-CoA dehydrogenase C10:0",
  "3HAD120" = "3-hydroxyacyl-CoA dehydrogenase C12:0",
  "3HAD140" = "3-hydroxyacyl-CoA dehydrogenase C14:0",
  "3HAD160" = "3-hydroxyacyl-CoA dehydrogenase C16:0",
  "3HAD161" = "3-hydroxyacyl-CoA dehydrogenase C16:1",
  "3HAD180" = "3-hydroxyacyl-CoA dehydrogenase C18:0",
  "3HAD40" = "3-hydroxyacyl-CoA dehydrogenase C4:0",
  "3HAD60" = "3-hydroxyacyl-CoA dehydrogenase C6:0",
  "3HAD80" = "3-hydroxyacyl-CoA dehydrogenase C8:0",
  "3MEACMPte" = "3-methyladenine cAMP transport",
  "3MOBDC2" = "3-methyl-2-oxobutanoate decarboxylase",
  "3MOPDC2" = "3-methyl-2-oxopentanoate decarboxylase",
  "3MOPLPAMO" = "3-methyl-2-oxopentanoate lipoamide oxidoreductase",
  "3OAACPR1" = "3-oxoacyl-ACP reductase 1",
  "3OADPCOAT" = "3-oxoadipyl-CoA thiolase",
  "3OPCPOOR" = "3-oxopropionyl-CoA oxidoreductase",
  "4ABUTtex" = "4-aminobutyrate extracellular transport",
  "4AHMMPtex" = "4-amino-5-hydroxymethyl-2-methylpyrimidine extracellular transport",
  "4FBZOR" = "4-fluorobenzoate oxidoreductase",
  "4FCHLOR" = "4-fluorocatechol chlorohydrolase",
  "4GBTAH" = "4-guanidinobutanoate amidinohydrolase",
  "4HBHYOXy" = "4-hydroxybenzoate hydroxylase",
  "4HBZCL" = "4-hydroxybenzoate CoA ligase"
)

# Fill or correct these manually as needed after checking:
# Figure6_MetaboliteLabels_StillNeedTranslation.csv
metabolite_label_map <- c(
  "galt[fe]" = "Galactitol exchange",
  "so4[fe]" = "Sulfate exchange",
  "actn R[u]" = "(R)-Acetoin",
  "cmp[fe]" = "CMP exchange",
  "glcn[fe]" = "Gluconate exchange",
  "h[fe]" = "Proton exchange",
  "co2528[ve]" = "Unknown Metabolite (co2528)",
  "glygly2[fe]" = "Glycylglycine exchange",
  "xylottr[u]" = "Xylotriose transport",
  "zn2[fe]" = "Zinc ion exchange",
  "dcholf[fe]" = "Deoxycholate exchange",
  "gbbtn[fe]" = "Gamma-butyrobetaine exchange",
  "mqn8[fe]" = "Menaquinone-8 exchange",
  "ca2[d]" = "Calcium ion exchange",
  "rblfrd[fe]" = "Riboflavin-derived metabolite exchange",
  "Tn antigen[fe]" = "Tn antigen exchange"
)

# ============================================================
# 4. HELPERS
# ============================================================

to_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "nd", "ND", "na", "n/a", "N/A", "<LOD", "<LOQ", "< LOQ")] <- NA
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- gsub("^<\\s*", "", x)
  suppressWarnings(as.numeric(x))
}

clean_label <- function(x, max_chars = Inf) {
  x <- as.character(x)
  x <- str_replace_all(x, "^Taxa\\.Taxa\\:\\:", "")
  x <- str_replace_all(x, "^Reaction\\.Reaction\\:\\:", "")
  x <- str_replace_all(x, "^Subsystem\\.Subsystem\\:\\:", "")
  x <- str_replace_all(x, "^Metabolite\\.Metabolite\\:\\:", "")
  x <- str_replace_all(x, "^Taxa\\:\\:", "")
  x <- str_replace_all(x, "^Reaction\\:\\:", "")
  x <- str_replace_all(x, "^Subsystem\\:\\:", "")
  x <- str_replace_all(x, "^Metabolite\\:\\:", "")
  x <- str_replace_all(x, "^(Taxa|OTU|Reaction|Subsystem|Metabolite)[\\.\\:\\_\\-\\s]+", "")
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_trim(x)
  
  if (is.finite(max_chars)) {
    x <- str_trunc(x, width = max_chars)
  }
  
  x
}

standardise_reaction_key <- function(x) {
  x <- clean_label(x, max_chars = Inf)
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "\\s+", " ")
  str_trim(x)
}

humanise_feature_label <- function(label, type = NULL, max_chars = Inf) {
  label_clean <- clean_label(label, max_chars = Inf)
  
  if (!is.null(type)) {
    type_chr <- as.character(type)
  } else {
    type_chr <- rep(NA_character_, length(label_clean))
  }
  
  out <- label_clean
  
  reaction_key <- standardise_reaction_key(label_clean)
  reaction_mapped <- reaction_label_map[reaction_key]
  
  reaction_idx <- is.na(type_chr) | type_chr == "Reaction" | str_detect(label_clean, "^Reaction")
  out[reaction_idx & !is.na(reaction_mapped)] <- unname(
    reaction_mapped[reaction_idx & !is.na(reaction_mapped)]
  )
  
  metabolite_idx <- !is.na(type_chr) & type_chr == "Metabolite"
  metabolite_key <- label_clean
  metabolite_mapped <- metabolite_label_map[metabolite_key]
  out[metabolite_idx & !is.na(metabolite_mapped)] <- unname(
    metabolite_mapped[metabolite_idx & !is.na(metabolite_mapped)]
  )
  
  out <- str_replace_all(out, "^EX ([A-Za-z0-9]+)\\(e\\)$", "Exchange of \\1")
  out <- str_replace_all(out, "^DM ([A-Za-z0-9 ]+)$", "Demand reaction \\1")
  out <- str_replace_all(out, "^sink ([A-Za-z0-9 ]+)$", "Sink reaction \\1")
  
  clean_label(out, max_chars = max_chars)
}

make_unique_safe <- function(x) {
  make.unique(as.character(x), sep = "_dup")
}

theme_clean <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", color = "black"),
      axis.text = element_text(face = "bold", color = "black"),
      legend.title = element_text(face = "bold", color = "black"),
      legend.text = element_text(color = "black"),
      strip.text = element_text(face = "bold", color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
      panel.grid.minor = element_blank()
    )
}

safe_print <- function(p) {
  if (isTRUE(show_plots_in_rstudio)) {
    tryCatch(
      print(p),
      error = function(e) {
        message("Plot was saved, but RStudio plot pane could not display it: ", e$message)
      }
    )
  }
}

get_taxa_labels <- function(physeq_obj, otu_rownames) {
  if (!is.null(tax_table(physeq_obj, errorIfNULL = FALSE))) {
    tax_df <- as.data.frame(tax_table(physeq_obj))
    tax_df$TaxonID <- rownames(tax_df)
    
    if ("Species" %in% colnames(tax_df)) {
      labels <- as.character(tax_df$Species)
      bad <- is.na(labels) | labels == "" | labels == "NA"
      labels[bad] <- tax_df$TaxonID[bad]
    } else {
      labels <- tax_df$TaxonID
    }
    
    names(labels) <- tax_df$TaxonID
    labels <- labels[otu_rownames]
    labels[is.na(labels)] <- otu_rownames[is.na(labels)]
    return(labels)
  }
  
  otu_rownames
}

is_vmh_like_label <- function(x) {
  x <- clean_label(x, max_chars = Inf)
  
  short_code <- nchar(x) <= 16 &
    str_detect(x, "^[A-Za-z0-9_\\-\\(\\)\\[\\]\\+]+$") &
    !str_detect(x, "\\s")
  
  vmh_exchange <- str_detect(x, "^EX[_\\s]")
  vmh_demand <- str_detect(x, "^DM[_\\s]")
  vmh_sink <- str_detect(x, "^sink[_\\s]")
  r_code <- str_detect(x, "^r[0-9]+$")
  epp_or_tex <- str_detect(x, "epp$|tex$")
  
  short_code | vmh_exchange | vmh_demand | vmh_sink | r_code | epp_or_tex
}

is_vmh_metabolite_like_label <- function(x) {
  x <- clean_label(x, max_chars = Inf)
  
  str_detect(x, "^[A-Za-z0-9_]+\\[[a-zA-Z]+\\]$") |
    str_detect(x, "^[A-Za-z0-9_]+\\([a-zA-Z]+\\)$") |
    (
      nchar(x) <= 14 &
        str_detect(x, "^[A-Za-z0-9_]+$") &
        !str_detect(x, "\\s")
    )
}

get_plotted_reaction_labels_needing_translation <- function(plotted_df, out_dir = NULL) {
  if (is.null(plotted_df) || nrow(plotted_df) == 0) {
    return(tibble())
  }
  
  out <- plotted_df %>%
    mutate(
      type = as.character(type),
      original_label = as.character(feature_label),
      human_label = humanise_feature_label(original_label, type = type, max_chars = Inf),
      still_cryptic = type == "Reaction" & is_vmh_like_label(human_label)
    ) %>%
    filter(type == "Reaction", still_cryptic) %>%
    distinct(source, type, model_feature, feature_key, original_label, human_label) %>%
    arrange(source, original_label)
  
  if (!is.null(out_dir)) {
    write_csv(
      out,
      file.path(out_dir, "Figure6_PLOTTED_ReactionLabels_StillNeedTranslation_WithMetabolites.csv")
    )
    
    if (nrow(out) > 0) {
      writeLines(
        paste0('"', out$original_label, '"', collapse = ",\n"),
        file.path(out_dir, "Figure6_PLOTTED_ReactionLabels_StillNeedTranslation_CopyPaste_WithMetabolites.txt")
      )
    }
  }
  
  out
}

compute_feature_cor <- function(feature_matrix, rates_vec, min_pairs = 4) {
  features <- rownames(feature_matrix)
  
  res_list <- lapply(features, function(f) {
    x <- to_numeric_safe(feature_matrix[f, ])
    y <- rates_vec[colnames(feature_matrix)]
    
    good <- !is.na(x) & !is.na(y)
    n_pairs <- sum(good)
    
    if (n_pairs < min_pairs) return(NULL)
    if (var(x[good], na.rm = TRUE) == 0) return(NULL)
    
    ct <- suppressWarnings(
      tryCatch(
        cor.test(x[good], y[good], method = "spearman", exact = FALSE),
        error = function(e) NULL
      )
    )
    
    if (is.null(ct)) return(NULL)
    
    data.frame(
      feature_key = f,
      rho = as.numeric(ct$estimate),
      p = ct$p.value,
      n_pairs = n_pairs,
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(res_list)
}

guess_id_col <- function(dt, preferred_names = NULL) {
  if (!is.null(preferred_names)) {
    hit <- preferred_names[preferred_names %in% colnames(dt)]
    if (length(hit) > 0) return(hit[1])
  }
  colnames(dt)[1]
}

prepare_feature_df <- function(dt, id_col = NULL, keep_samples, type_label, min_shared_samples = 4) {
  if (is.null(id_col) || !(id_col %in% colnames(dt))) {
    id_col <- colnames(dt)[1]
  }
  
  samples <- intersect(colnames(dt), keep_samples)
  
  if (length(samples) < min_shared_samples) {
    warning(
      "Skipping ", type_label,
      ": fewer than ", min_shared_samples,
      " shared samples found."
    )
    return(NULL)
  }
  
  df <- dt[, c(id_col, samples), drop = FALSE]
  
  feature_label_raw <- as.character(df[[id_col]])
  feature_label <- clean_label(feature_label_raw, max_chars = Inf)
  feature_label <- humanise_feature_label(feature_label, type = type_label, max_chars = Inf)
  
  feature_key <- paste0(type_label, "::", make_unique_safe(feature_label))
  
  rownames(df) <- feature_key
  df[[id_col]] <- NULL
  
  df <- as.data.frame(df, check.names = FALSE)
  df[] <- lapply(df, to_numeric_safe)
  
  lookup <- tibble(
    feature_key = feature_key,
    feature_label = feature_label,
    type = type_label
  )
  
  list(matrix = df[, samples, drop = FALSE], lookup = lookup)
}

median_impute_matrix <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  
  for (j in seq_len(ncol(x))) {
    v <- x[, j]
    if (all(is.na(v))) {
      v[] <- 0
    } else {
      med <- median(v, na.rm = TRUE)
      v[is.na(v)] <- med
    }
    x[, j] <- v
  }
  
  x
}

build_rf_matrix <- function(feature_mats, feature_lookup, Y, min_non_missing = 2) {
  
  feature_mats <- feature_mats[!vapply(feature_mats, is.null, logical(1))]
  
  feature_mats <- lapply(feature_mats, function(m) {
    m <- as.data.frame(m, check.names = FALSE)
    m[] <- lapply(m, to_numeric_safe)
    m
  })
  
  shared_samples <- Reduce(intersect, lapply(feature_mats, colnames))
  shared_samples <- intersect(shared_samples, names(Y))
  
  if (length(shared_samples) < 4) {
    stop(
      "Fewer than 4 shared samples across feature matrices and Y. ",
      "Check sample names in taxa/reaction/subsystem/metabolite matrices."
    )
  }
  
  message("\nRF build: shared samples across all feature classes and Y: ", length(shared_samples))
  
  feature_mats <- lapply(feature_mats, function(m) {
    m[, shared_samples, drop = FALSE]
  })
  
  features_df <- do.call(rbind, unname(feature_mats))
  features_df <- as.data.frame(features_df, check.names = FALSE)
  features_df[] <- lapply(features_df, to_numeric_safe)
  
  rf_lookup <- tibble(feature_key = rownames(features_df)) %>%
    left_join(feature_lookup, by = "feature_key") %>%
    mutate(
      feature_label = coalesce(feature_label, feature_key),
      type = coalesce(type, case_when(
        str_detect(feature_key, "^Taxa::") ~ "Taxa",
        str_detect(feature_key, "^Reaction::") ~ "Reaction",
        str_detect(feature_key, "^Subsystem::") ~ "Subsystem",
        str_detect(feature_key, "^Metabolite::") ~ "Metabolite",
        TRUE ~ "Unknown"
      )),
      feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
      model_feature = paste0("F", sprintf("%05d", row_number()))
    )
  
  rownames(features_df) <- rf_lookup$model_feature
  
  X_raw <- t(as.matrix(features_df))
  storage.mode(X_raw) <- "numeric"
  
  Y2 <- Y[rownames(X_raw)]
  
  valid_y <- !is.na(Y2) & is.finite(Y2)
  X_raw <- X_raw[valid_y, , drop = FALSE]
  Y2 <- Y2[valid_y]
  
  feature_non_missing <- colSums(!is.na(X_raw))
  feature_sd_pre <- apply(X_raw, 2, sd, na.rm = TRUE)
  
  keep_pre <- feature_non_missing >= min_non_missing
  X_pre <- X_raw[, keep_pre, drop = FALSE]
  
  if (ncol(X_pre) < 2) {
    warning("RF pre-filter retained fewer than 2 features. Relaxing to any non-missing feature.")
    keep_pre <- feature_non_missing >= 1
    X_pre <- X_raw[, keep_pre, drop = FALSE]
  }
  
  X_imp <- median_impute_matrix(X_pre)
  
  feature_sd_post <- apply(X_imp, 2, sd, na.rm = TRUE)
  keep_post <- is.finite(feature_sd_post) & feature_sd_post > 0
  
  X_imp <- X_imp[, keep_post, drop = FALSE]
  
  rf_lookup <- rf_lookup %>%
    filter(model_feature %in% colnames(X_imp)) %>%
    mutate(
      feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
      type = factor(type, levels = feature_type_levels_with_unknown)
    )
  
  filter_diag <- tibble(
    model_feature = colnames(X_raw),
    non_missing_n = feature_non_missing,
    sd_pre_imputation = feature_sd_pre,
    kept_pre_imputation = colnames(X_raw) %in% colnames(X_pre),
    kept_post_imputation = colnames(X_raw) %in% colnames(X_imp)
  ) %>%
    left_join(
      rf_lookup %>% select(model_feature, feature_key, feature_label, type),
      by = "model_feature"
    )
  
  list(
    X = as.data.frame(X_imp, check.names = FALSE),
    Y = Y2,
    rf_lookup = rf_lookup,
    filter_diag = filter_diag,
    n_features_raw = ncol(X_raw),
    n_features_after = ncol(X_imp),
    n_samples = nrow(X_imp),
    shared_samples = rownames(X_imp)
  )
}

# ============================================================
# 5. CHECK INPUT OBJECT
# ============================================================

if (!exists("vst_physeq")) {
  stop("vst_physeq not found. Load or compute it first.")
}

vst_physeq_subset <- subset_samples(
  vst_physeq,
  Condition %in% c("PFAS_Exposed", "Reference")
)

# ============================================================
# 6. COMPUTE AVERAGE PFAS ELIMINATION RATE
# ============================================================

pfas_cols <- c(
  "k_PFOA", "k_PFPeS", "k_PFHxS", "k_PFHpS",
  "k_LPFOS", "k_PFOS_MP1", "k_PFOS_MP345", "k_PFOS_MP26"
)

sdat <- as.data.frame(sample_data(vst_physeq_subset), check.names = FALSE)

present_pfas_cols <- intersect(pfas_cols, colnames(sdat))
if (length(present_pfas_cols) == 0) {
  stop("No PFAS elimination-rate columns found in sample_data().")
}

for (cc in present_pfas_cols) {
  sdat[[cc]] <- to_numeric_safe(sdat[[cc]])
}

sdat$average_k <- rowMeans(
  sdat[, present_pfas_cols, drop = FALSE],
  na.rm = TRUE
)

sdat$average_k[
  rowSums(!is.na(sdat[, present_pfas_cols, drop = FALSE])) == 0
] <- NA_real_

sample_data(vst_physeq_subset)$average_k <- sdat$average_k

pfas_rates <- as.numeric(sample_data(vst_physeq_subset)$average_k)
names(pfas_rates) <- sample_names(vst_physeq_subset)
pfas_rates <- pfas_rates[!is.na(pfas_rates)]

if (length(pfas_rates) < min_shared_samples) {
  stop("Too few samples with non-missing average_k.")
}

keep_samples <- names(pfas_rates)

message("\nSamples with non-missing average_k: ", length(keep_samples))

# ============================================================
# 7. FEATURE MATRICES
# ============================================================

otu_tab <- as.data.frame(otu_table(vst_physeq_subset), check.names = FALSE)

if (!taxa_are_rows(vst_physeq_subset)) {
  otu_tab <- t(otu_tab)
  otu_tab <- as.data.frame(otu_tab, check.names = FALSE)
}

otu_tab <- otu_tab[, keep_samples, drop = FALSE]
otu_tab[] <- lapply(otu_tab, to_numeric_safe)

taxa_labels <- get_taxa_labels(vst_physeq_subset, rownames(otu_tab))
taxa_labels <- humanise_feature_label(taxa_labels, type = "Taxa", max_chars = Inf)

taxa_feature_key <- paste0("Taxa::", make_unique_safe(taxa_labels))
rownames(otu_tab) <- taxa_feature_key

taxa_lookup <- tibble(
  feature_key = taxa_feature_key,
  feature_label = taxa_labels,
  type = "Taxa"
)

reaction_dt <- fread(reaction_path, data.table = FALSE, check.names = FALSE)
subsystem_dt <- fread(subsystem_path, data.table = FALSE, check.names = FALSE)
metabolite_dt <- fread(metabolic_data_path, data.table = FALSE, check.names = FALSE)

reaction_id_col <- guess_id_col(
  reaction_dt,
  preferred_names = c("Reaction", "Reactions", "reaction", "reaction_id")
)

subsystem_id_col <- guess_id_col(
  subsystem_dt,
  preferred_names = c("Subsystems", "Subsystem", "subsystem")
)

metabolite_id_col <- guess_id_col(
  metabolite_dt,
  preferred_names = c("Metabolite", "Metabolites", "metabolite", "metabolite_id", "Name", "ID")
)

reaction_prepared <- prepare_feature_df(
  reaction_dt,
  id_col = reaction_id_col,
  keep_samples = keep_samples,
  type_label = "Reaction",
  min_shared_samples = min_shared_samples
)

subsystem_prepared <- prepare_feature_df(
  subsystem_dt,
  id_col = subsystem_id_col,
  keep_samples = keep_samples,
  type_label = "Subsystem",
  min_shared_samples = min_shared_samples
)

metabolite_prepared <- prepare_feature_df(
  metabolite_dt,
  id_col = metabolite_id_col,
  keep_samples = keep_samples,
  type_label = "Metabolite",
  min_shared_samples = min_shared_samples
)

feature_mats <- list(Taxa = otu_tab)
feature_lookup <- taxa_lookup

if (!is.null(reaction_prepared)) {
  feature_mats$Reaction <- reaction_prepared$matrix
  feature_lookup <- bind_rows(feature_lookup, reaction_prepared$lookup)
}

if (!is.null(subsystem_prepared)) {
  feature_mats$Subsystem <- subsystem_prepared$matrix
  feature_lookup <- bind_rows(feature_lookup, subsystem_prepared$lookup)
}

if (!is.null(metabolite_prepared)) {
  feature_mats$Metabolite <- metabolite_prepared$matrix
  feature_lookup <- bind_rows(feature_lookup, metabolite_prepared$lookup)
}

message("\nFeature matrix sample-column diagnostics before RF build:")
for (nm in names(feature_mats)) {
  message("  ", nm, ": ", ncol(feature_mats[[nm]]), " sample columns")
}
message("  Shared across all feature matrices: ", length(Reduce(intersect, lapply(feature_mats, colnames))))

feature_lookup <- feature_lookup %>%
  distinct(feature_key, .keep_all = TRUE) %>%
  mutate(
    feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    type = factor(type, levels = feature_type_levels)
  )

message("\nFeatures loaded:")
message("  Taxa: ", nrow(feature_mats$Taxa))
if (!is.null(feature_mats$Reaction)) message("  Reactions: ", nrow(feature_mats$Reaction))
if (!is.null(feature_mats$Subsystem)) message("  Subsystems: ", nrow(feature_mats$Subsystem))
if (!is.null(feature_mats$Metabolite)) message("  Metabolites: ", nrow(feature_mats$Metabolite))

# ============================================================
# 8. SPEARMAN CORRELATION ANALYSIS
# ============================================================

cor_list <- lapply(names(feature_mats), function(tp) {
  mat_i <- feature_mats[[tp]]
  
  samples_i <- intersect(colnames(mat_i), names(pfas_rates))
  mat_i <- mat_i[, samples_i, drop = FALSE]
  
  compute_feature_cor(mat_i, pfas_rates, min_pairs = min_pairs_cor) %>%
    mutate(type = tp)
})

all_cor <- bind_rows(cor_list) %>%
  left_join(feature_lookup, by = c("feature_key", "type")) %>%
  mutate(
    feature_label = coalesce(feature_label, feature_key),
    feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    abs_rho = abs(rho),
    p_adj_BH = p.adjust(p, method = "BH"),
    direction = ifelse(rho > 0, "Positive", "Negative"),
    type = factor(type, levels = feature_type_levels)
  ) %>%
  arrange(type, desc(abs_rho))

write_csv(all_cor, out_csv_all)

sig_all <- all_cor %>%
  filter(abs_rho >= rho_threshold, p <= pval_threshold)

write_csv(sig_all, out_csv_sig)

if (nrow(sig_all) == 0) {
  stop("No features passed the rho/p-value thresholds for Figure 6.")
}

message("\nCorrelation features passing plot filter: ", nrow(sig_all))

# ============================================================
# 9. NETWORK NODES AND EDGES
# ============================================================

net_nodes <- sig_all %>%
  group_by(type) %>%
  slice_max(abs_rho, n = max_nodes_per_type, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(node_name = make_unique_safe(feature_label))

nodes <- distinct(bind_rows(
  tibble(
    name = net_nodes$node_name,
    label = net_nodes$feature_label,
    type = as.character(net_nodes$type)
  ),
  tibble(
    name = "PFAS_Elimination_Rate",
    label = "PFAS elimination rate",
    type = "PFAS"
  )
))

edges <- tibble(
  from = "PFAS_Elimination_Rate",
  to = net_nodes$node_name,
  weight = net_nodes$abs_rho,
  rho = net_nodes$rho,
  p = net_nodes$p,
  p_adj_BH = net_nodes$p_adj_BH,
  n_pairs = net_nodes$n_pairs,
  sign = ifelse(net_nodes$rho > 0, "pos", "neg")
)

write_csv(nodes, out_csv_nodes)
write_csv(edges, out_csv_edges)

node_stats <- edges %>%
  group_by(to) %>%
  summarise(
    mean_abs_rho = mean(abs(rho), na.rm = TRUE),
    mean_p = mean(p, na.rm = TRUE),
    node_size = -log10(mean_p + 1e-10),
    .groups = "drop"
  )

net <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

set.seed(123)
layout_network <- create_layout(net, layout = "fr") %>%
  left_join(node_stats, by = c("name" = "to")) %>%
  mutate(node_size = ifelse(is.na(node_size), 2.5, node_size))

# ============================================================
# 10. COLORS
# ============================================================

col_type <- c(
  "Taxa" = "purple3",
  "Reaction" = "deepskyblue3",
  "Subsystem" = "slateblue3",
  "Metabolite" = "darkorange2",
  "PFAS" = "red3",
  "Unknown" = "grey60"
)

col_edge_sign <- c(
  "pos" = "#B22222",
  "neg" = "#2E8B57"
)

# ============================================================
# 11. MAIN FIGURE 6A/B
# ============================================================

p_summary <- sig_all %>%
  group_by(type) %>%
  slice_max(abs_rho, n = top_summary_per_type, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    type = factor(type, levels = feature_type_levels),
    feature_label_plot = humanise_feature_label(feature_label, type = type, max_chars = 65),
    feature_label_plot = fct_reorder(feature_label_plot, abs_rho)
  ) %>%
  ggplot(aes(x = feature_label_plot, y = abs_rho, fill = type)) +
  geom_col(alpha = 0.95, colour = "black", width = 0.54) +
  coord_flip() +
  facet_wrap(~ type, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = col_type, name = "Feature class", drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    x = "Feature",
    y = expression("|Spearman " * rho * "|")
  ) +
  theme_clean(base_size = 15) +
  theme(
    axis.text.y = element_text(size = 11.2, face = "bold", margin = margin(r = 6)),
    axis.text.x = element_text(size = 12.5, face = "bold"),
    axis.title.y = element_text(size = 14, margin = margin(r = 12)),
    axis.title.x = element_text(size = 14),
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14),
    panel.grid.major.y = element_blank(),
    panel.spacing.y = unit(1.6, "lines"),
    plot.margin = margin(12, 10, 12, 10)
  )

p_network <- ggraph(layout_network) +
  geom_edge_link(
    aes(width = weight, colour = sign),
    alpha = 0.8,
    lineend = "round"
  ) +
  scale_edge_color_manual(
    values = col_edge_sign,
    name = "Correlation direction",
    labels = c(
      "neg" = "Negative (ρ < 0)",
      "pos" = "Positive (ρ > 0)"
    )
  ) +
  geom_node_point(
    aes(color = type, size = node_size),
    alpha = 0.95
  ) +
  geom_node_text(
    aes(label = label),
    repel = TRUE,
    size = 3.7,
    color = "black",
    fontface = "bold",
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = col_type,
    name = "Feature class",
    labels = c(
      "PFAS" = "PFAS elimination rate",
      "Reaction" = "Metabolic reaction",
      "Subsystem" = "Functional subsystem",
      "Taxa" = "Microbial taxon",
      "Metabolite" = "Model-predicted metabolite"
    )
  ) +
  scale_size_continuous(
    name = expression(-log[10](italic(p))),
    range = c(3.2, 9.0)
  ) +
  scale_edge_width(
    range = c(0.60, 3.4),
    name = "|ρ|"
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 5)),
    size = guide_legend(order = 2),
    edge_colour = guide_legend(order = 3),
    edge_width = guide_legend(order = 4)
  ) +
  theme_void() +
  theme(
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 11.5),
    legend.position = "right",
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================
# 12. BUILD RF MATRIX
# ============================================================

rf_build <- build_rf_matrix(
  feature_mats = feature_mats,
  feature_lookup = feature_lookup,
  Y = pfas_rates,
  min_non_missing = 2
)

X <- rf_build$X
Y <- rf_build$Y
rf_lookup <- rf_build$rf_lookup

write_csv(rf_build$filter_diag, csv_rf_filtering)

message("\nRF matrix diagnostic:")
message("  Samples: ", rf_build$n_samples)
message("  Raw features before RF filtering: ", rf_build$n_features_raw)
message("  Features after imputation and variance filtering: ", rf_build$n_features_after)
message("  RF refits planned: ", n_runs)

if (ncol(X) < 2) {
  stop(
    paste0(
      "RF still has fewer than 2 non-zero variance features. ",
      "Check diagnostic CSV: ", csv_rf_filtering
    )
  )
}

# ============================================================
# 13. RF RUN FUNCTION
# ============================================================

run_rf_once <- function(run_id, X, Y, train_fraction, seed,
                        cv_folds, cv_repeats,
                        num_trees, min_node_size) {
  
  set.seed(seed)
  
  train_idx <- createDataPartition(Y, p = train_fraction, list = FALSE)
  
  X_train <- X[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx]
  
  X_test <- X[-train_idx, , drop = FALSE]
  Y_test <- Y[-train_idx]
  
  train_ctrl <- trainControl(
    method = if (cv_repeats > 1) "repeatedcv" else "cv",
    number = cv_folds,
    repeats = cv_repeats
  )
  
  tgrid <- expand.grid(
    mtry = max(1, floor(sqrt(ncol(X_train)))),
    splitrule = "variance",
    min.node.size = min_node_size
  )
  
  rf_model <- train(
    x = X_train,
    y = Y_train,
    method = "ranger",
    trControl = train_ctrl,
    tuneGrid = tgrid,
    importance = "permutation",
    num.trees = num_trees,
    metric = "RMSE"
  )
  
  pred_test <- predict(rf_model, X_test)
  resid_test <- as.numeric(Y_test - pred_test)
  
  rmse_test <- sqrt(mean(resid_test^2, na.rm = TRUE))
  r2_test <- suppressWarnings(cor(Y_test, pred_test, use = "complete.obs")^2)
  
  sd_resid <- sd(resid_test, na.rm = TRUE)
  std_resid <- if (!is.finite(sd_resid) || sd_resid == 0) {
    rep(NA_real_, length(resid_test))
  } else {
    resid_test / sd_resid
  }
  
  shapiro_p <- if (length(resid_test) >= 3 && length(resid_test) <= 5000) {
    suppressWarnings(shapiro.test(resid_test)$p.value)
  } else {
    NA_real_
  }
  
  perf <- tibble(
    run = run_id,
    seed = seed,
    n_train = nrow(X_train),
    n_test = nrow(X_test),
    rmse = rmse_test,
    r2 = as.numeric(r2_test)
  )
  
  resid_stats <- tibble(
    run = run_id,
    resid_mean = mean(resid_test, na.rm = TRUE),
    resid_sd = sd(resid_test, na.rm = TRUE),
    shapiro_p = shapiro_p
  )
  
  var_imp <- varImp(rf_model, scale = TRUE)$importance
  
  imp_tbl <- tibble(
    run = run_id,
    model_feature = rownames(var_imp),
    importance = as.numeric(var_imp[[1]])
  ) %>%
    arrange(desc(importance))
  
  test_df <- tibble(
    run = run_id,
    sample = rownames(X_test),
    obs = as.numeric(Y_test),
    pred = as.numeric(pred_test),
    resid = resid_test,
    std_resid = std_resid
  )
  
  list(
    perf = perf,
    resid_stats = resid_stats,
    importance = imp_tbl,
    test_df = test_df
  )
}

# ============================================================
# 14. RUN RF REPEATED REFITS
# ============================================================

message("\n============================================================")
message("Starting RF repeated-refit analysis")
message("Total RF refits: ", n_runs)
message("============================================================")
flush.console()

run_results <- vector("list", n_runs)

for (i in seq_len(n_runs)) {
  message(sprintf("RF refit %03d / %03d", i, n_runs))
  flush.console()
  
  run_results[[i]] <- run_rf_once(
    run_id = i,
    X = X,
    Y = Y,
    train_fraction = train_fraction,
    seed = seed_base + i,
    cv_folds = cv_folds,
    cv_repeats = cv_repeats,
    num_trees = num_trees,
    min_node_size = min_node_size
  )
}

message("Finished RF repeated-refit analysis.")
flush.console()

perf_df <- bind_rows(lapply(run_results, `[[`, "perf"))
resid_stats_df <- bind_rows(lapply(run_results, `[[`, "resid_stats"))
test_all_df <- bind_rows(lapply(run_results, `[[`, "test_df"))

imp_all_df <- bind_rows(lapply(run_results, `[[`, "importance")) %>%
  left_join(
    rf_lookup %>%
      select(model_feature, feature_key, feature_label, type),
    by = "model_feature"
  ) %>%
  mutate(
    feature_label = coalesce(feature_label, model_feature),
    feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    type = coalesce(as.character(type), "Unknown"),
    type = factor(type, levels = feature_type_levels_with_unknown)
  )

# ============================================================
# 15. RF STABILITY SUMMARY
# ============================================================

imp_topk_df <- imp_all_df %>%
  group_by(run) %>%
  arrange(desc(importance), .by_group = TRUE) %>%
  slice_head(n = top_k) %>%
  ungroup() %>%
  mutate(selected_topk = 1L)

stability_df <- imp_all_df %>%
  group_by(model_feature, feature_key, feature_label, type) %>%
  summarise(
    mean_importance = mean(importance, na.rm = TRUE),
    sd_importance = sd(importance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    imp_topk_df %>%
      group_by(model_feature) %>%
      summarise(topk_count = sum(selected_topk), .groups = "drop"),
    by = "model_feature"
  ) %>%
  mutate(
    topk_count = replace_na(topk_count, 0L),
    topk_freq = topk_count / n_runs,
    feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    type = factor(type, levels = feature_type_levels_with_unknown)
  ) %>%
  arrange(desc(topk_freq), desc(mean_importance))

top_stable_per_type <- stability_df %>%
  group_by(type) %>%
  slice_max(order_by = topk_freq, n = 10, with_ties = FALSE) %>%
  ungroup()

top_display <- stability_df %>%
  filter(type %in% feature_type_levels) %>%
  group_by(type) %>%
  slice_max(order_by = topk_freq, n = top_rf_features_per_type_main, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(type, desc(topk_freq), desc(mean_importance))

stable_candidates <- stability_df %>%
  filter(topk_freq >= 0.25) %>%
  arrange(desc(topk_freq), desc(mean_importance))

# ============================================================
# 15B. METABOLITE TRANSLATION DIAGNOSTIC
# ============================================================

metabolite_labels_needing_translation <- stability_df %>%
  filter(type == "Metabolite") %>%
  mutate(
    current_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    needs_translation = is_vmh_metabolite_like_label(current_label)
  ) %>%
  filter(needs_translation) %>%
  distinct(
    model_feature,
    feature_key,
    current_label,
    topk_freq,
    mean_importance
  ) %>%
  arrange(desc(topk_freq), desc(mean_importance), current_label)

write_csv(
  metabolite_labels_needing_translation,
  csv_metabolite_untranslated
)

if (nrow(metabolite_labels_needing_translation) > 0) {
  writeLines(
    paste0('"', metabolite_labels_needing_translation$current_label, '" = "",'),
    txt_metabolite_template
  )
} else {
  writeLines(
    "No VMH-like metabolite labels detected after current mapping.",
    txt_metabolite_template
  )
}

cat("\n============================================================\n")
cat("METABOLITE LABELS STILL NEEDING TRANSLATION\n")
cat("============================================================\n")
print(metabolite_labels_needing_translation, n = Inf)

# ============================================================
# 16. WRITE RF CSV OUTPUTS
# ============================================================

write_csv(perf_df, csv_perf)
write_csv(resid_stats_df, csv_resid_stats)
write_csv(imp_all_df, csv_importance_all)
write_csv(imp_topk_df, csv_selected_topk)
write_csv(stability_df, csv_stability)
write_csv(top_stable_per_type, csv_topstable_type)
write_csv(top_display, csv_topdisplay)
write_csv(stable_candidates, csv_stablecand)

supp_stability_summary <- tibble(
  n_runs = n_runs,
  train_fraction = train_fraction,
  top_k = top_k,
  rmse_median = median(perf_df$rmse, na.rm = TRUE),
  rmse_iqr25 = quantile(perf_df$rmse, 0.25, na.rm = TRUE),
  rmse_iqr75 = quantile(perf_df$rmse, 0.75, na.rm = TRUE),
  r2_median = median(perf_df$r2, na.rm = TRUE),
  r2_iqr25 = quantile(perf_df$r2, 0.25, na.rm = TRUE),
  r2_iqr75 = quantile(perf_df$r2, 0.75, na.rm = TRUE),
  n_features_total = nrow(stability_df),
  n_features_selected_at_least_once = sum(stability_df$topk_count > 0, na.rm = TRUE),
  n_features_selected_in_ge10pct_runs = sum(stability_df$topk_freq >= 0.10, na.rm = TRUE),
  n_features_selected_in_ge25pct_runs = sum(stability_df$topk_freq >= 0.25, na.rm = TRUE),
  n_features_selected_in_ge50pct_runs = sum(stability_df$topk_freq >= 0.50, na.rm = TRUE)
)

write_csv(supp_stability_summary, csv_stab_overall)

supp_stability_by_type <- stability_df %>%
  group_by(type) %>%
  summarise(
    n_features = n(),
    n_selected_ge10pct = sum(topk_freq >= 0.10, na.rm = TRUE),
    n_selected_ge25pct = sum(topk_freq >= 0.25, na.rm = TRUE),
    n_selected_ge50pct = sum(topk_freq >= 0.50, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(supp_stability_by_type, csv_stab_bytype)

# ============================================================
# 17. MAIN FIGURE 6C — TOP STABLE RF PREDICTORS PER TYPE
# ============================================================

main_rf_plot_df <- top_display %>%
  mutate(
    feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    feature_wrapped = str_wrap(feature_label, width = 34),
    type = factor(type, levels = feature_type_levels)
  ) %>%
  group_by(type) %>%
  arrange(topk_freq, .by_group = TRUE) %>%
  mutate(feature_wrapped = factor(feature_wrapped, levels = unique(feature_wrapped))) %>%
  ungroup()

p_rf_stable_main <- ggplot(
  main_rf_plot_df,
  aes(x = topk_freq, y = feature_wrapped, fill = type)
) +
  geom_col(width = 0.70, color = "black", linewidth = 0.25) +
  facet_wrap(~ type, scales = "free_y", ncol = 1, drop = TRUE) +
  scale_fill_manual(values = col_type, name = "Feature class", drop = FALSE) +
  labs(
    x = sprintf("Selection frequency\n(top-%d across %d RF refits)", top_k, n_runs),
    y = "Predictive feature"
  ) +
  theme_clean(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 13.5),
    axis.text.y = element_text(size = 10.2, face = "bold", lineheight = 0.90),
    axis.text.x = element_text(size = 12.0, face = "bold"),
    axis.title.x = element_text(size = 13.5),
    axis.title.y = element_text(size = 13.5),
    panel.grid.major.y = element_blank(),
    panel.spacing.y = unit(1.0, "lines"),
    plot.margin = margin(8, 8, 8, 8)
  )

# ============================================================
# 18. MAIN FIGURE 6D — TOP STABLE FEATURE HEATMAP
# ============================================================

readable_stability <- stability_df %>%
  filter(type %in% feature_type_levels) %>%
  mutate(
    feature_label = humanise_feature_label(feature_label, type = type, max_chars = Inf),
    readable_label = !(type == "Reaction" & is_vmh_like_label(feature_label)),
    priority = case_when(
      type != "Reaction" ~ 2,
      readable_label ~ 1,
      TRUE ~ 0
    )
  ) %>%
  arrange(desc(priority), desc(topk_freq), desc(mean_importance))

top_heatmap_features <- readable_stability %>%
  slice_head(n = top_rf_features_heatmap) %>%
  pull(model_feature)

top_heatmap_features <- intersect(top_heatmap_features, colnames(X))

message("\nTop heatmap features retained: ", length(top_heatmap_features))
print(top_heatmap_features)

if (length(top_heatmap_features) >= 2) {
  
  X_heat <- X[, top_heatmap_features, drop = FALSE]
  
  cor_mat <- suppressWarnings(
    cor(
      X_heat,
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  )
  
  rownames(cor_mat) <- top_heatmap_features
  colnames(cor_mat) <- top_heatmap_features
  
  heatmap_labels <- rf_lookup %>%
    filter(model_feature %in% top_heatmap_features) %>%
    mutate(
      feature_label_full = humanise_feature_label(
        feature_label,
        type = type,
        max_chars = Inf
      ),
      feature_label_wrapped_raw = str_wrap(feature_label_full, width = 22),
      type = factor(type, levels = feature_type_levels),
      model_feature = factor(model_feature, levels = top_heatmap_features)
    ) %>%
    arrange(model_feature) %>%
    mutate(
      model_feature = as.character(model_feature),
      feature_label_wrapped = make.unique(
        as.character(feature_label_wrapped_raw),
        sep = " "
      )
    )
  
  if (nrow(heatmap_labels) != length(top_heatmap_features)) {
    warning(
      "Some top_heatmap_features were not found in rf_lookup. ",
      "Proceeding with matched features only."
    )
    
    top_heatmap_features <- intersect(
      top_heatmap_features,
      heatmap_labels$model_feature
    )
    
    cor_mat <- cor_mat[
      top_heatmap_features,
      top_heatmap_features,
      drop = FALSE
    ]
  }
  
  heatmap_order <- heatmap_labels$feature_label_wrapped[
    match(top_heatmap_features, heatmap_labels$model_feature)
  ]
  
  heatmap_df_raw <- as.data.frame(
    as.table(cor_mat),
    stringsAsFactors = FALSE
  )
  
  names(heatmap_df_raw)[seq_len(3)] <- c(
    "Feature1",
    "Feature2",
    "Spearman_rho"
  )
  
  heatmap_df <- heatmap_df_raw %>%
    mutate(
      Feature1 = as.character(Feature1),
      Feature2 = as.character(Feature2),
      Spearman_rho = as.numeric(Spearman_rho)
    ) %>%
    left_join(
      heatmap_labels %>%
        select(model_feature, Label1 = feature_label_wrapped),
      by = c("Feature1" = "model_feature")
    ) %>%
    left_join(
      heatmap_labels %>%
        select(model_feature, Label2 = feature_label_wrapped),
      by = c("Feature2" = "model_feature")
    ) %>%
    filter(!is.na(Label1), !is.na(Label2)) %>%
    mutate(
      Label1 = factor(Label1, levels = heatmap_order),
      Label2 = factor(Label2, levels = rev(heatmap_order)),
      x_num = as.numeric(Label1),
      y_num = as.numeric(Label2)
    ) %>%
    filter(!is.na(x_num), !is.na(y_num))
  
  label_position_df <- heatmap_labels %>%
    mutate(
      Label = factor(feature_label_wrapped, levels = heatmap_order),
      x_num = as.numeric(Label),
      y_num = length(heatmap_order) - x_num + 1
    ) %>%
    filter(!is.na(x_num), !is.na(y_num))
  
  row_type_boxes <- label_position_df %>%
    transmute(
      x = 0.25,
      y = y_num,
      type = type
    )
  
  col_type_boxes <- label_position_df %>%
    transmute(
      x = x_num,
      y = length(heatmap_order) + 0.75,
      type = type
    )
  
  write_csv(heatmap_df, csv_rf_feature_heatmap)
  
  p_rf_heatmap_main <- ggplot() +
    geom_tile(
      data = heatmap_df,
      aes(
        x = x_num,
        y = y_num,
        fill = Spearman_rho
      ),
      color = "white",
      linewidth = 0.25
    ) +
    geom_tile(
      data = row_type_boxes,
      aes(
        x = x,
        y = y,
        color = type
      ),
      fill = "white",
      width = 0.34,
      height = 0.78,
      linewidth = 1.35,
      show.legend = FALSE
    ) +
    geom_tile(
      data = col_type_boxes,
      aes(
        x = x,
        y = y,
        color = type
      ),
      fill = "white",
      width = 0.78,
      height = 0.34,
      linewidth = 1.35,
      show.legend = FALSE
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Spearman ρ"
    ) +
    scale_color_manual(
      values = col_type,
      guide = "none"
    ) +
    scale_x_continuous(
      breaks = seq_along(heatmap_order),
      labels = heatmap_order,
      expand = expansion(mult = c(0.08, 0.03))
    ) +
    scale_y_continuous(
      breaks = seq_along(rev(heatmap_order)),
      labels = rev(heatmap_order),
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Predictive feature",
      y = "Predictive feature"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(
        angle = 48,
        hjust = 1,
        vjust = 1,
        face = "bold",
        size = 9.2,
        lineheight = 0.86
      ),
      axis.text.y = element_text(
        face = "bold",
        size = 9.2,
        lineheight = 0.86
      ),
      axis.title = element_text(
        face = "bold",
        size = 14
      ),
      legend.title = element_text(
        face = "bold",
        size = 12
      ),
      legend.text = element_text(
        size = 11
      ),
      panel.grid = element_blank(),
      plot.margin = margin(12, 12, 12, 16)
    )
  
} else {
  
  p_rf_heatmap_main <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "Too few stable features for heatmap",
      fontface = "bold"
    ) +
    theme_void()
}

# ============================================================
# 19. CHECK ONLY PLOTTED FIGURE 6C/D LABELS FOR TRANSLATION
# ============================================================

plotted_rf_labels_for_translation <- bind_rows(
  main_rf_plot_df %>%
    mutate(source = "Figure6C") %>%
    select(source, type, model_feature, feature_key, feature_label),
  readable_stability %>%
    filter(model_feature %in% top_heatmap_features) %>%
    mutate(source = "Figure6D") %>%
    select(source, type, model_feature, feature_key, feature_label)
) %>%
  distinct(source, type, model_feature, feature_key, feature_label)

plotted_untranslated <- get_plotted_reaction_labels_needing_translation(
  plotted_rf_labels_for_translation,
  out_dir = out_dir
)

cat("\n============================================================\n")
cat("PLOTTED reaction labels still needing translation\n")
cat("============================================================\n")
print(plotted_untranslated, n = Inf)

# ============================================================
# 20. SAVE MAIN FIGURE 6 A-D
# ============================================================

figure6_final <- (p_summary | p_network) / (p_rf_stable_main | p_rf_heatmap_main) +
  plot_layout(
    widths = c(1.18, 1.22),
    heights = c(1.08, 1.18),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 30)
    )
  ) &
  theme(legend.position = "right")

ggsave(
  out_fig6_png,
  plot = figure6_final,
  width = 22,
  height = 19.2,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  out_fig6_pdf,
  plot = figure6_final,
  width = 22,
  height = 19.2,
  bg = "white",
  limitsize = FALSE
)

safe_print(figure6_final)

message("\nSaved main Figure 6:")
message("  ", out_fig6_png)
message("  ", out_fig6_pdf)

# ============================================================
# 21. SUPPLEMENTARY FIGURE 1 — RF QC
# ============================================================

rmse_med <- median(perf_df$rmse, na.rm = TRUE)
r2_med <- median(perf_df$r2, na.rm = TRUE)

resid_mean_med <- median(resid_stats_df$resid_mean, na.rm = TRUE)
resid_sd_med <- median(resid_stats_df$resid_sd, na.rm = TRUE)
shap_med <- median(resid_stats_df$shapiro_p, na.rm = TRUE)

std_resid_cor <- suppressWarnings(
  cor.test(
    test_all_df$obs,
    test_all_df$std_resid,
    method = "spearman",
    exact = FALSE
  )
)

std_resid_label <- paste0(
  "Spearman ρ = ", round(as.numeric(std_resid_cor$estimate), 2),
  "\np = ", signif(std_resid_cor$p.value, 2)
)

p_qc1 <- ggplot(test_all_df, aes(x = obs, y = pred)) +
  geom_point(alpha = 0.12, size = 1, colour = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_smooth(method = "lm", se = FALSE, colour = "red") +
  theme_bw(base_size = 12) +
  labs(
    title = "Observed vs predicted",
    x = "Observed PFAS elimination",
    y = "Predicted PFAS elimination"
  ) +
  annotate(
    "text",
    x = min(test_all_df$obs, na.rm = TRUE),
    y = max(test_all_df$pred, na.rm = TRUE),
    hjust = 0,
    vjust = 1.2,
    size = 4,
    label = paste0(
      "Median across runs:\n",
      "RMSE = ", round(rmse_med, 3), "\n",
      "R² = ", round(r2_med, 3)
    )
  )

p_qc2 <- ggplot(test_all_df, aes(x = pred, y = resid)) +
  geom_point(alpha = 0.12, size = 1, colour = "darkorange") +
  geom_smooth(method = "loess", se = FALSE, colour = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  theme_bw(base_size = 12) +
  labs(
    title = "Residuals vs predicted",
    x = "Predicted PFAS elimination",
    y = "Residual (observed − predicted)"
  ) +
  annotate(
    "text",
    x = min(test_all_df$pred, na.rm = TRUE),
    y = max(test_all_df$resid, na.rm = TRUE),
    hjust = 0,
    vjust = 1.2,
    size = 4,
    label = paste0(
      "Median across runs:\n",
      "Mean(resid) = ", round(resid_mean_med, 3), "\n",
      "SD(resid) = ", round(resid_sd_med, 3)
    )
  )

p_qc3 <- ggplot(test_all_df, aes(x = obs, y = std_resid)) +
  geom_point(alpha = 0.12, size = 1, colour = "purple3") +
  geom_smooth(method = "loess", se = FALSE, colour = "red") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  annotate(
    "text",
    x = min(test_all_df$obs, na.rm = TRUE),
    y = max(test_all_df$std_resid, na.rm = TRUE),
    hjust = 0,
    vjust = 1.2,
    size = 4,
    label = std_resid_label
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Standardized residuals vs observed",
    x = "Observed PFAS elimination",
    y = "Standardized residual"
  )

p_qc4 <- ggplot(test_all_df, aes(x = resid)) +
  geom_density(fill = "grey70", alpha = 0.7, colour = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
  theme_bw(base_size = 12) +
  labs(
    title = "Residual density",
    x = "Residual",
    y = "Density"
  ) +
  annotate(
    "text",
    x = min(test_all_df$resid, na.rm = TRUE),
    y = Inf,
    hjust = 0,
    vjust = 1.2,
    size = 4,
    label = paste0("Median Shapiro p = ", signif(shap_med, 3))
  )

qc_4panel <- (p_qc1 | p_qc2) / (p_qc3 | p_qc4) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 20)
    )
  )

ggsave(
  fig_qc_4panel,
  plot = qc_4panel,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  fig_qc_4panel_pdf,
  plot = qc_4panel,
  width = 12,
  height = 10,
  bg = "white",
  limitsize = FALSE
)

safe_print(qc_4panel)

# ============================================================
# 22. SUPPLEMENTARY FIGURE 2 — RF STABILITY DIAGNOSTICS
# ============================================================

rmse_iqr <- quantile(perf_df$rmse, c(0.25, 0.75), na.rm = TRUE)
r2_iqr <- quantile(perf_df$r2, c(0.25, 0.75), na.rm = TRUE)

n_selected_once <- sum(stability_df$topk_count > 0, na.rm = TRUE)
n_selected_10 <- sum(stability_df$topk_freq >= 0.10, na.rm = TRUE)
n_selected_25 <- sum(stability_df$topk_freq >= 0.25, na.rm = TRUE)

p_stab1 <- ggplot(perf_df, aes(x = "", y = rmse)) +
  geom_boxplot(outlier.alpha = 0.3, fill = "grey85", colour = "black") +
  geom_jitter(width = 0.08, alpha = 0.35, size = 1.2) +
  annotate(
    "text",
    x = 1.28,
    y = max(perf_df$rmse, na.rm = TRUE),
    hjust = 0,
    vjust = 1,
    size = 4,
    fontface = "bold",
    label = paste0(
      "Median = ", round(rmse_med, 3),
      "\nIQR = ", round(rmse_iqr[1], 3), "–", round(rmse_iqr[2], 3),
      "\nRuns = ", n_runs
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "RMSE across refits",
    x = NULL,
    y = "RMSE"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_stab2 <- ggplot(perf_df, aes(x = "", y = r2)) +
  geom_boxplot(outlier.alpha = 0.3, fill = "grey85", colour = "black") +
  geom_jitter(width = 0.08, alpha = 0.35, size = 1.2) +
  annotate(
    "text",
    x = 1.28,
    y = max(perf_df$r2, na.rm = TRUE),
    hjust = 0,
    vjust = 1,
    size = 4,
    fontface = "bold",
    label = paste0(
      "Median = ", round(r2_med, 3),
      "\nIQR = ", round(r2_iqr[1], 3), "–", round(r2_iqr[2], 3),
      "\nRuns = ", n_runs
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "R² across refits",
    x = NULL,
    y = "R²"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_stab3 <- ggplot(stability_df, aes(x = topk_freq)) +
  geom_histogram(bins = 25, color = "black", fill = "grey80") +
  geom_vline(xintercept = 0.10, linetype = "dashed", colour = "grey35") +
  geom_vline(xintercept = 0.25, linetype = "dashed", colour = "red") +
  annotate(
    "text",
    x = 0.25,
    y = Inf,
    label = "25%",
    hjust = -0.1,
    vjust = 1.4,
    size = 3.5,
    fontface = "bold",
    colour = "red"
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    hjust = 1.05,
    vjust = 1.25,
    size = 4,
    fontface = "bold",
    label = paste0(
      "Selected ≥1×: ", n_selected_once,
      "\n≥10% runs: ", n_selected_10,
      "\n≥25% runs: ", n_selected_25
    )
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = sprintf("Selection frequency distribution (top-%d)", top_k),
    x = "Selection frequency",
    y = "Number of features"
  )

type_composition <- stability_df %>%
  mutate(
    stability_group = case_when(
      topk_freq >= 0.25 ~ "≥25% of refits",
      topk_freq >= 0.10 ~ "10–24% of refits",
      topk_freq > 0 ~ "Selected at least once",
      TRUE ~ "Never selected"
    ),
    stability_group = factor(
      stability_group,
      levels = c("≥25% of refits", "10–24% of refits", "Selected at least once", "Never selected")
    )
  ) %>%
  filter(stability_group != "Never selected") %>%
  dplyr::count(stability_group, type, name = "n") %>%
  group_by(stability_group) %>%
  mutate(
    total_in_group = sum(n),
    prop = n / total_in_group,
    label = ifelse(n > 0, paste0(n, "\n", round(100 * prop), "%"), "")
  ) %>%
  ungroup()

p_stab4 <- ggplot(
  type_composition,
  aes(x = stability_group, y = n, fill = type)
) +
  geom_col(color = "black", linewidth = 0.25, width = 0.72) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.3,
    fontface = "bold"
  ) +
  scale_fill_manual(values = col_type, name = "Feature class", drop = FALSE) +
  labs(
    title = "Feature classes among selected predictors",
    x = NULL,
    y = "Number of features"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

stability_4panel <- (p_stab1 | p_stab2) / (p_stab3 | p_stab4) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 20)
    )
  )

ggsave(
  fig_stability_4panel,
  plot = stability_4panel,
  width = 13,
  height = 10.5,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  fig_stability_4panel_pdf,
  plot = stability_4panel,
  width = 13,
  height = 10.5,
  bg = "white",
  limitsize = FALSE
)

safe_print(stability_4panel)

# ============================================================
# 23. ADDITIONAL SUPPLEMENTARY FIGURE — TOP-STABLE FEATURES
# ============================================================

extra_top_plot_df <- stability_df %>%
  filter(type %in% feature_type_levels) %>%
  group_by(type) %>%
  slice_max(order_by = topk_freq, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    feature_wrapped = str_wrap(humanise_feature_label(feature_label, type = type, max_chars = Inf), width = 34),
    type = factor(type, levels = feature_type_levels)
  ) %>%
  group_by(type) %>%
  arrange(topk_freq, .by_group = TRUE) %>%
  mutate(feature_wrapped = factor(feature_wrapped, levels = unique(feature_wrapped))) %>%
  ungroup()

p_extra_top <- ggplot(
  extra_top_plot_df,
  aes(x = topk_freq, y = feature_wrapped, fill = type)
) +
  geom_col(width = 0.78, color = "black", linewidth = 0.25) +
  facet_wrap(~ type, scales = "free_y", ncol = 2, drop = TRUE) +
  scale_fill_manual(values = col_type, name = "Feature class", drop = FALSE) +
  labs(
    title = "Top stable RF predictors by class",
    x = sprintf("Selection frequency (top-%d across %d RF refits)", top_k, n_runs),
    y = "Predictive feature"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10.5),
    axis.text.y = element_text(size = 8.5, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

p_extra_heat <- p_rf_heatmap_main +
  labs(title = "Spearman correlations among top stable RF predictors") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))

extra_rf_figure <- p_extra_top / p_extra_heat +
  plot_layout(heights = c(1.05, 1.00)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 20)
    )
  )

ggsave(
  fig_rf_extra,
  plot = extra_rf_figure,
  width = 14.5,
  height = 13,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  fig_rf_extra_pdf,
  plot = extra_rf_figure,
  width = 14.5,
  height = 13,
  bg = "white",
  limitsize = FALSE
)

safe_print(extra_rf_figure)

# ============================================================
# 24. FINAL MESSAGES
# ============================================================

message("\n============================================================")
message("DONE")
message("Outputs written to:")
message(out_dir)
message("============================================================")

message("\nMain Figure 6:")
message(" - ", basename(out_fig6_png))
message(" - ", basename(out_fig6_pdf))
message(" - ", basename(out_csv_all))
message(" - ", basename(out_csv_sig))
message(" - ", basename(out_csv_nodes))
message(" - ", basename(out_csv_edges))

message("\nMetabolite translation diagnostics:")
message(" - ", basename(csv_metabolite_untranslated))
message(" - ", basename(txt_metabolite_template))

message("\nRF supplementary figures:")
message(" - ", basename(fig_qc_4panel))
message(" - ", basename(fig_qc_4panel_pdf))
message(" - ", basename(fig_stability_4panel))
message(" - ", basename(fig_stability_4panel_pdf))
message(" - ", basename(fig_rf_extra))
message(" - ", basename(fig_rf_extra_pdf))

message("\nRF CSVs:")
message(" - ", basename(csv_perf))
message(" - ", basename(csv_resid_stats))
message(" - ", basename(csv_importance_all))
message(" - ", basename(csv_selected_topk))
message(" - ", basename(csv_stability))
message(" - ", basename(csv_topstable_type))
message(" - ", basename(csv_topdisplay))
message(" - ", basename(csv_stablecand))
message(" - ", basename(csv_stab_overall))
message(" - ", basename(csv_stab_bytype))
message(" - ", basename(csv_rf_filtering))
message(" - ", basename(csv_rf_feature_heatmap))

message("\nRF refits completed: ", n_runs)


##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
############################################################
# Diagnostic analysis for apparent PFAS-exposed subclusters
# in Figure 4A subsystem heatmap
#
# Fixes included:
# - Duplicate taxonomy row names handled with make.unique()
# - metadata_exp already contains Sample_Name, so rownames_to_column()
#   is NOT used when joining metadata to cluster assignments
# - Duplicate metadata column names cleaned with make.unique()
#
# This script:
# 1. Recreates exposed-sample subsystem clustering
# 2. Assigns exposed samples to two subsystem-derived clusters
# 3. Tests whether clusters associate with metadata variables
# 4. Identifies subsystem/reaction features driving the split
# 5. Checks whether the same split is visible at the species level
# 6. Saves diagnostic plots and CSV outputs
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(phyloseq)
  library(vegan)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# ============================================================
# 1. Paths
# ============================================================

subsystem_path <- file.path("data", "model_outputs", "SubsystemAbundance.csv")
reaction_path  <- file.path("data", "model_outputs", "ReactionAbundance.csv")
metadata_path  <- file.path("data", "metadata", "Sample_Metadata_New.tsv")

otu_file <- file.path("data", "processed", "otu_table.tsv")
tax_file <- file.path("data", "processed", "taxonomy_table.tsv")

output_dir <- file.path("results", "supplementary_figures", "cluster_diagnostics")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 2. Settings
# ============================================================

condition_of_interest <- "PFAS_Exposed"
k_clusters <- 2
top_n_subsystems <- 35
top_n_species <- 100
presence_threshold <- 1e-9

candidate_technical_vars <- c(
  "Sequencing_Run",
  "Run",
  "Lane",
  "Flowcell",
  "Library_Batch",
  "Extraction_Batch",
  "DNA_Extraction_Batch",
  "Sequencing_Batch",
  "Plate",
  "Date",
  "Read_Depth",
  "Total_Reads",
  "Mapped_Reads",
  "Host_Reads",
  "N50",
  "Assembly_N50",
  "Number_of_Contigs",
  "Contigs"
)

candidate_biological_vars <- c(
  "Age",
  "Sex",
  "BMI",
  "weight",
  "PFOS",
  "PFHxS",
  "PFOA",
  "PFNA",
  "sum_PFAS",
  "average_k",
  "Average_k",
  "k_PFOA",
  "k_PFPeS",
  "k_PFHxS",
  "k_PFHpS",
  "k_LPFOS",
  "k_PFOS_MP1",
  "k_PFOS_MP345",
  "k_PFOS_MP26"
)

# ============================================================
# 3. Helper functions
# ============================================================

to_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c(
    "", "NA", "NaN", "nd", "ND", "na", "n/a", "N/A",
    "<LOD", "<LOQ", "< LOQ", "LOD", "LOQ"
  )] <- NA
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- gsub("^<\\s*", "", x)
  suppressWarnings(as.numeric(x))
}

read_clean_delim <- function(path, sep = "\t") {
  df <- data.table::fread(
    path,
    sep = sep,
    data.table = FALSE,
    check.names = FALSE
  )
  colnames(df) <- make.unique(trimws(colnames(df)))
  df
}

find_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

zscore_rows <- function(mat) {
  mat_z <- t(scale(t(mat)))
  mat_z[is.na(mat_z)] <- 0
  mat_z[is.nan(mat_z)] <- 0
  mat_z[is.infinite(mat_z)] <- 0
  mat_z
}

matrix_from_feature_csv <- function(path, feature_col_name = NULL) {
  df <- data.table::fread(
    path,
    data.table = FALSE,
    check.names = FALSE
  )
  
  colnames(df) <- make.unique(trimws(colnames(df)))
  
  if (is.null(feature_col_name)) {
    feature_col_name <- colnames(df)[1]
  }
  
  feature_col_index <- which(colnames(df) == feature_col_name)[1]
  if (is.na(feature_col_index)) {
    feature_col_index <- 1
  }
  
  colnames(df)[feature_col_index] <- "Feature"
  
  df$Feature <- make.unique(as.character(df$Feature))
  
  numeric_cols <- setdiff(colnames(df), "Feature")
  df[numeric_cols] <- lapply(df[numeric_cols], to_numeric_safe)
  
  mat <- df %>%
    column_to_rownames("Feature") %>%
    as.matrix()
  
  storage.mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  
  mat
}

safe_wilcox <- function(df, value_col, group_col = "Cluster") {
  dat <- df %>%
    select(all_of(c(value_col, group_col))) %>%
    mutate(
      value = to_numeric_safe(.data[[value_col]]),
      group = as.factor(.data[[group_col]])
    ) %>%
    filter(!is.na(value), !is.na(group))
  
  if (n_distinct(dat$group) != 2 || nrow(dat) < 4) {
    return(tibble(
      variable = value_col,
      p_value = NA_real_,
      test = "wilcox"
    ))
  }
  
  res <- tryCatch(
    wilcox.test(value ~ group, data = dat),
    error = function(e) NULL
  )
  
  if (is.null(res)) {
    return(tibble(
      variable = value_col,
      p_value = NA_real_,
      test = "wilcox"
    ))
  }
  
  tibble(
    variable = value_col,
    p_value = res$p.value,
    test = "wilcox"
  )
}

safe_fisher_or_chisq <- function(df, value_col, group_col = "Cluster") {
  dat <- df %>%
    select(all_of(c(value_col, group_col))) %>%
    mutate(
      value = as.factor(.data[[value_col]]),
      group = as.factor(.data[[group_col]])
    ) %>%
    filter(!is.na(value), !is.na(group))
  
  if (n_distinct(dat$group) < 2 || n_distinct(dat$value) < 2) {
    return(tibble(
      variable = value_col,
      p_value = NA_real_,
      test = "categorical"
    ))
  }
  
  tab <- table(dat$group, dat$value)
  
  res <- tryCatch(
    {
      if (any(tab < 5)) fisher.test(tab) else chisq.test(tab)
    },
    error = function(e) NULL
  )
  
  if (is.null(res)) {
    return(tibble(
      variable = value_col,
      p_value = NA_real_,
      test = "categorical"
    ))
  }
  
  tibble(
    variable = value_col,
    p_value = res$p.value,
    test = ifelse(any(tab < 5), "fisher", "chisq")
  )
}

# ============================================================
# 4. Read metadata
# ============================================================

metadata <- read_clean_delim(metadata_path, sep = "\t")

sample_col <- find_first_existing_col(
  metadata,
  c("Sample_Name", "SampleID", "Sample_ID", "Sample")
)

condition_col <- find_first_existing_col(
  metadata,
  c("Condition", "Group")
)

if (is.null(sample_col) || is.null(condition_col)) {
  stop("Could not find sample or condition column in metadata.")
}

metadata <- metadata %>%
  mutate(
    Sample_Name = trimws(as.character(.data[[sample_col]])),
    Condition = trimws(as.character(.data[[condition_col]]))
  )

# Remove duplicate Sample_Name columns if make.unique created Sample_Name.1 etc.
# Keep the standardized Sample_Name created above.
metadata <- metadata %>%
  select(-matches("^Sample_Name\\.\\d+$"))

k_cols <- c(
  "k_PFOA",
  "k_PFPeS",
  "k_PFHxS",
  "k_PFHpS",
  "k_LPFOS",
  "k_PFOS_MP1",
  "k_PFOS_MP345",
  "k_PFOS_MP26"
)

available_k_cols <- intersect(k_cols, colnames(metadata))

if (length(available_k_cols) > 0) {
  metadata$Average_k <- rowMeans(
    as.data.frame(
      lapply(metadata[, available_k_cols, drop = FALSE], to_numeric_safe)
    ),
    na.rm = TRUE
  )
} else {
  metadata$Average_k <- NA_real_
}

cat("\nMetadata dimensions:\n")
print(dim(metadata))

cat("\nMetadata sample column used:\n")
print(sample_col)

cat("\nMetadata condition column used:\n")
print(condition_col)

# ============================================================
# 5. Read subsystem and reaction matrices
# ============================================================

subsystem_mat <- matrix_from_feature_csv(subsystem_path)
reaction_mat <- matrix_from_feature_csv(reaction_path)

cat("\nSubsystem matrix dimensions:\n")
print(dim(subsystem_mat))

cat("\nReaction matrix dimensions:\n")
print(dim(reaction_mat))

# ============================================================
# 6. Build species matrix from OTU + taxonomy files
# ============================================================

otu_df <- read.delim(
  otu_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

otu_mat <- as.matrix(otu_df)
storage.mode(otu_mat) <- "numeric"
otu_mat[is.na(otu_mat)] <- 0
rownames(otu_mat) <- trimws(rownames(otu_mat))

cat("\nOTU matrix dimensions:\n")
print(dim(otu_mat))

tax_df_raw <- read.delim(
  tax_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(tax_df_raw) <- make.unique(trimws(colnames(tax_df_raw)))

tax_id_col <- colnames(tax_df_raw)[1]
tax_ids <- trimws(as.character(tax_df_raw[[tax_id_col]]))
tax_ids_unique <- make.unique(tax_ids)

tax_df <- tax_df_raw
rownames(tax_df) <- tax_ids_unique
tax_df[[tax_id_col]] <- tax_ids_unique

cat("\nTaxonomy table dimensions:\n")
print(dim(tax_df))

cat("\nNumber of duplicate raw taxonomy identifiers detected:\n")
print(sum(duplicated(tax_ids)))

otu_taxa_raw <- rownames(otu_mat)
otu_taxa_unique <- make.unique(otu_taxa_raw)
rownames(otu_mat) <- otu_taxa_unique

common_taxa <- intersect(rownames(otu_mat), rownames(tax_df))

if (length(common_taxa) == 0) {
  stop(
    "No matching taxa between OTU table and taxonomy table after make.unique(). ",
    "Check whether the first taxonomy column matches OTU rownames."
  )
}

otu_mat <- otu_mat[common_taxa, , drop = FALSE]
tax_df <- tax_df[common_taxa, , drop = FALSE]

cat("\nCommon taxa between OTU and taxonomy:\n")
print(length(common_taxa))

species_col <- find_first_existing_col(
  tax_df,
  c("Species", "species", "s", "Taxon", "taxon")
)

if (is.null(species_col)) {
  warning("No Species column found in taxonomy table. Using matched taxon IDs.")
  species_labels <- rownames(otu_mat)
} else {
  species_labels <- as.character(tax_df[[species_col]])
  species_labels <- trimws(species_labels)
  species_labels[is.na(species_labels) | species_labels == ""] <- rownames(otu_mat)[
    is.na(species_labels) | species_labels == ""
  ]
}

species_df <- as.data.frame(otu_mat) %>%
  rownames_to_column("Original_Taxon") %>%
  mutate(Species = species_labels) %>%
  select(-Original_Taxon) %>%
  group_by(Species) %>%
  summarise(
    across(everything(), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  )

species_mat <- species_df %>%
  column_to_rownames("Species") %>%
  as.matrix()

storage.mode(species_mat) <- "numeric"
species_mat[is.na(species_mat)] <- 0

sample_totals <- colSums(species_mat, na.rm = TRUE)
sample_totals[sample_totals == 0] <- NA

species_rel <- sweep(species_mat, 2, sample_totals, FUN = "/")
species_rel[is.na(species_rel)] <- 0

cat("\nSpecies relative-abundance matrix dimensions:\n")
print(dim(species_rel))

# ============================================================
# 7. Align samples across subsystem, reaction, species, metadata
# ============================================================

common_samples <- Reduce(
  intersect,
  list(
    colnames(subsystem_mat),
    colnames(reaction_mat),
    colnames(species_rel),
    metadata$Sample_Name
  )
)

cat("\nCommon samples across subsystem, reaction, species, metadata:\n")
print(length(common_samples))

if (length(common_samples) < 10) {
  stop("Too few common samples across subsystem, reaction, species, and metadata.")
}

metadata_aligned <- metadata %>%
  filter(Sample_Name %in% common_samples) %>%
  distinct(Sample_Name, .keep_all = TRUE) %>%
  arrange(match(Sample_Name, common_samples))

rownames(metadata_aligned) <- metadata_aligned$Sample_Name

subsystem_mat <- subsystem_mat[, metadata_aligned$Sample_Name, drop = FALSE]
reaction_mat  <- reaction_mat[, metadata_aligned$Sample_Name, drop = FALSE]
species_rel   <- species_rel[, metadata_aligned$Sample_Name, drop = FALSE]

exposed_samples <- metadata_aligned %>%
  filter(Condition == condition_of_interest) %>%
  pull(Sample_Name)

subsystem_exp <- subsystem_mat[, exposed_samples, drop = FALSE]
reaction_exp  <- reaction_mat[, exposed_samples, drop = FALSE]
species_exp   <- species_rel[, exposed_samples, drop = FALSE]

metadata_exp <- metadata_aligned %>%
  filter(Sample_Name %in% exposed_samples) %>%
  arrange(match(Sample_Name, exposed_samples))

rownames(metadata_exp) <- metadata_exp$Sample_Name

cat("\nNumber of PFAS-exposed samples used:\n")
print(length(exposed_samples))

# ============================================================
# 8. Filter variable features and select top subsystem rows
# ============================================================

subsystem_exp <- subsystem_exp[
  apply(subsystem_exp, 1, function(x) sd(x, na.rm = TRUE) > 0),
  ,
  drop = FALSE
]

reaction_exp <- reaction_exp[
  apply(reaction_exp, 1, function(x) sd(x, na.rm = TRUE) > 0),
  ,
  drop = FALSE
]

species_exp <- species_exp[
  apply(species_exp, 1, function(x) sd(x, na.rm = TRUE) > 0),
  ,
  drop = FALSE
]

top_subsystems <- subsystem_exp %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  mutate(total_abundance = rowSums(across(-Feature), na.rm = TRUE)) %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = min(top_n_subsystems, nrow(.))) %>%
  pull(Feature)

subsystem_top <- subsystem_exp[top_subsystems, , drop = FALSE]
subsystem_top_z <- zscore_rows(subsystem_top)

# ============================================================
# 9. Derive two clusters from exposed subsystem heatmap
# ============================================================

sample_dist_subsystem <- dist(t(subsystem_top_z), method = "euclidean")
sample_hclust_subsystem <- hclust(sample_dist_subsystem, method = "ward.D2")

cluster_assignments <- cutree(sample_hclust_subsystem, k = k_clusters)

cluster_df <- tibble(
  Sample_Name = names(cluster_assignments),
  Cluster = paste0("Cluster_", cluster_assignments)
) %>%
  left_join(
    metadata_exp,
    by = "Sample_Name"
  )

write.csv(
  cluster_df,
  file.path(output_dir, "Exposed_SubsystemCluster_Assignments.csv"),
  row.names = FALSE
)

cat("\nCluster sizes:\n")
print(table(cluster_df$Cluster))

# ============================================================
# 10. Heatmap with subsystem-derived clusters
# ============================================================

annotation_col <- cluster_df %>%
  select(
    Sample_Name,
    Cluster,
    any_of(c("Average_k", "average_k", "Sex", "BMI", "Age"))
  ) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$Sample_Name
annotation_col$Sample_Name <- NULL

pheatmap(
  subsystem_top_z,
  cluster_rows = TRUE,
  cluster_cols = sample_hclust_subsystem,
  annotation_col = annotation_col,
  fontsize_row = 7,
  fontsize_col = 6,
  main = "PFAS-exposed samples: subsystem-derived clusters",
  filename = file.path(output_dir, "Diagnostic_Exposed_SubsystemHeatmap_WithClusters.png"),
  width = 12,
  height = 8
)

# ============================================================
# 11. Test cluster association with metadata variables
# ============================================================

available_technical_vars <- intersect(candidate_technical_vars, colnames(cluster_df))
available_biological_vars <- intersect(candidate_biological_vars, colnames(cluster_df))

all_candidate_vars <- unique(c(available_technical_vars, available_biological_vars))

metadata_tests <- purrr::map_dfr(
  all_candidate_vars,
  function(v) {
    x <- cluster_df[[v]]
    x_num <- to_numeric_safe(x)
    is_numeric_like <- sum(!is.na(x_num)) >= max(5, floor(0.5 * length(x)))
    
    if (is_numeric_like) {
      cluster_df[[v]] <- x_num
      safe_wilcox(cluster_df, v, group_col = "Cluster")
    } else {
      safe_fisher_or_chisq(cluster_df, v, group_col = "Cluster")
    }
  }
) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    variable_type = case_when(
      variable %in% available_technical_vars ~ "technical",
      variable %in% available_biological_vars ~ "biological",
      TRUE ~ "other"
    )
  ) %>%
  arrange(p_value)

write.csv(
  metadata_tests,
  file.path(output_dir, "Cluster_Association_MetadataTests.csv"),
  row.names = FALSE
)

cat("\nTop metadata associations with subsystem-derived clusters:\n")
print(metadata_tests %>% slice_head(n = 20))

# ============================================================
# 12. Plot top numeric metadata associations
# ============================================================

top_numeric_vars <- metadata_tests %>%
  filter(test == "wilcox", !is.na(p_value)) %>%
  slice_head(n = 8) %>%
  pull(variable)

if (length(top_numeric_vars) > 0) {
  
  numeric_plot_df <- cluster_df %>%
    select(Sample_Name, Cluster, all_of(top_numeric_vars)) %>%
    pivot_longer(
      cols = all_of(top_numeric_vars),
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    mutate(Value = to_numeric_safe(Value)) %>%
    filter(!is.na(Value))
  
  p_numeric <- ggplot(
    numeric_plot_df,
    aes(x = Cluster, y = Value, fill = Cluster)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 2.2, alpha = 0.85) +
    facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    ) +
    labs(x = NULL, y = "Value")
  
  ggsave(
    file.path(output_dir, "Cluster_Metadata_NumericVariables.png"),
    p_numeric,
    width = 12,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    file.path(output_dir, "Cluster_Metadata_NumericVariables.pdf"),
    p_numeric,
    width = 12,
    height = 7,
    bg = "white"
  )
}

# ============================================================
# 13. Identify subsystem features driving clusters
# ============================================================

subsystem_cluster_long <- subsystem_exp %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  pivot_longer(
    cols = -Feature,
    names_to = "Sample_Name",
    values_to = "Value"
  ) %>%
  left_join(cluster_df %>% select(Sample_Name, Cluster), by = "Sample_Name")

subsystem_cluster_tests <- subsystem_cluster_long %>%
  group_by(Feature) %>%
  summarise(
    mean_cluster_1 = mean(Value[Cluster == "Cluster_1"], na.rm = TRUE),
    mean_cluster_2 = mean(Value[Cluster == "Cluster_2"], na.rm = TRUE),
    median_cluster_1 = median(Value[Cluster == "Cluster_1"], na.rm = TRUE),
    median_cluster_2 = median(Value[Cluster == "Cluster_2"], na.rm = TRUE),
    log2FC_cluster2_vs_1 = log2((mean_cluster_2 + 1e-9) / (mean_cluster_1 + 1e-9)),
    p_value = tryCatch(
      wilcox.test(Value ~ Cluster)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    abs_log2FC = abs(log2FC_cluster2_vs_1)
  ) %>%
  arrange(p_value)

write.csv(
  subsystem_cluster_tests,
  file.path(output_dir, "Subsystem_Features_Driving_ExposedClusters.csv"),
  row.names = FALSE
)

cat("\nTop subsystem features driving exposed clusters:\n")
print(subsystem_cluster_tests %>% slice_head(n = 25))

# ============================================================
# 14. Identify reaction features driving clusters
# ============================================================

reaction_cluster_long <- reaction_exp %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  pivot_longer(
    cols = -Feature,
    names_to = "Sample_Name",
    values_to = "Value"
  ) %>%
  left_join(cluster_df %>% select(Sample_Name, Cluster), by = "Sample_Name")

reaction_cluster_tests <- reaction_cluster_long %>%
  group_by(Feature) %>%
  summarise(
    mean_cluster_1 = mean(Value[Cluster == "Cluster_1"], na.rm = TRUE),
    mean_cluster_2 = mean(Value[Cluster == "Cluster_2"], na.rm = TRUE),
    prevalence_cluster_1 = mean(Value[Cluster == "Cluster_1"] > presence_threshold, na.rm = TRUE),
    prevalence_cluster_2 = mean(Value[Cluster == "Cluster_2"] > presence_threshold, na.rm = TRUE),
    log2FC_cluster2_vs_1 = log2((mean_cluster_2 + 1e-9) / (mean_cluster_1 + 1e-9)),
    p_value = tryCatch(
      wilcox.test(Value ~ Cluster)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    abs_log2FC = abs(log2FC_cluster2_vs_1),
    prevalence_difference = prevalence_cluster_2 - prevalence_cluster_1,
    abs_prevalence_difference = abs(prevalence_difference)
  ) %>%
  arrange(p_value)

write.csv(
  reaction_cluster_tests,
  file.path(output_dir, "Reaction_Features_Driving_ExposedClusters.csv"),
  row.names = FALSE
)

cat("\nTop reaction features driving exposed clusters:\n")
print(reaction_cluster_tests %>% slice_head(n = 25))

# ============================================================
# 15. Identify species associated with subsystem-derived clusters
# ============================================================

species_cluster_long <- species_exp %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  pivot_longer(
    cols = -Feature,
    names_to = "Sample_Name",
    values_to = "Value"
  ) %>%
  left_join(cluster_df %>% select(Sample_Name, Cluster), by = "Sample_Name")

species_cluster_tests <- species_cluster_long %>%
  group_by(Feature) %>%
  summarise(
    prevalence_cluster_1 = mean(Value[Cluster == "Cluster_1"] > 0, na.rm = TRUE),
    prevalence_cluster_2 = mean(Value[Cluster == "Cluster_2"] > 0, na.rm = TRUE),
    mean_cluster_1 = mean(Value[Cluster == "Cluster_1"], na.rm = TRUE),
    mean_cluster_2 = mean(Value[Cluster == "Cluster_2"], na.rm = TRUE),
    log2FC_cluster2_vs_1 = log2((mean_cluster_2 + 1e-9) / (mean_cluster_1 + 1e-9)),
    p_value = tryCatch(
      wilcox.test(Value ~ Cluster)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    abs_log2FC = abs(log2FC_cluster2_vs_1),
    prevalence_difference = prevalence_cluster_2 - prevalence_cluster_1,
    abs_prevalence_difference = abs(prevalence_difference)
  ) %>%
  arrange(p_value)

write.csv(
  species_cluster_tests,
  file.path(output_dir, "Species_Features_Associated_With_SubsystemClusters.csv"),
  row.names = FALSE
)

cat("\nTop species associated with subsystem-derived clusters:\n")
print(species_cluster_tests %>% slice_head(n = 25))

# ============================================================
# 16. Compare subsystem and species sample structure
# ============================================================

top_species <- species_exp %>%
  as.data.frame() %>%
  rownames_to_column("Feature") %>%
  mutate(row_sd = apply(species_exp, 1, sd, na.rm = TRUE)) %>%
  arrange(desc(row_sd)) %>%
  slice_head(n = min(top_n_species, nrow(.))) %>%
  pull(Feature)

species_top <- species_exp[top_species, , drop = FALSE]
species_top_z <- zscore_rows(species_top)

sample_dist_species <- dist(t(species_top_z), method = "euclidean")

mantel_res <- vegan::mantel(
  as.dist(sample_dist_subsystem),
  as.dist(sample_dist_species),
  method = "spearman",
  permutations = 999
)

capture.output(
  mantel_res,
  file = file.path(output_dir, "Mantel_Subsystem_vs_Species_Distances.txt")
)

cat("\nMantel test: subsystem vs species sample distances\n")
print(mantel_res)

# ============================================================
# 17. PCA comparison: subsystem vs species
# ============================================================

subsystem_pca <- prcomp(t(subsystem_top_z), scale. = FALSE)
species_pca <- prcomp(t(species_top_z), scale. = FALSE)

subsystem_pca_df <- as.data.frame(subsystem_pca$x[, 1:2]) %>%
  rownames_to_column("Sample_Name") %>%
  left_join(cluster_df %>% select(Sample_Name, Cluster), by = "Sample_Name")

species_pca_df <- as.data.frame(species_pca$x[, 1:2]) %>%
  rownames_to_column("Sample_Name") %>%
  left_join(cluster_df %>% select(Sample_Name, Cluster), by = "Sample_Name")

subsystem_var <- summary(subsystem_pca)$importance[2, 1:2] * 100
species_var <- summary(species_pca)$importance[2, 1:2] * 100

p_subsystem_pca <- ggplot(
  subsystem_pca_df,
  aes(x = PC1, y = PC2, color = Cluster, label = Sample_Name)
) +
  geom_point(size = 3.5, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  theme_bw(base_size = 13) +
  labs(
    title = "Subsystem-level PCA",
    x = paste0("PC1 (", round(subsystem_var[1], 1), "%)"),
    y = paste0("PC2 (", round(subsystem_var[2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

p_species_pca <- ggplot(
  species_pca_df,
  aes(x = PC1, y = PC2, color = Cluster, label = Sample_Name)
) +
  geom_point(size = 3.5, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3, max.overlaps = 50) +
  theme_bw(base_size = 13) +
  labs(
    title = "Species-level PCA colored by subsystem cluster",
    x = paste0("PC1 (", round(species_var[1], 1), "%)"),
    y = paste0("PC2 (", round(species_var[2], 1), "%)")
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

p_pca_combined <- p_subsystem_pca | p_species_pca

ggsave(
  file.path(output_dir, "Subsystem_vs_Species_PCA_ClusterComparison.png"),
  p_pca_combined,
  width = 14,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(output_dir, "Subsystem_vs_Species_PCA_ClusterComparison.pdf"),
  p_pca_combined,
  width = 14,
  height = 6,
  bg = "white"
)

# ============================================================
# 18. Presence/absence diagnostics for possible stochastic features
# ============================================================

subsystem_presence_tests <- subsystem_cluster_long %>%
  mutate(Present = Value > presence_threshold) %>%
  group_by(Feature) %>%
  summarise(
    present_cluster_1 = mean(Present[Cluster == "Cluster_1"], na.rm = TRUE),
    present_cluster_2 = mean(Present[Cluster == "Cluster_2"], na.rm = TRUE),
    prevalence_difference = present_cluster_2 - present_cluster_1,
    p_value = tryCatch(
      fisher.test(table(Cluster, Present))$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    abs_prevalence_difference = abs(prevalence_difference)
  ) %>%
  arrange(desc(abs_prevalence_difference), p_value)

write.csv(
  subsystem_presence_tests,
  file.path(output_dir, "Subsystem_PresenceAbsence_Driving_Clusters.csv"),
  row.names = FALSE
)

reaction_presence_tests <- reaction_cluster_long %>%
  mutate(Present = Value > presence_threshold) %>%
  group_by(Feature) %>%
  summarise(
    present_cluster_1 = mean(Present[Cluster == "Cluster_1"], na.rm = TRUE),
    present_cluster_2 = mean(Present[Cluster == "Cluster_2"], na.rm = TRUE),
    prevalence_difference = present_cluster_2 - present_cluster_1,
    p_value = tryCatch(
      fisher.test(table(Cluster, Present))$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    abs_prevalence_difference = abs(prevalence_difference)
  ) %>%
  arrange(desc(abs_prevalence_difference), p_value)

write.csv(
  reaction_presence_tests,
  file.path(output_dir, "Reaction_PresenceAbsence_Driving_Clusters.csv"),
  row.names = FALSE
)

cat("\nTop subsystem presence/absence differences:\n")
print(subsystem_presence_tests %>% slice_head(n = 25))

cat("\nTop reaction presence/absence differences:\n")
print(reaction_presence_tests %>% slice_head(n = 25))

# ============================================================
# 19. Compact text report
# ============================================================

sink(file.path(output_dir, "Cluster_Diagnostic_Report.txt"))

cat("Diagnostic report for PFAS-exposed subsystem clustering\n")
cat("=====================================================\n\n")

cat("Paths used:\n")
cat("Subsystem: ", subsystem_path, "\n")
cat("Reaction:  ", reaction_path, "\n")
cat("Metadata:  ", metadata_path, "\n")
cat("OTU:       ", otu_file, "\n")
cat("Taxonomy:  ", tax_file, "\n\n")

cat("Number of exposed samples:\n")
print(length(exposed_samples))

cat("\nNumber of duplicate raw taxonomy identifiers detected:\n")
print(sum(duplicated(tax_ids)))

cat("\nNumber of common taxa after make.unique alignment:\n")
print(length(common_taxa))

cat("\nCluster sizes:\n")
print(table(cluster_df$Cluster))

cat("\nTop metadata associations:\n")
print(metadata_tests %>% slice_head(n = 20))

cat("\nTop subsystem feature differences:\n")
print(subsystem_cluster_tests %>% slice_head(n = 25))

cat("\nTop reaction feature differences:\n")
print(reaction_cluster_tests %>% slice_head(n = 25))

cat("\nTop subsystem presence/absence differences:\n")
print(subsystem_presence_tests %>% slice_head(n = 25))

cat("\nTop reaction presence/absence differences:\n")
print(reaction_presence_tests %>% slice_head(n = 25))

cat("\nTop species differences using subsystem-derived cluster labels:\n")
print(species_cluster_tests %>% slice_head(n = 25))

cat("\nMantel test: subsystem vs species sample distances:\n")
print(mantel_res)

sink()

cat("\nDone.\n")
cat("Diagnostic outputs saved to:\n", output_dir, "\n")