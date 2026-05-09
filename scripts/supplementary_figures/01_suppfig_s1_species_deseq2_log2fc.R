############################################################
# 01_suppfig_s1_species_deseq2_log2fc.R
#
# Supplementary Figure S1: species-level DESeq2 log2 fold-change plots.
#
############################################################

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
