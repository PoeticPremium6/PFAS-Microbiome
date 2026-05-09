# Figure 2
# A. Alpha diversity indices
# B. Bray-Curtis PCoA with average PFAS elimination rate
# C. Beta-dispersion distance to centroid
#
# Full working script:
# - Rebuilds phyloseq object
# - Removes samples with missing Condition
# - Avoids "NA" appearing on x-axis
# - Uses PFAS-exposed = purple, Reference = grey/black
# - Adds PERMANOVA and betadisper analyses
# - Saves individual panels and combined Figure 2
# ==========================================================

# -------------------
# 0. Load libraries
# -------------------
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(vegan)
library(cowplot)

# -------------------
# 1. File paths
# -------------------
otu_file  <- file.path("data", "processed", "otu_table.tsv")
tax_file  <- file.path("data", "processed", "taxonomy_table.tsv")
meta_file <- file.path("data", "metadata", "PFAS_Metadata.csv")

fig_dir <- file.path("results", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# -------------------
# 2. Import OTU table
# -------------------
otu_df <- read.delim(
  otu_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if ("Taxa" %in% colnames(otu_df)) {
  otu_df$Taxa <- trimws(otu_df$Taxa)
  rownames(otu_df) <- make.unique(otu_df$Taxa)
  otu_df$Taxa <- NULL
} else {
  rownames(otu_df) <- make.unique(trimws(rownames(otu_df)))
}

otu_mat <- as.matrix(otu_df)
otu <- otu_table(otu_mat, taxa_are_rows = TRUE)

# -------------------
# 3. Import taxonomy table
# -------------------
tax_df <- read.delim(
  tax_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if ("Taxa" %in% colnames(tax_df)) {
  tax_df$Taxa <- trimws(tax_df$Taxa)
  rownames(tax_df) <- make.unique(tax_df$Taxa)
  tax_df$Taxa <- NULL
} else {
  rownames(tax_df) <- make.unique(trimws(tax_df[, 1]))
  tax_df <- tax_df[, -1, drop = FALSE]
}

tax <- tax_table(as.matrix(tax_df))

# -------------------
# 4. Import metadata
# -------------------
meta_df <- read.csv(
  meta_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if ("Sample_Name" %in% colnames(meta_df)) {
  meta_df$Sample_Name <- trimws(meta_df$Sample_Name)
  rownames(meta_df) <- meta_df$Sample_Name
  meta_df$Sample_Name <- NULL
} else if ("SampleID" %in% colnames(meta_df)) {
  meta_df$SampleID <- trimws(meta_df$SampleID)
  rownames(meta_df) <- meta_df$SampleID
  meta_df$SampleID <- NULL
} else if ("Sample" %in% colnames(meta_df)) {
  meta_df$Sample <- trimws(meta_df$Sample)
  rownames(meta_df) <- meta_df$Sample
  meta_df$Sample <- NULL
} else {
  rownames(meta_df) <- trimws(rownames(meta_df))
}

samp <- sample_data(meta_df)

# -------------------
# 5. Match taxa and samples
# -------------------
common_taxa <- intersect(taxa_names(otu), rownames(tax_df))
otu <- prune_taxa(common_taxa, otu)
tax <- tax[common_taxa, ]

common_samples <- intersect(sample_names(otu), rownames(meta_df))
otu <- prune_samples(common_samples, otu)
samp <- samp[common_samples, ]

# -------------------
# 6. Build phyloseq object
# -------------------
physeq <- phyloseq(otu, tax, samp)

cat("\nInitial phyloseq object:\n")
print(physeq)

# -------------------
# 7. Subset to PFAS_Exposed and Reference only
# -------------------
physeq_subset <- subset_samples(
  physeq,
  !is.na(Condition) & Condition %in% c("PFAS_Exposed", "Reference")
)

physeq_subset <- prune_taxa(
  taxa_sums(physeq_subset) > 0,
  physeq_subset
)

physeq_subset <- prune_samples(
  sample_sums(physeq_subset) > 0,
  physeq_subset
)

sample_data(physeq_subset)$Condition <- droplevels(
  factor(
    sample_data(physeq_subset)$Condition,
    levels = c("PFAS_Exposed", "Reference")
  )
)

cat("\nSamples by condition after subsetting:\n")
print(table(sample_data(physeq_subset)$Condition, useNA = "ifany"))

# -------------------
# 8. Shared colours and labels
# -------------------
condition_fill <- c(
  "PFAS_Exposed" = "#8E63D7",
  "Reference" = "#8A8A8A"
)

condition_outline <- c(
  "PFAS_Exposed" = "#7B1FA2",
  "Reference" = "black"
)

condition_labels <- c(
  "PFAS_Exposed" = "PFAS-exposed",
  "Reference" = "Reference"
)

# -------------------
# 9. Helper function for panel labels
# -------------------
add_panel_label <- function(plot_obj, label_text, label_size = 24) {
  cowplot::ggdraw(plot_obj) +
    cowplot::draw_label(
      label_text,
      x = 0.015,
      y = 0.985,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = label_size
    )
}

# ==========================================================
# Figure 2
# A. Alpha diversity indices
# B. Bray-Curtis PCoA with average PFAS elimination rate
#
# Updated version:
# - Only includes Figure 2A and Figure 2B
# - Removes beta-dispersion panel to avoid redundancy with Figure 1D
# - Removes samples with missing Condition
# - Avoids "NA" appearing on x-axis
# - Uses PFAS-exposed = purple, Reference = grey/black
# - Saves individual panels and combined Figure 2AB
# ==========================================================

# -------------------
# 0. Load libraries
# -------------------
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(vegan)
library(cowplot)

# -------------------
# 1. File paths
# -------------------
otu_file  <- file.path("data", "processed", "otu_table.tsv")
tax_file  <- file.path("data", "processed", "taxonomy_table.tsv")
meta_file <- file.path("data", "metadata", "PFAS_Metadata.csv")

fig_dir <- file.path("results", "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# -------------------
# 2. Import OTU table
# -------------------
otu_df <- read.delim(
  otu_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if ("Taxa" %in% colnames(otu_df)) {
  otu_df$Taxa <- trimws(otu_df$Taxa)
  rownames(otu_df) <- make.unique(otu_df$Taxa)
  otu_df$Taxa <- NULL
} else {
  rownames(otu_df) <- make.unique(trimws(rownames(otu_df)))
}

otu_mat <- as.matrix(otu_df)
otu <- otu_table(otu_mat, taxa_are_rows = TRUE)

# -------------------
# 3. Import taxonomy table
# -------------------
tax_df <- read.delim(
  tax_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if ("Taxa" %in% colnames(tax_df)) {
  tax_df$Taxa <- trimws(tax_df$Taxa)
  rownames(tax_df) <- make.unique(tax_df$Taxa)
  tax_df$Taxa <- NULL
} else {
  rownames(tax_df) <- make.unique(trimws(tax_df[, 1]))
  tax_df <- tax_df[, -1, drop = FALSE]
}

tax <- tax_table(as.matrix(tax_df))

# -------------------
# 4. Import metadata
# -------------------
meta_df <- read.csv(
  meta_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if ("Sample_Name" %in% colnames(meta_df)) {
  meta_df$Sample_Name <- trimws(meta_df$Sample_Name)
  rownames(meta_df) <- meta_df$Sample_Name
  meta_df$Sample_Name <- NULL
} else if ("SampleID" %in% colnames(meta_df)) {
  meta_df$SampleID <- trimws(meta_df$SampleID)
  rownames(meta_df) <- meta_df$SampleID
  meta_df$SampleID <- NULL
} else if ("Sample" %in% colnames(meta_df)) {
  meta_df$Sample <- trimws(meta_df$Sample)
  rownames(meta_df) <- meta_df$Sample
  meta_df$Sample <- NULL
} else {
  rownames(meta_df) <- trimws(rownames(meta_df))
}

samp <- sample_data(meta_df)

# -------------------
# 5. Match taxa and samples
# -------------------
common_taxa <- intersect(taxa_names(otu), rownames(tax_df))
otu <- prune_taxa(common_taxa, otu)
tax <- tax[common_taxa, ]

common_samples <- intersect(sample_names(otu), rownames(meta_df))
otu <- prune_samples(common_samples, otu)
samp <- samp[common_samples, ]

# -------------------
# 6. Build phyloseq object
# -------------------
physeq <- phyloseq(otu, tax, samp)

cat("\nInitial phyloseq object:\n")
print(physeq)

# -------------------
# 7. Subset to PFAS_Exposed and Reference only
# -------------------
physeq_subset <- subset_samples(
  physeq,
  !is.na(Condition) & Condition %in% c("PFAS_Exposed", "Reference")
)

physeq_subset <- prune_taxa(
  taxa_sums(physeq_subset) > 0,
  physeq_subset
)

physeq_subset <- prune_samples(
  sample_sums(physeq_subset) > 0,
  physeq_subset
)

sample_data(physeq_subset)$Condition <- droplevels(
  factor(
    sample_data(physeq_subset)$Condition,
    levels = c("PFAS_Exposed", "Reference")
  )
)

cat("\nSamples by condition after subsetting:\n")
print(table(sample_data(physeq_subset)$Condition, useNA = "ifany"))

# -------------------
# 8. Shared colours and labels
# -------------------
condition_fill <- c(
  "PFAS_Exposed" = "#8E63D7",
  "Reference" = "#8A8A8A"
)

condition_outline <- c(
  "PFAS_Exposed" = "#7B1FA2",
  "Reference" = "black"
)

condition_labels <- c(
  "PFAS_Exposed" = "PFAS-exposed",
  "Reference" = "Reference"
)

# -------------------
# 9. Helper function for panel labels
# -------------------
add_panel_label <- function(plot_obj, label_text, label_size = 24) {
  cowplot::ggdraw(plot_obj) +
    cowplot::draw_label(
      label_text,
      x = 0.015,
      y = 0.985,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = label_size
    )
}

# ==========================================================
# Figure 2A: Alpha diversity
# ==========================================================

# -------------------
# 10. Calculate alpha diversity
# -------------------
alpha_div <- estimate_richness(
  physeq_subset,
  measures = c("Shannon", "Simpson", "Observed", "Chao1")
)

alpha_div$SampleID <- rownames(alpha_div)

alpha_meta <- as.data.frame(sample_data(physeq_subset)) %>%
  rownames_to_column("SampleID")

alpha_div_meta <- alpha_div %>%
  left_join(alpha_meta, by = "SampleID") %>%
  filter(!is.na(Condition)) %>%
  mutate(
    Condition = factor(Condition, levels = c("PFAS_Exposed", "Reference")),
    Condition_label = recode(
      as.character(Condition),
      "PFAS_Exposed" = "PFAS-exposed",
      "Reference" = "Reference"
    ),
    Condition_label = factor(Condition_label, levels = c("PFAS-exposed", "Reference"))
  )

alpha_long <- alpha_div_meta %>%
  pivot_longer(
    cols = c("Shannon", "Simpson", "Observed", "Chao1"),
    names_to = "Index",
    values_to = "Value"
  ) %>%
  filter(!is.na(Value), !is.na(Condition_label)) %>%
  mutate(
    Index = factor(Index, levels = c("Shannon", "Simpson", "Observed", "Chao1"))
  )

# -------------------
# 11. Compute t-tests
# -------------------
alpha_stats <- alpha_long %>%
  group_by(Index) %>%
  summarise(
    p = t.test(
      Value[Condition == "PFAS_Exposed"],
      Value[Condition == "Reference"]
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = paste0("p = ", signif(p, 3))
  )

print(alpha_stats)

alpha_long <- alpha_long %>%
  left_join(alpha_stats, by = "Index")

# -------------------
# 12. Plot one alpha diversity panel
# -------------------
plot_alpha_index <- function(index_name) {
  
  df_sub <- alpha_long %>%
    filter(Index == index_name)
  
  p_lab <- unique(df_sub$p_label)
  
  y_max <- max(df_sub$Value, na.rm = TRUE)
  y_min <- min(df_sub$Value, na.rm = TRUE)
  y_range <- y_max - y_min
  
  if (is.na(y_range) || y_range == 0) {
    y_range <- y_max * 0.2
  }
  
  ggplot(
    df_sub,
    aes(x = Condition_label, y = Value, fill = Condition)
  ) +
    stat_summary(
      fun = mean,
      geom = "bar",
      width = 0.62,
      color = "black",
      linewidth = 0.5
    ) +
    stat_summary(
      fun.data = mean_se,
      geom = "errorbar",
      width = 0.18,
      linewidth = 0.7
    ) +
    geom_jitter(
      aes(color = Condition),
      width = 0.08,
      size = 1.6,
      alpha = 0.55,
      show.legend = FALSE
    ) +
    annotate(
      "text",
      x = 1.5,
      y = y_max + 0.15 * y_range,
      label = p_lab,
      size = 4.4,
      fontface = "bold"
    ) +
    scale_fill_manual(values = condition_fill, labels = condition_labels) +
    scale_color_manual(values = condition_outline, labels = condition_labels) +
    scale_x_discrete(drop = TRUE) +
    labs(
      title = index_name,
      x = NULL,
      y = paste0(index_name, " diversity"),
      fill = "Condition"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 17, face = "bold", hjust = 0),
      axis.title.y = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12.5, face = "bold", color = "grey25"),
      axis.text.y = element_text(size = 12, face = "bold", color = "grey25"),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    coord_cartesian(
      ylim = c(0, y_max + 0.25 * y_range),
      clip = "off"
    )
}

p_shannon  <- plot_alpha_index("Shannon")
p_simpson  <- plot_alpha_index("Simpson")
p_observed <- plot_alpha_index("Observed")
p_chao1    <- plot_alpha_index("Chao1")

p_alpha_raw <- (p_shannon | p_simpson) / (p_observed | p_chao1)
p_alpha <- add_panel_label(p_alpha_raw, "A")

print(p_alpha)

ggsave(
  filename = file.path(fig_dir, "Figure2A_alpha_diversity.png"),
  plot = p_alpha,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# ==========================================================
# Figure 2B: Bray-Curtis PCoA with average elimination rate
# ==========================================================

# -------------------
# 13. Compute average elimination rate constant
# -------------------
pfas_columns <- c(
  "k_PFOA", "k_PFPeS", "k_PFHxS", "k_PFHpS",
  "k_LPFOS", "k_PFOS_MP1", "k_PFOS_MP345", "k_PFOS_MP26"
)

available_cols <- intersect(
  pfas_columns,
  colnames(sample_data(physeq_subset))
)

if (length(available_cols) == 0) {
  stop("No PFAS elimination rate constant columns found in metadata.")
} else if (length(available_cols) < length(pfas_columns)) {
  warning(
    "Some PFAS elimination columns are missing. Using only: ",
    paste(available_cols, collapse = ", ")
  )
}

sample_data(physeq_subset)$average_k <- rowMeans(
  as.data.frame(sample_data(physeq_subset))[, available_cols, drop = FALSE],
  na.rm = TRUE
)

average_k_vec <- sample_data(physeq_subset)$average_k
average_k_vec[is.nan(average_k_vec)] <- NA
sample_data(physeq_subset)$average_k <- average_k_vec

# -------------------
# 14. Relative abundance transform for Bray-Curtis
# -------------------
physeq_rel <- transform_sample_counts(
  physeq_subset,
  function(x) x / sum(x)
)

# -------------------
# 15. Bray-Curtis distance and PCoA
# -------------------
bray_dist <- phyloseq::distance(
  physeq_rel,
  method = "bray"
)

ordination_bray <- ordinate(
  physeq_rel,
  method = "PCoA",
  distance = bray_dist
)

axis_labels <- ordination_bray$values$Relative_eig * 100

ord_df <- plot_ordination(
  physeq_rel,
  ordination_bray,
  justDF = TRUE
)

ord_df <- ord_df %>%
  filter(!is.na(Condition)) %>%
  mutate(
    Condition = factor(Condition, levels = c("PFAS_Exposed", "Reference")),
    Condition_label = recode(
      as.character(Condition),
      "PFAS_Exposed" = "PFAS-exposed",
      "Reference" = "Reference"
    ),
    Condition_label = factor(Condition_label, levels = c("PFAS-exposed", "Reference"))
  )

# -------------------
# 16. PERMANOVA
# -------------------
meta_beta <- data.frame(
  sample_data(physeq_rel),
  check.names = FALSE
)

meta_beta <- meta_beta[labels(bray_dist), , drop = FALSE]
meta_beta <- as.data.frame(meta_beta)
meta_beta$Condition <- factor(meta_beta$Condition, levels = c("PFAS_Exposed", "Reference"))

stopifnot(all(rownames(meta_beta) == labels(bray_dist)))

set.seed(123)

permanova_res <- vegan::adonis2(
  bray_dist ~ Condition,
  data = meta_beta,
  permutations = 999
)

print(permanova_res)

permanova_R2 <- permanova_res$R2[1]
permanova_p  <- permanova_res$`Pr(>F)`[1]

permanova_label <- paste0(
  "PERMANOVA: R² = ",
  signif(permanova_R2, 2),
  ", p = ",
  signif(permanova_p, 3)
)

# -------------------
# 17. Plot PCoA
# -------------------
p_pcoa_raw <- ggplot(ord_df, aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    data = ord_df %>% filter(Condition == "Reference"),
    fill = "#8A8A8A",
    color = "black",
    shape = 21,
    size = 5,
    stroke = 0.8,
    alpha = 0.9
  ) +
  geom_point(
    data = ord_df %>% filter(Condition == "PFAS_Exposed"),
    aes(fill = average_k),
    color = "#7B1FA2",
    shape = 21,
    size = 5,
    stroke = 0.9,
    alpha = 0.95
  ) +
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = permanova_label,
    hjust = 1.02,
    vjust = 1.2,
    size = 3.9,
    fontface = "bold",
    label.size = 0.25,
    fill = "white"
  ) +
  scale_fill_gradient(
    low = "#D9D9D9",
    high = "#9C27FF",
    name = "Average elimination rate",
    na.value = "#8A8A8A"
  ) +
  labs(
    x = paste0("PCoA1 (", round(axis_labels[1], 2), "%)"),
    y = paste0("PCoA2 (", round(axis_labels[2], 2), "%)")
  ) +
  guides(
    fill = guide_colorbar(
      title = "Average elimination rate",
      barwidth = 1.2,
      barheight = 5,
      title.position = "top"
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13, face = "bold", color = "grey25"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )

p_pcoa <- add_panel_label(p_pcoa_raw, "B")

print(p_pcoa)

ggsave(
  filename = file.path(fig_dir, "Figure2B_Bray_Curtis_PCoA_with_Elimination_Rate.png"),
  plot = p_pcoa,
  width = 8.5,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

# ==========================================================
# Combined Figure 2: A/B only
# ==========================================================

combined_figure2 <- p_alpha / p_pcoa +
  plot_layout(heights = c(1.05, 1.00))

print(combined_figure2)

ggsave(
  filename = file.path(fig_dir, "Figure2_Combined_AB.png"),
  plot = combined_figure2,
  width = 10,
  height = 13,
  dpi = 300,
  bg = "white"
)

cat("\nFigure 2 completed and saved to:\n", fig_dir, "\n")
cat("PERMANOVA label: ", permanova_label, "\n")



#########################################
#########################################
