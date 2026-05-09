############################################################
# 03_suppfig_s3_compound_specific_taxa_correlations.R
#
# Supplementary compound-specific PFAS/taxa correlation figure; this script also contains the paired Figure 3 workflow.
############################################################

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

