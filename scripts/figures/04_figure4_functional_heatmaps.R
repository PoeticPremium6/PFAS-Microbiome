# FIGURE 4 - FINAL TIGHT LAYOUT SCRIPT
#
# Fixes in this version:
# - Moves heatmap y-axis title further left so it does not overlap dendrograms.
# - Moves "Samples" axis title upward.
# - Removes the awkward "Heatmap fill" descriptor from inside the figure.
# - Makes heatmap x/y tick labels larger and bold.
# - Keeps A-E panel labels closer to panels.
# - Reduces white space between heatmaps and barplots.
# - Saves combined figure plus separate A/B and C-E figures.
############################################################
############################################################

# =========================
# 0. LOAD LIBRARIES
# =========================
library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(cowplot)
library(pheatmap)
library(ggplotify)
library(grid)

# =========================
# 1. FILE PATHS
# =========================
subsystem_path <- file.path("data", "model_outputs", "SubsystemAbundance.csv")
metadata_path  <- file.path("data", "metadata", "Sample_Metadata_New.tsv")

metabolic_data_path <- file.path("data", "model_outputs", "Metabolite_Clean_Summarized_Test.csv")
reaction_data_path  <- file.path("data", "model_outputs", "ReactionAbundance.csv")

output_dir <- file.path("results", "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# 2. USER SETTINGS
# =========================

# For main figure, 35 rows is much cleaner than 50.
# Keep 50 for a supplementary full heatmap if needed.
top_n_subsystems <- 35

alpha_raw <- 0.05
alpha_fdr <- 0.05
flux_threshold <- 1e-6

# Heatmap text/cell sizing
heatmap_font_row    <- 7.0
heatmap_font_col    <- 7.0
heatmap_font_main   <- 13
heatmap_font_legend <- 9
heatmap_cellheight  <- 8
heatmap_cellwidth   <- 8

# External heatmap axis labels and panel labels
heatmap_axis_title_size  <- 16
heatmap_panel_label_size <- 24

# Bar plot font settings
bar_axis_y_size <- 9.8
bar_axis_x_size <- 10.8
bar_title_size  <- 13.5
bar_star_size   <- 3.1

# =========================
# 3. HELPER FUNCTIONS
# =========================

read_clean_delim <- function(path, sep = "\t") {
  df <- fread(path, sep = sep, data.table = FALSE, check.names = FALSE)
  colnames(df) <- make.unique(trimws(colnames(df)))
  df
}

find_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

to_numeric_safe <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "NaN", "nd", "ND", "na", "n/a", "N/A",
             "<LOD", "<LOQ", "< LOQ", "LOD", "LOQ")] <- NA
  x <- gsub(",", ".", x, fixed = TRUE)
  x <- gsub("^<\\s*", "", x)
  suppressWarnings(as.numeric(x))
}

p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p)   ~ "",
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE       ~ ""
  )
}

wrap_axis_labels <- function(x, width = 24) {
  stringr::str_wrap(x, width = width)
}

# ------------------------------------------------------------
# Safely style pheatmap text grobs.
# This avoids: "must specify only one of 'font' and 'fontface'"
# ------------------------------------------------------------
safe_text_style <- function(grob_i,
                            fontsize = NULL,
                            fontface = "bold",
                            col = "black") {
  
  if (!is.null(grob_i$gp$font)) {
    grob_i$gp$font <- NULL
  }
  
  if (!is.null(fontsize)) {
    grob_i$gp$fontsize <- fontsize
  }
  
  if (!is.null(fontface)) {
    grob_i$gp$fontface <- fontface
  }
  
  if (!is.null(col)) {
    grob_i$gp$col <- col
  }
  
  grob_i
}

style_pheatmap_text <- function(ph,
                                row_fontsize = heatmap_font_row,
                                col_fontsize = heatmap_font_col,
                                main_fontsize = heatmap_font_main,
                                legend_fontsize = heatmap_font_legend) {
  
  g <- ph$gtable
  
  for (i in seq_along(g$grobs)) {
    
    grob_name <- g$layout$name[i]
    grob_i <- g$grobs[[i]]
    
    if (inherits(grob_i, "text")) {
      
      if (grob_name == "row_names") {
        g$grobs[[i]] <- safe_text_style(
          grob_i,
          fontsize = row_fontsize,
          fontface = "bold",
          col = "black"
        )
      } else if (grob_name == "col_names") {
        g$grobs[[i]] <- safe_text_style(
          grob_i,
          fontsize = col_fontsize,
          fontface = "bold",
          col = "black"
        )
      } else if (grob_name == "main") {
        g$grobs[[i]] <- safe_text_style(
          grob_i,
          fontsize = main_fontsize,
          fontface = "bold",
          col = "black"
        )
      } else if (grepl("annotation|legend", grob_name, ignore.case = TRUE)) {
        g$grobs[[i]] <- safe_text_style(
          grob_i,
          fontsize = legend_fontsize,
          fontface = "bold",
          col = "black"
        )
      }
    }
  }
  
  ph$gtable <- g
  ph
}

# ------------------------------------------------------------
# Wrap pheatmap with external axis labels.
# Key spacing improvements:
# - heatmap is shifted right to make room for "Subsystems"
# - "Samples" is moved upward
# - panel tag is not added here; cowplot adds A/B outside
# ------------------------------------------------------------
wrap_heatmap_plot <- function(pheatmap_obj,
                              x_label = "Samples",
                              y_label = "Subsystems") {
  
  hm_grob <- ggplotify::as.grob(pheatmap_obj)
  
  cowplot::ggdraw() +
    cowplot::draw_grob(
      hm_grob,
      x = 0.085,
      y = 0.145,
      width = 0.890,
      height = 0.815
    ) +
    cowplot::draw_label(
      x_label,
      x = 0.54,
      y = 0.075,
      fontface = "bold",
      size = heatmap_axis_title_size,
      color = "black"
    ) +
    cowplot::draw_label(
      y_label,
      x = 0.020,
      y = 0.535,
      angle = 90,
      fontface = "bold",
      size = heatmap_axis_title_size,
      color = "black"
    )
}

# ------------------------------------------------------------
# Bar plot builder
# ------------------------------------------------------------
make_bar_plot <- function(df,
                          label_col,
                          title_text,
                          fill_color,
                          xlab_text = "Mean difference\n(PFAS-exposed \u2212 Reference)") {
  
  df <- df %>%
    arrange(Difference) %>%
    mutate(Label = wrap_axis_labels(.data[[label_col]], width = 22))
  
  max_abs <- max(abs(df$Difference), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  
  ggplot(df, aes(x = reorder(Label, Difference), y = Difference)) +
    geom_col(
      fill = fill_color,
      color = "black",
      width = 0.72,
      linewidth = 0.45
    ) +
    geom_hline(
      yintercept = 0,
      color = "grey25",
      linewidth = 0.35
    ) +
    geom_text(
      aes(
        label = significance_raw,
        y = ifelse(
          Difference > 0,
          Difference + 0.045 * max_abs,
          Difference - 0.045 * max_abs
        )
      ),
      size = bar_star_size,
      fontface = "bold"
    ) +
    coord_flip(clip = "off") +
    labs(
      title = title_text,
      x = NULL,
      y = xlab_text
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.10, 0.14))) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(
        size = bar_title_size,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 5)
      ),
      axis.title.x = element_text(
        size = bar_axis_x_size,
        face = "bold",
        color = "black",
        margin = margin(t = 4)
      ),
      axis.text.x = element_text(
        size = bar_axis_x_size,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(
        size = bar_axis_y_size,
        face = "bold",
        color = "black",
        lineheight = 0.88
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 3, r = 10, b = 3, l = 3)
    )
}

# =========================
# 4. READ METADATA
# =========================

condition_data <- read_clean_delim(metadata_path, sep = "\t")

sample_col <- find_first_existing_col(
  condition_data,
  c("Sample_Name", "SampleID", "Sample_ID", "Sample")
)

condition_col <- find_first_existing_col(
  condition_data,
  c("Condition", "Group")
)

if (is.null(sample_col) || is.null(condition_col)) {
  stop("Could not find sample or condition column in metadata.")
}

condition_data[[sample_col]] <- trimws(as.character(condition_data[[sample_col]]))
condition_data[[condition_col]] <- trimws(as.character(condition_data[[condition_col]]))

k_cols <- c(
  "k_PFOA", "k_PFPeS", "k_PFHxS", "k_PFHpS",
  "k_LPFOS", "k_PFOS_MP1", "k_PFOS_MP345", "k_PFOS_MP26"
)

available_k_cols <- intersect(k_cols, colnames(condition_data))

if (length(available_k_cols) > 0) {
  condition_data$Average_k <- rowMeans(
    as.data.frame(lapply(condition_data[, available_k_cols, drop = FALSE], to_numeric_safe)),
    na.rm = TRUE
  )
} else {
  condition_data$Average_k <- NA_real_
}

condition_data <- condition_data %>%
  mutate(
    Condition = .data[[condition_col]],
    Sample_Name = .data[[sample_col]]
  )

# =========================
# 5. READ SUBSYSTEM DATA
# =========================

subsystem_data <- fread(subsystem_path, data.table = FALSE, check.names = FALSE)
colnames(subsystem_data) <- make.unique(trimws(colnames(subsystem_data)))

colnames(subsystem_data)[1] <- "Subsystem"

numeric_cols <- setdiff(colnames(subsystem_data), "Subsystem")

subsystem_data[numeric_cols] <- lapply(
  subsystem_data[numeric_cols],
  to_numeric_safe
)

subsystem_data$total_abundance <- rowSums(
  subsystem_data[, numeric_cols, drop = FALSE],
  na.rm = TRUE
)

top_subsystems <- subsystem_data %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = top_n_subsystems) %>%
  select(-total_abundance)

# =========================
# 6. BUILD HEATMAP MATRICES
# =========================

pfas_samples <- condition_data %>%
  filter(Condition == "PFAS_Exposed") %>%
  pull(Sample_Name)

ref_samples <- condition_data %>%
  filter(Condition == "Reference") %>%
  pull(Sample_Name)

pfas_samples <- intersect(pfas_samples, colnames(top_subsystems))
ref_samples  <- intersect(ref_samples, colnames(top_subsystems))

if (length(pfas_samples) == 0) stop("No PFAS-exposed samples matched subsystem table.")
if (length(ref_samples) == 0)  stop("No reference samples matched subsystem table.")

pfas_data <- top_subsystems[, c("Subsystem", pfas_samples), drop = FALSE]
ref_data  <- top_subsystems[, c("Subsystem", ref_samples), drop = FALSE]

pfas_matrix <- as.matrix(pfas_data[, -1, drop = FALSE])
ref_matrix  <- as.matrix(ref_data[, -1, drop = FALSE])

rownames(pfas_matrix) <- pfas_data$Subsystem
rownames(ref_matrix)  <- ref_data$Subsystem

pfas_matrix[is.na(pfas_matrix)] <- 0
ref_matrix[is.na(ref_matrix)]   <- 0

combined_matrix <- cbind(pfas_matrix, ref_matrix)

scaled_combined <- t(scale(t(combined_matrix)))
scaled_combined[is.na(scaled_combined)] <- 0
scaled_combined[is.nan(scaled_combined)] <- 0
scaled_combined[is.infinite(scaled_combined)] <- 0

scaled_pfas <- scaled_combined[, colnames(pfas_matrix), drop = FALSE]
scaled_ref  <- scaled_combined[, colnames(ref_matrix), drop = FALSE]

# =========================
# 7. HEATMAP ANNOTATIONS
# =========================

pfas_ann <- condition_data %>%
  filter(Sample_Name %in% colnames(scaled_pfas)) %>%
  select(Sample_Name, Condition, Average_k)

pfas_ann <- pfas_ann[match(colnames(scaled_pfas), pfas_ann$Sample_Name), ]
rownames(pfas_ann) <- pfas_ann$Sample_Name
pfas_ann$Sample_Name <- NULL

ref_ann <- condition_data %>%
  filter(Sample_Name %in% colnames(scaled_ref)) %>%
  select(Sample_Name, Condition)

ref_ann <- ref_ann[match(colnames(scaled_ref), ref_ann$Sample_Name), ]
rownames(ref_ann) <- ref_ann$Sample_Name
ref_ann$Sample_Name <- NULL

ann_colors <- list(
  Condition = c("PFAS_Exposed" = "#7B2CBF", "Reference" = "#2D5BCE")
)

heat_colors <- colorRampPalette(
  c("#313695", "#74add1", "#ffffbf", "#fdae61", "#a50026")
)(100)

heat_breaks <- seq(-3, 3, length.out = 101)

# =========================
# 8. BUILD HEATMAPS
# =========================

pheatmap_pfas <- pheatmap(
  scaled_pfas,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = pfas_ann,
  annotation_colors = ann_colors,
  color = heat_colors,
  breaks = heat_breaks,
  main = "PFAS-exposed",
  fontsize_row = heatmap_font_row,
  fontsize_col = heatmap_font_col,
  fontsize = heatmap_font_legend,
  angle_col = 90,
  cellheight = heatmap_cellheight,
  cellwidth = heatmap_cellwidth,
  border_color = "grey65",
  legend = TRUE,
  legend_breaks = c(-3, -1.5, 0, 1.5, 3),
  legend_labels = c("-3", "-1.5", "0", "1.5", "3"),
  silent = TRUE
)

pheatmap_pfas <- style_pheatmap_text(
  pheatmap_pfas,
  row_fontsize = heatmap_font_row,
  col_fontsize = heatmap_font_col,
  main_fontsize = heatmap_font_main,
  legend_fontsize = heatmap_font_legend
)

pheatmap_ref <- pheatmap(
  scaled_ref,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = ref_ann,
  annotation_colors = ann_colors,
  color = heat_colors,
  breaks = heat_breaks,
  main = "Reference",
  fontsize_row = heatmap_font_row,
  fontsize_col = heatmap_font_col,
  fontsize = heatmap_font_legend,
  angle_col = 90,
  cellheight = heatmap_cellheight,
  cellwidth = heatmap_cellwidth,
  border_color = "grey65",
  legend = TRUE,
  legend_breaks = c(-3, -1.5, 0, 1.5, 3),
  legend_labels = c("-3", "-1.5", "0", "1.5", "3"),
  silent = TRUE
)

pheatmap_ref <- style_pheatmap_text(
  pheatmap_ref,
  row_fontsize = heatmap_font_row,
  col_fontsize = heatmap_font_col,
  main_fontsize = heatmap_font_main,
  legend_fontsize = heatmap_font_legend
)

heatmap_A <- wrap_heatmap_plot(
  pheatmap_pfas,
  x_label = "Samples",
  y_label = "Subsystems"
)

heatmap_B <- wrap_heatmap_plot(
  pheatmap_ref,
  x_label = "Samples",
  y_label = "Subsystems"
)

heatmap_row <- cowplot::plot_grid(
  heatmap_A,
  heatmap_B,
  ncol = 2,
  rel_widths = c(1.33, 1.07),
  labels = c("A", "B"),
  label_size = heatmap_panel_label_size,
  label_fontface = "bold",
  label_x = c(0.006, 0.006),
  label_y = c(0.995, 0.995),
  hjust = 0,
  vjust = 1,
  align = "h",
  axis = "tb",
  greedy = TRUE
)

# =========================
# 9. METABOLITE ANALYSIS
# =========================

metabolic_data <- read.csv(
  metabolic_data_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(metabolic_data) <- make.unique(trimws(colnames(metabolic_data)))

metabolic_long <- metabolic_data %>%
  pivot_longer(
    -Metabolite,
    names_to = "Sample_Name",
    values_to = "Value"
  ) %>%
  mutate(Value = to_numeric_safe(Value))

biomass_values <- metabolic_long %>%
  filter(Metabolite == "microbeBiomass[fe]") %>%
  select(Sample_Name, Biomass = Value)

normalized_metabolites <- metabolic_long %>%
  inner_join(biomass_values, by = "Sample_Name") %>%
  mutate(Normalized_Value = Value / Biomass) %>%
  filter(Metabolite != "microbeBiomass[fe]") %>%
  inner_join(
    condition_data %>% select(Sample_Name, Condition),
    by = "Sample_Name"
  )

metabolite_tests <- normalized_metabolites %>%
  group_by(Metabolite) %>%
  summarise(
    p_value = tryCatch(
      t.test(Normalized_Value ~ Condition)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(p_value, method = "BH"),
    significance_raw = p_to_stars(p_value)
  )

metabolite_summary <- normalized_metabolites %>%
  group_by(Metabolite, Condition) %>%
  summarise(
    Mean_Value = mean(Normalized_Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Condition, values_from = Mean_Value) %>%
  mutate(Difference = PFAS_Exposed - Reference) %>%
  inner_join(metabolite_tests, by = "Metabolite") %>%
  filter(
    !is.na(Difference),
    abs(Difference) >= flux_threshold,
    p_value < alpha_raw
  )

top_pos_met <- metabolite_summary %>%
  filter(Difference > 0) %>%
  arrange(desc(Difference)) %>%
  slice_head(n = 10)

top_pos_met$Metabolite <- recode(
  top_pos_met$Metabolite,
  "glc_D[fe]" = "D-Glucose (Extracellular)",
  "no[fe]" = "Nitric oxide (Extracellular)",
  "ichor[u]" = "Isochorismate (Cytosol)",
  "4abz[fe]" = "4-Aminobenzoate (Extracellular)",
  "cmp[u]" = "Cytidine monophosphate",
  "ichor[fe]" = "Isochorismate (Extracellular)",
  "glygn5[fe]" = "Glycogen Structure 5 (Glycogenin)",
  "srtn[fe]" = "Serotonin (Extracellular)",
  "5htrp[fe]" = "5-hydroxy-L-tryptophan (Extracellular)",
  "5htrp[u]" = "5-hydroxy-L-tryptophan (Cytosol)"
)

# =========================
# 10. REACTION ANALYSIS
# =========================

reaction_data <- read.csv(
  reaction_data_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(reaction_data) <- make.unique(trimws(colnames(reaction_data)))

reaction_long <- reaction_data %>%
  pivot_longer(
    -Reaction,
    names_to = "Sample_Name",
    values_to = "Abundance"
  ) %>%
  mutate(Abundance = to_numeric_safe(Abundance)) %>%
  inner_join(
    condition_data %>% select(Sample_Name, Condition),
    by = "Sample_Name"
  )

reaction_tests <- reaction_long %>%
  group_by(Reaction) %>%
  summarise(
    p_value = tryCatch(
      t.test(Abundance ~ Condition)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(p_value, method = "BH"),
    significance_raw = p_to_stars(p_value)
  )

reaction_summary <- reaction_long %>%
  group_by(Reaction, Condition) %>%
  summarise(
    Mean_Abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Condition, values_from = Mean_Abundance) %>%
  mutate(Difference = PFAS_Exposed - Reference) %>%
  inner_join(reaction_tests, by = "Reaction") %>%
  filter(
    !is.na(Difference),
    abs(Difference) >= flux_threshold,
    p_value < alpha_raw
  )

top_neg_rxn <- reaction_summary %>%
  filter(Difference < 0) %>%
  arrange(Difference) %>%
  slice_head(n = 10)

top_pos_rxn <- reaction_summary %>%
  filter(Difference > 0) %>%
  arrange(desc(Difference)) %>%
  slice_head(n = 10)

pretty_rxn_names <- c(
  "ALATA_L" = "L-Alanine Transaminase",
  "SPMS" = "Spermidine Synthase",
  "LDH_L" = "L-Lactate Dehydrogenase",
  "EX_4abz(e)" = "4-Aminobenzoate Exchange",
  "L_LACD" = "L-Lactate Dehydrogenase",
  "L_LACD2" = "L-Lactate Dehydrogenase (ubiquinone)",
  "L_LACD3" = "L-Lactate Dehydrogenase (menaquinone)",
  "ADMDC" = "Adenosylmethionine Decarboxylase",
  "EX_ala_L(e)" = "Exchange of L-Alanine",
  "EAR11M12" = "nad[c]P+ oxidoreductase (A-specific)",
  "RMK2" = "Rhamnulokinase 2",
  "RMK" = "Rhamnulokinase",
  "EX_ppi(e)" = "Exchange of PPI",
  "EX_r406(e)" = "Unknown Exchange Reaction (r406)",
  "EX_r788(e)" = "Unknown Exchange Reaction (r788)",
  "CDPGHL" = "CDPglucose 4,6-Hydro-Lyase",
  "G3PD8i" = "Glycerol 3 phosphate dehydrogenase",
  "BUTKr" = "Butyrate Kinase",
  "COBCOMOX" = "Coenzyme B-M-Phenazine Oxidoreductase",
  "PBUTT" = "Phosphate Butyryltransferase"
)

top_neg_rxn$Reaction <- dplyr::recode(top_neg_rxn$Reaction, !!!pretty_rxn_names)
top_pos_rxn$Reaction <- dplyr::recode(top_pos_rxn$Reaction, !!!pretty_rxn_names)

# =========================
# 11. BUILD BAR PLOTS
# =========================

plot_C <- make_bar_plot(
  top_pos_met,
  label_col = "Metabolite",
  title_text = "Top Positive Metabolites",
  fill_color = "#6A5ACD"
)

plot_D <- make_bar_plot(
  top_neg_rxn,
  label_col = "Reaction",
  title_text = "Top Negative Reactions",
  fill_color = "#FA8072"
)

plot_E <- make_bar_plot(
  top_pos_rxn,
  label_col = "Reaction",
  title_text = "Top Positive Reactions",
  fill_color = "#6A5ACD"
)

bar_row <- cowplot::plot_grid(
  plot_C,
  plot_D,
  plot_E,
  ncol = 3,
  rel_widths = c(1, 1, 1),
  labels = c("C", "D", "E"),
  label_size = 22,
  label_fontface = "bold",
  label_x = c(0.006, 0.006, 0.006),
  label_y = c(0.995, 0.995, 0.995),
  hjust = 0,
  vjust = 1,
  align = "h",
  axis = "tb",
  greedy = TRUE
)

# =========================
# 12. FINAL FIGURE COMBINATION
# =========================

final_figure4 <- cowplot::plot_grid(
  heatmap_row,
  bar_row,
  ncol = 1,
  rel_heights = c(1.02, 0.90),
  align = "v",
  axis = "lr",
  greedy = TRUE
)

# =========================
# 13. SAVE OUTPUTS
# =========================

ggsave(
  file.path(output_dir, "Figure4_Combined_Final_TightSpacing.png"),
  plot = final_figure4,
  width = 18,
  height = 12.8,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(output_dir, "Figure4AB_SubsystemHeatmaps_TightSpacing.png"),
  plot = heatmap_row,
  width = 18,
  height = 7.2,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  file.path(output_dir, "Figure4CDE_Barplots_TightSpacing.png"),
  plot = bar_row,
  width = 18,
  height = 5.9,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)

# =========================
# 14. OPTIONAL EXPORT TABLES
# =========================

write.csv(
  top_subsystems,
  file.path(output_dir, "Figure4_TopSubsystems_Input.csv"),
  row.names = FALSE
)

write.csv(
  top_pos_met,
  file.path(output_dir, "Figure4_TopPositiveMetabolites.csv"),
  row.names = FALSE
)

write.csv(
  top_neg_rxn,
  file.path(output_dir, "Figure4_TopNegativeReactions.csv"),
  row.names = FALSE
)

write.csv(
  top_pos_rxn,
  file.path(output_dir, "Figure4_TopPositiveReactions.csv"),
  row.names = FALSE
)

# =========================
# 15. PRINT AND SUMMARY
# =========================

print(final_figure4)

cat("\nSaved:\n",
    file.path(output_dir, "Figure4_Combined_Final_TightSpacing.png"), "\n",
    file.path(output_dir, "Figure4AB_SubsystemHeatmaps_TightSpacing.png"), "\n",
    file.path(output_dir, "Figure4CDE_Barplots_TightSpacing.png"), "\n")

cat("\nHeatmap fill represents row z-score abundance.\n")
cat("Bar plot stars indicate raw p-values; BH-FDR values are available in exported summaries if needed.\n")


#####################
#####################
