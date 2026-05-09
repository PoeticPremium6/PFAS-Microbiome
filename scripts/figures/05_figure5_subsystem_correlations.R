# =====================================================================
# Main text figure-generation code for PFAS microbiome manuscript
#
# This public version uses repository-relative paths only.
# Expected input files should be placed under data/processed,
# data/metadata, and data/model_outputs as described in docs/figure_workflow.md.
# Outputs are written to results/figures.
# =====================================================================

# Figure 5 — final revised script with larger fonts/dots
#
# Subplot A:
#   subsystem abundance ~ log10 blood:feces PFAS ratio
#   log10((blood + pseudocount) / (feces + pseudocount))
#
# Subplot B:
#   subsystem abundance ~ PFAS elimination rate constants
#
# Visual strategy:
# - Lollipop plots instead of bars.
# - X-axis shows Spearman rho.
# - Point shape shows raw p-value category.
# - Unused p-value categories are dropped from legends.
# - Panel B color indicates PFAS type.
# - Larger font, larger dots, thicker lollipop stems.
# - Reduced white space in combined figure.
#
# Statistical note:
# - Displayed features are selected using raw p <= 0.05.
# - FDR-adjusted p-values are exported but not used for figure inclusion.
#####################
#####################

# ================================================================
# 0. LIBRARIES
# ================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(data.table)
library(stringr)
library(cowplot)
library(readr)
library(tibble)
library(forcats)
library(grid)

# ================================================================
# 1. OUTPUT DIRECTORY
# ================================================================

out_dir <- file.path("results", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ================================================================
# 2. INPUT FILES
# ================================================================

subsystem_file <- file.path("data", "model_outputs", "SubsystemAbundance.csv")
metadata_file  <- file.path("data", "metadata", "Sample_Metadata_New.tsv")

# ================================================================
# 3. SETTINGS
# ================================================================

alpha_raw <- 0.05
top_n_B_per_k <- 5

pfas_palette <- c(
  "k_PFOA"       = "#000000",
  "k_PFPeS"      = "#E69F00",
  "k_PFHxS"      = "#56B4E9",
  "k_PFHpS"      = "#009E73",
  "k_LPFOS"      = "#F0E442",
  "k_PFOS_MP1"   = "#0072B2",
  "k_PFOS_MP345" = "#D55E00",
  "k_PFOS_MP26"  = "#CC79A7"
)

pretty_k_labels <- c(
  "k_PFOA"       = "k_PFOA",
  "k_PFPeS"      = "k_PFPeS",
  "k_PFHxS"      = "k_PFHxS",
  "k_PFHpS"      = "k_PFHpS",
  "k_LPFOS"      = "k_LPFOS",
  "k_PFOS_MP1"   = "k_PFOS_MP1",
  "k_PFOS_MP345" = "k_PFOS_MP345",
  "k_PFOS_MP26"  = "k_PFOS_MP26"
)

direction_palette <- c(
  "Positive" = "#4D4D4D",
  "Negative" = "#B3B3B3"
)

sig_shape_values <- c(
  "p ≤ 0.05"  = 16,  # circle
  "p ≤ 0.01"  = 17,  # triangle
  "p ≤ 0.001" = 18   # diamond
)

# Manual vertical dodging for Panel B.
# Increase slightly if PFAS-specific points are still too close.
b_dodge_step <- 0.165

# Larger dots and thicker lines
dot_size_A <- 5.2
dot_size_B <- 4.6

line_width_A <- 1.65
line_width_B <- 1.45

base_font_size <- 15

# ================================================================
# 4. HELPER FUNCTIONS
# ================================================================

repair_names <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x)] <- "X"
  make.unique(x, sep = "_dup")
}

read_tsv_repaired <- function(path) {
  df <- data.table::fread(
    path,
    sep = "\t",
    check.names = FALSE,
    data.table = FALSE
  )
  names(df) <- repair_names(names(df))
  tibble::as_tibble(df, .name_repair = "minimal")
}

read_csv_repaired <- function(path) {
  df <- data.table::fread(
    path,
    check.names = FALSE,
    data.table = FALSE
  )
  names(df) <- repair_names(names(df))
  tibble::as_tibble(df, .name_repair = "minimal")
}

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

sig_label <- function(p) {
  dplyr::case_when(
    is.na(p)    ~ "",
    p <= 0.001  ~ "***",
    p <= 0.01   ~ "**",
    p <= 0.05   ~ "*",
    TRUE        ~ ""
  )
}

sig_shape_label <- function(p) {
  dplyr::case_when(
    is.na(p)    ~ NA_character_,
    p <= 0.001  ~ "p ≤ 0.001",
    p <= 0.01   ~ "p ≤ 0.01",
    p <= 0.05   ~ "p ≤ 0.05",
    TRUE        ~ NA_character_
  )
}

safe_spearman <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 4) {
    return(tibble(Correlation = NA_real_, P_Value = NA_real_, N = length(x)))
  }
  
  if (sd(x, na.rm = TRUE) == 0 || sd(y, na.rm = TRUE) == 0) {
    return(tibble(Correlation = NA_real_, P_Value = NA_real_, N = length(x)))
  }
  
  test <- suppressWarnings(
    tryCatch(
      cor.test(x, y, method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
  )
  
  if (is.null(test)) {
    return(tibble(Correlation = NA_real_, P_Value = NA_real_, N = length(x)))
  }
  
  tibble(
    Correlation = unname(test$estimate),
    P_Value = test$p.value,
    N = length(x)
  )
}

theme_pub <- function(base_size = base_font_size) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      axis.title = element_text(face = "bold", color = "black", size = base_size + 2),
      axis.text.x = element_text(face = "bold", color = "black", size = base_size),
      axis.text.y = element_text(face = "bold", color = "black", size = base_size),
      strip.text = element_text(face = "bold", size = base_size + 2, color = "black"),
      strip.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.40),
      panel.grid.major.x = element_line(color = "grey88", linewidth = 0.40),
      legend.title = element_text(face = "bold", color = "black", size = base_size + 1),
      legend.text = element_text(face = "bold", color = "black", size = base_size),
      plot.margin = margin(4, 6, 6, 4)
    )
}

# ================================================================
# 5. LOAD DATA
# ================================================================

subsystem_data <- read_csv_repaired(subsystem_file)
names(subsystem_data)[1] <- "Subsystems"

pfas_metadata <- read_tsv_repaired(metadata_file)

cat("\nMetadata column names repaired with make.unique().\n")
cat("Example columns containing BAN_feces:\n")
print(grep("BAN_feces", names(pfas_metadata), value = TRUE))

if ("Sample_Name" %in% names(pfas_metadata)) {
  pfas_metadata <- pfas_metadata %>%
    mutate(Sample = trimws(as.character(Sample_Name)))
} else if ("Sample" %in% names(pfas_metadata)) {
  pfas_metadata <- pfas_metadata %>%
    mutate(Sample = trimws(as.character(Sample)))
} else if ("SampleID" %in% names(pfas_metadata)) {
  pfas_metadata <- pfas_metadata %>%
    mutate(Sample = trimws(as.character(SampleID)))
} else {
  pfas_metadata <- pfas_metadata %>%
    rename(Sample = 1) %>%
    mutate(Sample = trimws(as.character(Sample)))
}

sample_names <- intersect(names(subsystem_data)[-1], pfas_metadata$Sample)

if (length(sample_names) == 0) {
  stop("No matching sample names between subsystem data columns and metadata Sample column.")
}

subsystem_data <- subsystem_data %>%
  dplyr::select(Subsystems, all_of(sample_names)) %>%
  mutate(across(-Subsystems, to_numeric_safe))

pfas_metadata <- pfas_metadata %>%
  dplyr::filter(Sample %in% sample_names)

# ================================================================
# 6. LONG DATA
# ================================================================

data_long <- subsystem_data %>%
  pivot_longer(
    -Subsystems,
    names_to = "Sample",
    values_to = "Abundance"
  ) %>%
  left_join(pfas_metadata, by = "Sample") %>%
  mutate(
    Subsystems = as.character(Subsystems),
    Abundance = to_numeric_safe(Abundance)
  )

# ================================================================
# 7. SUBPLOT A — LOG10 BLOOD:FECES PFAS RATIO
# ================================================================

pfas_pairs <- list(
  PFHxS = c("PFHxS1", "fecal_PFHxS"),
  PFOA  = c("PFOA1",  "fecal_PFOA"),
  PFNA  = c("PFNA1",  "fecal_PFNA"),
  PFOS  = c("PFOS1",  "fecal_L_PFOS")
)

missing_pair_cols <- unique(unlist(pfas_pairs))[!unique(unlist(pfas_pairs)) %in% names(data_long)]

if (length(missing_pair_cols) > 0) {
  warning(
    "These blood/feces PFAS columns are missing and will be skipped: ",
    paste(missing_pair_cols, collapse = ", ")
  )
}

plotA_all <- purrr::map_dfr(names(pfas_pairs), function(pfas) {
  
  blood <- pfas_pairs[[pfas]][1]
  feces <- pfas_pairs[[pfas]][2]
  
  if (!all(c(blood, feces) %in% names(data_long))) {
    return(tibble())
  }
  
  dat <- data_long %>%
    mutate(
      Blood = to_numeric_safe(.data[[blood]]),
      Feces = to_numeric_safe(.data[[feces]])
    ) %>%
    filter(!is.na(Blood), !is.na(Feces), !is.na(Abundance))
  
  nonzero_values <- c(dat$Blood[dat$Blood > 0], dat$Feces[dat$Feces > 0])
  
  if (length(nonzero_values) == 0) {
    warning("No positive non-zero blood/feces values found for ", pfas, ". Skipping.")
    return(tibble())
  }
  
  pseudocount <- min(nonzero_values, na.rm = TRUE) / 2
  
  dat <- dat %>%
    mutate(
      Log_Blood_Feces_Ratio = log10((Blood + pseudocount) / (Feces + pseudocount))
    ) %>%
    filter(is.finite(Log_Blood_Feces_Ratio))
  
  dat %>%
    group_by(Subsystems) %>%
    group_modify(~ safe_spearman(.x$Abundance, .x$Log_Blood_Feces_Ratio)) %>%
    ungroup() %>%
    mutate(
      PFAS = pfas,
      Ratio_Metric = "log10((blood + pseudocount) / (feces + pseudocount))",
      Pseudocount = pseudocount
    )
}) %>%
  filter(!is.na(P_Value), !is.na(Correlation)) %>%
  group_by(PFAS) %>%
  mutate(
    P_adj_BH = p.adjust(P_Value, method = "BH"),
    Sig_raw = sig_label(P_Value),
    Sig_fdr = sig_label(P_adj_BH),
    Sig_shape = sig_shape_label(P_Value)
  ) %>%
  ungroup()

plotA_data <- plotA_all %>%
  filter(P_Value <= alpha_raw) %>%
  group_by(PFAS) %>%
  arrange(desc(abs(Correlation)), .by_group = TRUE) %>%
  ungroup()

write_csv(
  plotA_all %>%
    transmute(
      PFAS,
      Subsystem = Subsystems,
      Correlation,
      P_Value,
      P_adj_BH,
      Sig_raw,
      Sig_fdr,
      N,
      Ratio_Metric,
      Pseudocount
    ),
  file.path(out_dir, "Figure5A_ALL_LogBloodFecesRatio_Subsystems_withFDR.csv")
)

write_csv(
  plotA_data %>%
    transmute(
      PFAS,
      Subsystem = Subsystems,
      Correlation,
      P_Value,
      P_adj_BH,
      Sig_raw,
      Sig_fdr,
      N,
      Ratio_Metric,
      Pseudocount
    ),
  file.path(out_dir, "Figure5A_PLOTTED_LogBloodFecesRatio_Subsystems_RAWp_le0.05_withFDR.csv")
)

if (nrow(plotA_data) == 0) {
  
  plot_A <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No subsystem correlations at raw p ≤ 0.05",
      fontface = "bold"
    ) +
    theme_void()
  
} else {
  
  plotA_plot <- plotA_data %>%
    mutate(
      Subsystems_wrapped = str_wrap(Subsystems, 34),
      Direction = ifelse(Correlation >= 0, "Positive", "Negative"),
      Sig_shape = factor(Sig_shape)
    ) %>%
    droplevels() %>%
    group_by(PFAS) %>%
    mutate(
      Subsystems_wrapped = fct_reorder(Subsystems_wrapped, Correlation)
    ) %>%
    ungroup()
  
  shape_values_A <- sig_shape_values[
    names(sig_shape_values) %in% levels(plotA_plot$Sig_shape)
  ]
  
  A_lim <- max(abs(plotA_plot$Correlation), na.rm = TRUE)
  A_lim <- max(A_lim, 0.05)
  
  plot_A <- ggplot(
    plotA_plot,
    aes(y = Subsystems_wrapped)
  ) +
    geom_segment(
      aes(
        x = 0,
        xend = Correlation,
        yend = Subsystems_wrapped,
        color = Direction
      ),
      linewidth = line_width_A,
      alpha = 0.85
    ) +
    geom_point(
      aes(
        x = Correlation,
        color = Direction,
        shape = Sig_shape
      ),
      size = dot_size_A,
      alpha = 0.95
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.55
    ) +
    scale_color_manual(
      values = direction_palette,
      name = "Direction"
    ) +
    scale_shape_manual(
      values = shape_values_A,
      name = "Raw p-value"
    ) +
    scale_x_continuous(
      limits = c(-A_lim * 1.16, A_lim * 1.16),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    facet_wrap(~ PFAS, ncol = 2, scales = "free_y") +
    labs(
      x = "Spearman correlation with log10 blood:feces PFAS ratio",
      y = "Metabolic subsystem"
    ) +
    guides(
      color = "none",
      shape = guide_legend(
        order = 1,
        override.aes = list(color = "grey35", size = 5.2)
      )
    ) +
    theme_pub(base_size = base_font_size) +
    theme(
      legend.position = "right",
      axis.title.x = element_text(margin = margin(t = 6)),
      axis.title.y = element_text(margin = margin(r = 6)),
      panel.spacing = unit(0.45, "lines"),
      plot.margin = margin(2, 4, 3, 2)
    )
}

# ================================================================
# 8. SUBPLOT B — ELIMINATION RATE CONSTANTS
# ================================================================

elimination_cols <- grep("^k_", names(pfas_metadata), value = TRUE)
elimination_cols <- elimination_cols[elimination_cols %in% names(pfas_palette)]

if (length(elimination_cols) == 0) {
  stop("No k_* elimination-rate columns found in metadata.")
}

plotB_all <- purrr::map_dfr(elimination_cols, function(k) {
  
  dat <- data_long %>%
    mutate(K_value = to_numeric_safe(.data[[k]])) %>%
    filter(!is.na(K_value), !is.na(Abundance))
  
  dat %>%
    group_by(Subsystems) %>%
    group_modify(~ safe_spearman(.x$Abundance, .x$K_value)) %>%
    ungroup() %>%
    mutate(PFAS = k)
}) %>%
  filter(!is.na(P_Value), !is.na(Correlation)) %>%
  group_by(PFAS) %>%
  mutate(
    P_adj_BH = p.adjust(P_Value, method = "BH"),
    Sig_raw = sig_label(P_Value),
    Sig_fdr = sig_label(P_adj_BH),
    Sig_shape = sig_shape_label(P_Value)
  ) %>%
  ungroup()

plotB_data <- plotB_all %>%
  filter(P_Value <= alpha_raw) %>%
  group_by(PFAS) %>%
  slice_max(
    order_by = abs(Correlation),
    n = top_n_B_per_k,
    with_ties = FALSE
  ) %>%
  ungroup()

write_csv(
  plotB_all %>%
    transmute(
      PFAS,
      Subsystem = Subsystems,
      Correlation,
      P_Value,
      P_adj_BH,
      Sig_raw,
      Sig_fdr,
      N
    ),
  file.path(out_dir, "Figure5B_ALL_EliminationRate_Subsystems_withFDR.csv")
)

write_csv(
  plotB_data %>%
    transmute(
      PFAS,
      Subsystem = Subsystems,
      Correlation,
      P_Value,
      P_adj_BH,
      Sig_raw,
      Sig_fdr,
      N
    ),
  file.path(out_dir, "Figure5B_PLOTTED_EliminationRate_Subsystems_Top5_RAWp_le0.05_withFDR.csv")
)

if (nrow(plotB_data) == 0) {
  
  plot_B <- ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0,
      label = "No subsystem correlations at raw p ≤ 0.05",
      fontface = "bold"
    ) +
    theme_void()
  
} else {
  
  subsystem_order <- plotB_data %>%
    group_by(Subsystems) %>%
    summarise(
      max_abs = max(abs(Correlation), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(max_abs)) %>%
    pull(Subsystems)
  
  plotB_plot <- plotB_data %>%
    mutate(
      PFAS = factor(PFAS, levels = elimination_cols),
      PFAS_index = as.numeric(PFAS),
      PFAS_offset = (PFAS_index - mean(seq_along(elimination_cols))) * b_dodge_step,
      Subsystems_wrapped = str_wrap(Subsystems, 30),
      Subsystems_wrapped = factor(
        Subsystems_wrapped,
        levels = str_wrap(rev(subsystem_order), 30)
      ),
      y_base = as.numeric(Subsystems_wrapped),
      y_pos = y_base + PFAS_offset,
      Sig_shape = factor(Sig_shape)
    ) %>%
    droplevels()
  
  shape_values_B <- sig_shape_values[
    names(sig_shape_values) %in% levels(plotB_plot$Sig_shape)
  ]
  
  B_lim <- max(abs(plotB_plot$Correlation), na.rm = TRUE)
  B_lim <- max(B_lim, 0.10)
  
  y_breaks <- seq_along(levels(plotB_plot$Subsystems_wrapped))
  y_labels <- levels(plotB_plot$Subsystems_wrapped)
  
  plot_B <- ggplot(plotB_plot) +
    geom_segment(
      aes(
        x = 0,
        xend = Correlation,
        y = y_pos,
        yend = y_pos,
        color = PFAS
      ),
      linewidth = line_width_B,
      alpha = 0.72
    ) +
    geom_point(
      aes(
        x = Correlation,
        y = y_pos,
        color = PFAS,
        shape = Sig_shape
      ),
      size = dot_size_B,
      alpha = 0.95
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.55
    ) +
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      expand = expansion(mult = c(0.020, 0.030))
    ) +
    scale_color_manual(
      values = pfas_palette,
      labels = pretty_k_labels[names(pretty_k_labels) %in% elimination_cols],
      name = "PFAS type"
    ) +
    scale_shape_manual(
      values = shape_values_B,
      name = "Raw p-value"
    ) +
    scale_x_continuous(
      limits = c(-B_lim * 1.16, B_lim * 1.16),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      x = "Spearman correlation with elimination rate",
      y = "Metabolic subsystem"
    ) +
    guides(
      color = guide_legend(
        order = 1,
        override.aes = list(size = 4.8, linewidth = 1.4, shape = 16)
      ),
      shape = guide_legend(
        order = 2,
        override.aes = list(color = "grey35", size = 5.2)
      )
    ) +
    theme_pub(base_size = base_font_size) +
    theme(
      legend.position = "right",
      axis.title.x = element_text(margin = margin(t = 7)),
      axis.title.y = element_text(margin = margin(r = 6)),
      panel.spacing = unit(0.4, "lines"),
      plot.margin = margin(2, 4, 4, 2),
      legend.box.spacing = unit(0.15, "cm"),
      legend.spacing.y = unit(0.08, "cm")
    )
}

# ================================================================
# 9. PANEL LABELS AND FINAL FIGURE
# ================================================================

plot_A_labeled <- ggdraw(plot_A) +
  draw_label(
    "A",
    x = 0.005,
    y = 0.995,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 30
  )

plot_B_labeled <- ggdraw(plot_B) +
  draw_label(
    "B",
    x = 0.005,
    y = 0.995,
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 30
  )

final_figure <- plot_grid(
  plot_A_labeled,
  plot_B_labeled,
  ncol = 1,
  rel_heights = c(0.88, 1.45),
  align = "v",
  axis = "lr"
)

print(final_figure)

# ================================================================
# 10. SAVE OUTPUTS
# ================================================================

ggsave(
  filename = file.path(out_dir, "Figure5_Final_LargerTextDots.png"),
  plot = final_figure,
  width = 10.5,
  height = 12.4,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(out_dir, "Figure5_Final_LargerTextDots.pdf"),
  plot = final_figure,
  width = 10.5,
  height = 12.4,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(out_dir, "Figure5A_Final_LargerTextDots.png"),
  plot = plot_A_labeled,
  width = 10.5,
  height = 4.8,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(out_dir, "Figure5B_Final_LargerTextDots.png"),
  plot = plot_B_labeled,
  width = 10.5,
  height = 7.6,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

# ================================================================
# 11. SUMMARY
# ================================================================

message("Wrote outputs to: ", out_dir)
message(" - Figure5_Final_LargerTextDots.png / Figure5_Final_LargerTextDots.pdf")
message(" - Figure5A_Final_LargerTextDots.png")
message(" - Figure5B_Final_LargerTextDots.png")
message(" - Figure5A_ALL_LogBloodFecesRatio_Subsystems_withFDR.csv")
message(" - Figure5A_PLOTTED_LogBloodFecesRatio_Subsystems_RAWp_le0.05_withFDR.csv")
message(" - Figure5B_ALL_EliminationRate_Subsystems_withFDR.csv")
message(" - Figure5B_PLOTTED_EliminationRate_Subsystems_Top5_RAWp_le0.05_withFDR.csv")

##########################
##########################
##########################
############################################################
############################################################
