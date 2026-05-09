############################################################
# 04_suppfig_s4_species_relative_abundance_dotplot.R
#
# Supplementary Figure S4: species relative-abundance dot plot.
#
# Paths have been converted to repo-relative locations:
#   data/processed/
#   data/metadata/
#   data/model_outputs/
#   results/supplementary_figures/
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
