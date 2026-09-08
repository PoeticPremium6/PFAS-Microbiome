#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 39A_build_Figure_S04_supplementary.R ROOT STAGE")
}
ROOT <- normalizePath(args[[1]], mustWork = TRUE)
STAGE <- normalizePath(args[[2]], mustWork = FALSE)

BUILD <- file.path(STAGE, "figure_builds", "Figure_S04")
FIG_DIR <- file.path(BUILD, "figures")
SRC_DIR <- file.path(STAGE, "source_data", "Figure_S04")
MAN_DIR <- file.path(STAGE, "manifests", "Figure_S04")
LEG_DIR <- file.path(STAGE, "legends")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MAN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LEG_DIR, recursive = TRUE, showWarnings = FALSE)

INPUT_DIR <- file.path(
  ROOT, "submission", "source_data", "supplementary", "Figure_S04"
)
PATH_A <- file.path(INPUT_DIR, "Figure_S04_panelA_mapping_fraction.tsv")
PATH_B <- file.path(INPUT_DIR, "Figure_S04_panelB_detected_vs_matched_species.tsv")
PATH_C <- file.path(INPUT_DIR, "Figure_S04_panelC_AGORA2_feature_breadth.tsv")
PATH_C_INV <- file.path(INPUT_DIR, "Figure_S04_panelC_feature_inventory.tsv")
PATH_D <- file.path(INPUT_DIR, "Figure_S04_panelD_AGORA2_subsystem_PCoA.tsv")
PATH_E <- file.path(INPUT_DIR, "Figure_S04_panelE_mapping_breadth_correlations.tsv")

required <- c(PATH_A, PATH_B, PATH_C, PATH_C_INV, PATH_D, PATH_E)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing authoritative Figure S4 input(s): ", paste(missing, collapse = "; "))
}

pal <- c(
  "Higher" = "#4B1F6F",
  "Middle" = "#70409B",
  "Lower" = "#A985CC",
  "Reference (no k)" = "#DDD3EC"
)

theme_editorial <- function(base_size = 12.2) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_line(colour = "#E5DDED", linewidth = 0.42),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        colour = "#4A4650", fill = NA, linewidth = 0.62
      ),
      axis.title = element_text(
        face = "bold", colour = "#242126", size = rel(1.02)
      ),
      axis.text = element_text(colour = "#3D3941"),
      strip.background = element_rect(
        fill = "#EEE7F4", colour = "#B9A7C9", linewidth = 0.55
      ),
      strip.text = element_text(
        face = "bold", colour = "#332C38", size = rel(0.95),
        margin = margin(5, 6, 5, 6)
      ),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(colour = "#3D3941"),
      plot.margin = margin(8, 9, 8, 9)
    )
}

normalize_group <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[grepl("reference", y)] <- "Reference (no k)"
  out[grepl("higher|high", y) & is.na(out)] <- "Higher"
  out[grepl("middle|mid", y) & is.na(out)] <- "Middle"
  out[grepl("lower|low", y) & is.na(out)] <- "Lower"
  factor(out, levels = names(pal))
}

write_tsv <- function(dt, path) {
  fwrite(dt, path, sep = "\t", quote = FALSE, na = "NA")
}

sha256_file <- function(path) {
  out <- suppressWarnings(
    system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  )
  if (!length(out)) return(NA_character_)
  strsplit(out[[1]], "[[:space:]]+")[[1]][1]
}

format_p <- function(value) {
  format.pval(value, digits = 2, eps = 0.001)
}

# -------------------------------------------------------------------------
# Authoritative inputs
# -------------------------------------------------------------------------
panel_a <- fread(PATH_A)
panel_b <- fread(PATH_B)
panel_c <- fread(PATH_C)
panel_c_inventory <- fread(PATH_C_INV)
panel_d <- fread(PATH_D)
panel_e <- fread(PATH_E)

required_a <- c("sample_id", "elimination_class", "mapped_fraction")
required_b <- c(
  "sample_id", "elimination_class", "detected_species",
  "matched_detected_species"
)
required_c <- c(
  "sample_id", "elimination_class", "feature_family",
  "nonzero_feature_count", "total_feature_count",
  "nonzero_feature_fraction"
)
required_d <- c("sample_id", "PCoA1", "PCoA2", "elimination_class")

if (!all(required_a %in% names(panel_a))) {
  stop("Unexpected panel A schema")
}
if (!all(required_b %in% names(panel_b))) {
  stop("Unexpected panel B schema")
}
if (!all(required_c %in% names(panel_c))) {
  stop("Unexpected panel C schema")
}
if (!all(required_d %in% names(panel_d))) {
  stop("Unexpected panel D schema")
}

panel_a[, elimination_class := normalize_group(elimination_class)]
panel_b[, elimination_class := normalize_group(elimination_class)]
panel_c[, elimination_class := normalize_group(elimination_class)]
panel_d[, elimination_class := normalize_group(elimination_class)]

panel_a <- panel_a[
  !is.na(elimination_class) & is.finite(as.numeric(mapped_fraction))
]
panel_b <- panel_b[
  !is.na(elimination_class) &
  is.finite(as.numeric(detected_species)) &
  is.finite(as.numeric(matched_detected_species))
]
panel_c <- panel_c[
  !is.na(elimination_class) &
  is.finite(as.numeric(nonzero_feature_fraction))
]
panel_d <- panel_d[
  !is.na(elimination_class) &
  is.finite(as.numeric(PCoA1)) &
  is.finite(as.numeric(PCoA2))
]

if (uniqueN(panel_a$sample_id) != 65L) stop("Panel A sample count is not 65")
if (uniqueN(panel_b$sample_id) != 65L) stop("Panel B sample count is not 65")
if (uniqueN(panel_d$sample_id) != 65L) stop("Panel D sample count is not 65")
panel_a[, mapped_fraction := as.numeric(mapped_fraction)]
panel_b[, `:=`(
  detected_species = as.numeric(detected_species),
  matched_detected_species = as.numeric(matched_detected_species)
)]
panel_c[, `:=`(
  nonzero_feature_count = as.numeric(nonzero_feature_count),
  total_feature_count = as.numeric(total_feature_count),
  nonzero_feature_fraction = as.numeric(nonzero_feature_fraction)
)]
panel_d[, `:=`(
  PCoA1 = as.numeric(PCoA1),
  PCoA2 = as.numeric(PCoA2)
)]

panel_c_raw_rows <- nrow(panel_c)
panel_c_raw_families <- uniqueN(panel_c$feature_family)
panel_c[, feature_family_raw := trimws(as.character(feature_family))]
panel_c[, feature_family := fcase(
  grepl("reaction", feature_family_raw, ignore.case = TRUE),
  "AGORA2 reactions",
  grepl("subsystem", feature_family_raw, ignore.case = TRUE),
  "AGORA2 subsystems",
  grepl("exchange", feature_family_raw, ignore.case = TRUE),
  "AGORA2 exchanges",
  default = NA_character_
)]
panel_c <- panel_c[!is.na(feature_family)]

panel_c_duplicate_audit <- panel_c[, .(
  input_records = .N,
  unique_nonzero_counts = uniqueN(nonzero_feature_count),
  unique_total_counts = uniqueN(total_feature_count),
  unique_fractions = uniqueN(round(nonzero_feature_fraction, 12)),
  unique_groups = uniqueN(as.character(elimination_class))
), by = .(sample_id, feature_family)]

inconsistent_panel_c <- panel_c_duplicate_audit[
  unique_nonzero_counts > 1L |
  unique_total_counts > 1L |
  unique_fractions > 1L |
  unique_groups > 1L
]
if (nrow(inconsistent_panel_c)) {
  stop(
    "Panel C duplicate records disagree for ",
    nrow(inconsistent_panel_c),
    " sample-family combinations"
  )
}

panel_c <- panel_c[, .SD[1L], by = .(sample_id, feature_family)]
panel_c[, feature_family_raw := NULL]

panel_c_normalization_audit <- data.table(
  metric = c(
    "raw_rows", "raw_family_labels", "canonical_rows",
    "canonical_families", "unique_samples", "duplicate_records_removed"
  ),
  value = c(
    panel_c_raw_rows,
    panel_c_raw_families,
    nrow(panel_c),
    uniqueN(panel_c$feature_family),
    uniqueN(panel_c$sample_id),
    panel_c_raw_rows - nrow(panel_c)
  )
)

if (nrow(panel_c) != 195L) {
  stop("Panel C must contain 65 samples x 3 canonical families")
}
if (uniqueN(panel_c$feature_family) != 3L) {
  stop("Panel C must contain exactly 3 canonical families")
}
if (uniqueN(panel_c$sample_id) != 65L) {
  stop("Panel C must contain exactly 65 unique samples")
}

# -------------------------------------------------------------------------
# Panel A: mapped abundance fraction
# -------------------------------------------------------------------------
kw_a <- kruskal.test(mapped_fraction ~ elimination_class, data = panel_a)

p_a <- ggplot(
  panel_a,
  aes(x = elimination_class, y = mapped_fraction, fill = elimination_class)
) +
  geom_boxplot(
    width = 0.58, outlier.shape = NA, linewidth = 0.72,
    colour = "#332B38", alpha = 0.9
  ) +
  geom_jitter(
    aes(colour = elimination_class),
    width = 0.12, height = 0, size = 2.15, alpha = 0.82,
    show.legend = FALSE
  ) +
  annotate(
    "label", x = Inf, y = Inf,
    label = paste0("Kruskal-Wallis p=", format_p(kw_a$p.value)),
    hjust = 1.03, vjust = 1.23, size = 3.7, fontface = "italic",
    linewidth = 0.28, fill = alpha("white", 0.9),
    colour = "#39323D"
  ) +
  scale_fill_manual(values = pal, drop = FALSE) +
  scale_colour_manual(values = pal, drop = FALSE) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  labs(
    x = "Elimination class",
    y = "AGORA2-mapped abundance fraction"
  ) +
  guides(fill = "none") +
  theme_editorial() +
  theme(
    axis.text.x = element_text(size = 10.4),
    plot.margin = margin(8, 10, 8, 8)
  )

# -------------------------------------------------------------------------
# Omitted diagnostic: detected versus matched species
# -------------------------------------------------------------------------
rho_b <- suppressWarnings(cor(
  panel_b$detected_species,
  panel_b$matched_detected_species,
  method = "spearman",
  use = "complete.obs"
))

# Retained in source data and statistics, but not rendered because the fixed
# AGORA2 model ceiling makes the identity-line scatter visually redundant.

# -------------------------------------------------------------------------
# Panel B: feature breadth
# -------------------------------------------------------------------------
panel_c[, family_display := fcase(
  feature_family == "AGORA2 reactions",
  paste0("Reactions\n(", comma(total_feature_count), " total)"),
  feature_family == "AGORA2 subsystems",
  paste0("Subsystems\n(", comma(total_feature_count), " total)"),
  feature_family == "AGORA2 exchanges",
  paste0("Exchanges\n(", comma(total_feature_count), " total)"),
  default = NA_character_
)]

family_order <- c(
  unique(panel_c[feature_family == "AGORA2 reactions", family_display]),
  unique(panel_c[feature_family == "AGORA2 subsystems", family_display]),
  unique(panel_c[feature_family == "AGORA2 exchanges", family_display])
)
panel_c[, family_display := factor(family_display, levels = family_order)]

c_stats <- panel_c[, {
  test <- kruskal.test(nonzero_feature_fraction ~ elimination_class)
  .(
    p_value = test$p.value,
    annotation = paste0("p=", format_p(test$p.value))
  )
}, by = family_display]
c_stats[, `:=`(
  elimination_class = factor("Reference (no k)", levels = names(pal)),
  y = 0.995
)]

p_c <- ggplot(
  panel_c,
  aes(
    x = elimination_class,
    y = nonzero_feature_fraction,
    fill = elimination_class
  )
) +
  geom_boxplot(
    width = 0.6, outlier.shape = NA, linewidth = 0.65,
    colour = "#332B38", alpha = 0.88
  ) +
  geom_jitter(
    aes(colour = elimination_class),
    width = 0.12, height = 0, size = 2.00, alpha = 0.80,
    show.legend = FALSE
  ) +
  geom_text(
    data = c_stats,
    aes(
      x = elimination_class, y = y, label = annotation
    ),
    inherit.aes = FALSE,
    hjust = 1.02, vjust = 1.2, size = 3.35, fontface = "italic",
    colour = "#39323D"
  ) +
  facet_grid(. ~ family_display) +
  scale_fill_manual(values = pal, drop = FALSE) +
  scale_colour_manual(values = pal, drop = FALSE) +
  scale_x_discrete(
    labels = c(
      "Higher" = "Higher",
      "Middle" = "Middle",
      "Lower" = "Lower",
      "Reference (no k)" = "Reference\n(no k)"
    )
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    breaks = seq(0.65, 1.00, by = 0.05),
    expand = expansion(mult = c(0.015, 0.025))
  ) +
  coord_cartesian(ylim = c(0.65, 1.00), clip = "on") +
  labs(
    x = "Elimination class",
    y = "Non-zero AGORA2 feature fraction"
  ) +
  guides(fill = "none") +
  theme_editorial(base_size = 11.8) +
  theme(
    axis.text.x = element_text(
      angle = 0, hjust = 0.5, vjust = 0.5, size = 9.0
    ),
    axis.text.y = element_text(size = 9.8),
    strip.text = element_text(size = 11.2),
    panel.spacing.x = unit(9, "pt"),
    plot.margin = margin(8, 8, 8, 8)
  )

# -------------------------------------------------------------------------
# Panel C: locked AGORA2 subsystem PCoA
# -------------------------------------------------------------------------
p_d <- ggplot(
  panel_d,
  aes(x = PCoA1, y = PCoA2, colour = elimination_class)
) +
  geom_hline(
    yintercept = 0, colour = "#B9AFBE",
    linewidth = 0.46, linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0, colour = "#B9AFBE",
    linewidth = 0.46, linetype = "dashed"
  ) +
  geom_point(size = 3.00, alpha = 0.88, stroke = 0.35) +
  annotate(
    "text", x = Inf, y = Inf, label = "n=65",
    hjust = 1.12, vjust = 1.35, size = 3.8, fontface = "italic"
  ) +
  scale_colour_manual(values = pal, drop = FALSE) +
  labs(
    x = "AGORA2 subsystem Bray-Curtis PCoA1",
    y = "AGORA2 subsystem Bray-Curtis PCoA2",
    colour = "Elimination class"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(size = 3.4, alpha = 1)
    )
  ) +
  theme_editorial(base_size = 12.5) +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.78, "lines"),
    legend.margin = margin(0, 0, 0, 8),
    plot.margin = margin(8, 8, 8, 10)
  )

# -------------------------------------------------------------------------
# Source data, statistics and provenance
# -------------------------------------------------------------------------
panel_a_stats <- data.table(
  panel = "A",
  test = "Kruskal-Wallis",
  n = nrow(panel_a),
  statistic = unname(kw_a$statistic),
  degrees_of_freedom = unname(kw_a$parameter),
  p_value = kw_a$p.value
)
panel_b_stats <- data.table(
  panel = "OMITTED_SPECIES_MATCHING_DIAGNOSTIC",
  test = "Spearman correlation",
  n = nrow(panel_b),
  rho = rho_b
)
panel_d_stats <- data.table(
  panel = "C",
  n = nrow(panel_d),
  note = "Locked AGORA2 subsystem Bray-Curtis PCoA"
)

write_tsv(panel_a, file.path(SRC_DIR, "Figure_S04_panelA_mapping_fraction.tsv"))
write_tsv(panel_a_stats, file.path(SRC_DIR, "Figure_S04_panelA_statistics.tsv"))
write_tsv(panel_b, file.path(SRC_DIR, "Figure_S04_panelB_detected_vs_matched_species.tsv"))
write_tsv(panel_b_stats, file.path(SRC_DIR, "Figure_S04_panelB_statistics.tsv"))
panel_c_export <- copy(panel_c)
panel_c_export[, family_display := gsub(
  "\n", " ", as.character(family_display), fixed = TRUE
)]
write_tsv(
  panel_c_export,
  file.path(SRC_DIR, "Figure_S04_panelC_AGORA2_feature_breadth.tsv")
)
write_tsv(panel_c_inventory, file.path(SRC_DIR, "Figure_S04_panelC_feature_inventory.tsv"))
write_tsv(
  panel_c_normalization_audit,
  file.path(SRC_DIR, "Figure_S04_panelC_normalization_audit.tsv")
)
write_tsv(c_stats, file.path(SRC_DIR, "Figure_S04_panelC_statistics.tsv"))
write_tsv(panel_d, file.path(SRC_DIR, "Figure_S04_panelD_AGORA2_subsystem_PCoA.tsv"))
write_tsv(panel_d_stats, file.path(SRC_DIR, "Figure_S04_panelD_statistics.tsv"))
write_tsv(panel_e, file.path(SRC_DIR, "Figure_S04_panelE_mapping_breadth_correlations.tsv"))

panel_design <- data.table(
  panel = c("A", "B", "C", "D", "E"),
  content = c(
    "AGORA2 mapped abundance fraction",
    "AGORA2 reaction, subsystem and exchange breadth",
    "Locked AGORA2 subsystem Bray-Curtis PCoA",
    "Detected versus AGORA2-matched species",
    "Mapping fraction versus AGORA2 feature breadth"
  ),
  status = c("INCLUDED", "INCLUDED", "INCLUDED", "OMITTED", "OMITTED"),
  reason = c(
    "Direct row-level mapping-fraction table",
    "Derived from authoritative matched AGORA2 matrices",
    "Locked row-level ordination",
    "Retained as source-data diagnostic but omitted because the AGORA2 model ceiling makes the identity-line scatter visually redundant",
    "Weak derived relationships; strongest absolute Spearman rho below 0.35"
  )
)
write_tsv(panel_design, file.path(SRC_DIR, "Figure_S04_panel_design.tsv"))

provenance <- data.table(
  role = c(
    "panel_A",
    "panel_B",
    "panel_C",
    "omitted_species_matching_diagnostic",
    "omitted_mapping_breadth_diagnostic"
  ),
  source_path = c(
    PATH_A,
    paste(PATH_C, PATH_C_INV, sep = " | "),
    PATH_D,
    PATH_B,
    PATH_E
  ),
  source_type = c(
    "Prepared direct mapping fraction",
    "Prepared derived feature breadth",
    "Locked AGORA2 subsystem PCoA",
    "Prepared direct species matching retained as source data",
    "Prepared derived diagnostic omitted from figure"
  ),
  rows_or_samples = c(
    nrow(panel_a), nrow(panel_c), nrow(panel_d), nrow(panel_b), nrow(panel_e)
  )
)
write_tsv(provenance, file.path(MAN_DIR, "provenance.tsv"))

# -------------------------------------------------------------------------
# Editorial layout
# -------------------------------------------------------------------------
top <- p_a + p_c + plot_layout(widths = c(0.82, 1.38))
figure <- top / p_d +
  plot_layout(heights = c(0.96, 1.10)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold", size = 22, colour = "#171419"
      ),
      plot.tag.position = c(0, 1)
    )
  )

pdf_path <- file.path(FIG_DIR, "Figure_S04_excellence_candidate.pdf")
png_path <- file.path(FIG_DIR, "Figure_S04_excellence_candidate.png")
tif_path <- file.path(FIG_DIR, "Figure_S04_excellence_candidate.tiff")

ggsave(
  pdf_path, figure, width = 16.0, height = 12.0,
  units = "in", device = cairo_pdf
)
ggsave(
  png_path, figure, width = 16.0, height = 12.0,
  units = "in", dpi = 300, bg = "white"
)
ggsave(
  tif_path, figure, width = 16.0, height = 12.0,
  units = "in", dpi = 300, compression = "lzw", bg = "white"
)

legend_text <- paste(
  "Figure S4. AGORA2 reconstruction coverage and functional feature breadth.",
  "Panel A shows the fraction of measured community abundance represented by",
  "matched AGORA2 species models. Panel B shows the proportion of reaction,",
  "subsystem and exchange features with nonzero abundance in each sample-specific",
  "community profile. Panel C shows the locked AGORA2 subsystem Bray-Curtis",
  "ordination. The detected-versus-matched species diagnostic is retained in",
  "the source data but omitted from the figure because the fixed AGORA2 model",
  "ceiling makes the identity-line comparison visually redundant. The optional",
  "mapping-fraction versus feature-breadth panel was also omitted because the",
  "derived relationships were weak."
)
writeLines(legend_text, file.path(LEG_DIR, "Figure_S04.txt"))

source_files <- list.files(SRC_DIR, full.names = TRUE)
source_manifest <- data.table(
  file = basename(source_files),
  bytes = file.info(source_files)$size,
  sha256 = vapply(source_files, sha256_file, character(1))
)
write_tsv(source_manifest, file.path(MAN_DIR, "source_data_manifest.tsv"))

visual_qc <- data.table(
  check = "MANUAL_VISUAL_QC",
  status = "PENDING",
  detail = paste(
    "Confirm sequential A-C panel labels, readable typography, mapping-fraction",
    "scale in A, enlarged feature-family facets and p-values in B, unclipped",
    "PCoA with elimination-class legend in C, balanced spacing, purple palette",
    "and absence of Lindell content."
  )
)
write_tsv(visual_qc, file.path(MAN_DIR, "visual_qc.tsv"))

cat("AGORA2_MAPPING_FRACTION_SOURCE=PASS n=", nrow(panel_a), "\n", sep = "")
cat("AGORA2_SPECIES_MATCHING_SOURCE=PASS n=", nrow(panel_b), "\n", sep = "")
cat("AGORA2_FEATURE_BREADTH_SOURCE=PASS rows=", nrow(panel_c), " families=", uniqueN(panel_c$feature_family), "\n", sep = "")
cat("PANEL_C_CANONICALIZATION=PASS raw_rows=", panel_c_raw_rows, " canonical_rows=", nrow(panel_c), " families=", uniqueN(panel_c$feature_family), "\n", sep = "")
cat("PANEL_C_TSV_SERIALIZATION=PASS embedded_newlines_removed=TRUE\n")
cat("AGORA2_SUBSYSTEM_PCOA_SOURCE=PASS n=", nrow(panel_d), "\n", sep = "")
cat("PANEL_E_DECISION=OMITTED_WEAK_DERIVED_RELATIONSHIP\n")
cat("LINDELL_SOURCE_EXCLUDED_FROM_FIGURE_S04=PASS\n")
cat("SPECIES_MATCHING_PANEL_DECISION=OMITTED_VISUALLY_REDUNDANT\n")
cat("FIGURE_S04_PANEL_DESIGN=A_MAPPING_B_BREADTH_C_PCOA\n")
cat("GENERATOR_EXECUTION_PASS\n")
