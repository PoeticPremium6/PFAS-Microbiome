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
  stop("Usage: 36A_build_Figure_S03_supplementary.R ROOT STAGE")
}
ROOT <- normalizePath(args[[1]], mustWork = TRUE)
STAGE <- normalizePath(args[[2]], mustWork = FALSE)

BUILD <- file.path(STAGE, "figure_builds", "Figure_S03")
FIG_DIR <- file.path(BUILD, "figures")
SRC_DIR <- file.path(STAGE, "source_data", "Figure_S03")
MAN_DIR <- file.path(STAGE, "manifests", "Figure_S03")
LEG_DIR <- file.path(STAGE, "legends")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MAN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LEG_DIR, recursive = TRUE, showWarnings = FALSE)

PATH_MATRIX <- file.path(
  ROOT,
  "08_microbiome_statistics/00_inputs/corrected_feature_matrices/",
  "humann_pathways_sample_matched.tsv"
)
PATH_PCOA <- file.path(
  ROOT,
  "11_figures_LOCKED/figure_set_locked_20260714_090127/source_data/",
  "Figure_3C_HUMAnN_scaled_PCoA.tsv"
)
PATH_HUMANN_DA <- file.path(
  ROOT,
  "20_item2_DA_diversity_audit/elimination_corrected_rerun_20260730_203216/",
  "results/task_00908_CLR_LM_humann_pathways_k_LPFOS_tertile/result.tsv"
)
PATH_AGORA_DA <- file.path(
  ROOT,
  "20_item2_DA_diversity_audit/elimination_corrected_rerun_20260730_203216/",
  "results/task_00620_CLR_LM_agora2_subsystem_k_LPFOS_tertile/result.tsv"
)

required <- c(PATH_MATRIX, PATH_PCOA, PATH_HUMANN_DA, PATH_AGORA_DA)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop("Missing authoritative Figure S3 input(s): ", paste(missing, collapse = "; "))
}

pal <- c(
  "Higher" = "#4B1F6F",
  "Middle" = "#70409B",
  "Lower" = "#A985CC",
  "Reference (no k)" = "#DDD3EC"
)
support_fill <- c(
  "FDR q<0.05" = "#4B1F6F",
  "Nominal p<0.05" = "#7650A0",
  "No nominal support" = "white"
)
support_shape <- c(
  "FDR q<0.05" = 21,
  "Nominal p<0.05" = 22,
  "No nominal support" = 23
)

theme_editorial <- function(base_size = 10.5) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_line(colour = "#E5DDED", linewidth = 0.38),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "#4A4650", fill = NA, linewidth = 0.55),
      axis.title = element_text(face = "bold", colour = "#242126", size = rel(1.05)),
      axis.text = element_text(colour = "#3D3941"),
      strip.background = element_rect(fill = "#EEE7F4", colour = "#B9A7C9", linewidth = 0.5),
      strip.text = element_text(face = "bold", colour = "#332C38", margin = margin(4, 5, 4, 5)),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(colour = "#3D3941"),
      plot.margin = margin(6, 7, 6, 7)
    )
}

normalize_group <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[grepl("reference", y)] <- "Reference (no k)"
  out[grepl("higher|high", y) & is.na(out)] <- "Higher"
  out[grepl("middle|mid", y) & is.na(out)] <- "Middle"
  out[grepl("lower|low", y) & is.na(out)] <- "Lower"
  out
}

safe_min <- function(x, default = NA_real_) {
  x <- x[is.finite(x)]
  if (length(x)) min(x) else default
}

sha256_file <- function(path) {
  out <- suppressWarnings(system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE))
  if (!length(out)) return(NA_character_)
  strsplit(out[[1]], "[[:space:]]+")[[1]][1]
}

write_tsv <- function(dt, path) {
  fwrite(dt, path, sep = "\t", quote = FALSE, na = "NA")
}

# -------------------------------------------------------------------------
# Panel A: functional Shannon diversity derived from HUMAnN pathway matrix
# -------------------------------------------------------------------------
pcoa <- fread(PATH_PCOA)
required_pcoa <- c("sample_id", "PCoA1", "PCoA2", "elimination_class")
if (!all(required_pcoa %in% names(pcoa))) {
  stop("Unexpected locked HUMAnN PCoA schema: ", paste(names(pcoa), collapse = ","))
}
pcoa[, elimination_class := normalize_group(elimination_class)]
pcoa <- pcoa[!is.na(elimination_class)]
pcoa[, elimination_class := factor(
  elimination_class,
  levels = c("Higher", "Middle", "Lower", "Reference (no k)")
)]
if (uniqueN(pcoa$sample_id) < 60L) {
  stop("Locked HUMAnN PCoA cohort unexpectedly small")
}

mat <- fread(PATH_MATRIX, check.names = FALSE)
feature_col <- names(mat)[1]
sample_cols <- intersect(as.character(pcoa$sample_id), names(mat))
if (length(sample_cols) < 60L) {
  stop("Insufficient matched HUMAnN matrix samples: ", length(sample_cols))
}

features <- as.character(mat[[feature_col]])
keep_feature <- !grepl(
  "UNMAPPED|UNINTEGRATED|UNCLASSIFIED",
  features,
  ignore.case = TRUE
)
mat_use <- mat[keep_feature]
x <- as.matrix(mat_use[, ..sample_cols])
suppressWarnings(storage.mode(x) <- "numeric")
x[!is.finite(x) | x < 0] <- 0

sample_totals <- colSums(x)
observed <- colSums(x > 0)
shannon <- vapply(seq_along(sample_cols), function(i) {
  values <- x[, i]
  total <- sum(values)
  if (!is.finite(total) || total <= 0) return(NA_real_)
  proportions <- values[values > 0] / total
  -sum(proportions * log(proportions))
}, numeric(1))

alpha <- data.table(
  sample_id = sample_cols,
  shannon_diversity = shannon,
  observed_pathway_richness = observed,
  pathway_abundance_sum = sample_totals
)
alpha <- merge(
  alpha,
  pcoa[, .(sample_id, elimination_class)],
  by = "sample_id",
  all.x = TRUE
)
alpha <- alpha[is.finite(shannon_diversity) & !is.na(elimination_class)]
setcolorder(alpha, c(
  "sample_id", "elimination_class", "shannon_diversity",
  "observed_pathway_richness", "pathway_abundance_sum"
))
if (nrow(alpha) < 60L) stop("Functional alpha derivation retained too few samples")

kw <- kruskal.test(shannon_diversity ~ elimination_class, data = alpha)
alpha_stats <- data.table(
  test = "Kruskal-Wallis",
  metric = "HUMAnN pathway Shannon diversity",
  n = nrow(alpha),
  groups = uniqueN(alpha$elimination_class),
  statistic = unname(kw$statistic),
  degrees_of_freedom = unname(kw$parameter),
  p_value = kw$p.value,
  source_basis = "Derived from matched HUMAnN pathway abundance matrix"
)

p_a <- ggplot(alpha, aes(x = elimination_class, y = shannon_diversity, fill = elimination_class)) +
  geom_boxplot(
    width = 0.58, outlier.shape = NA, linewidth = 0.65,
    colour = "#332B38", alpha = 0.88
  ) +
  geom_jitter(
    aes(colour = elimination_class),
    width = 0.12, height = 0, size = 1.9, alpha = 0.82,
    show.legend = FALSE
  ) +
  annotate(
    "label", x = Inf, y = Inf,
    label = sprintf("Kruskal-Wallis p=%s", format.pval(kw$p.value, digits = 2, eps = 0.001)),
    hjust = 1.03, vjust = 1.25, size = 3.15, fontface = "italic",
    label.size = 0.25, fill = alpha("white", 0.88), colour = "#39323D"
  ) +
  scale_fill_manual(values = pal, drop = FALSE) +
  scale_colour_manual(values = pal, drop = FALSE) +
  labs(
    x = "Elimination class",
    y = "HUMAnN pathway Shannon diversity"
  ) +
  guides(fill = "none") +
  theme_editorial() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# -------------------------------------------------------------------------
# Panel B: locked HUMAnN PCoA
# -------------------------------------------------------------------------
p_b <- ggplot(pcoa, aes(x = PCoA1, y = PCoA2, colour = elimination_class)) +
  geom_hline(yintercept = 0, colour = "#B9AFBE", linewidth = 0.42, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "#B9AFBE", linewidth = 0.42, linetype = "dashed") +
  geom_point(size = 2.45, alpha = 0.88, stroke = 0.35) +
  annotate(
    "text", x = Inf, y = Inf, label = paste0("n=", nrow(pcoa)),
    hjust = 1.12, vjust = 1.35, size = 3.35, fontface = "italic"
  ) +
  scale_colour_manual(values = pal, drop = FALSE) +
  labs(
    x = "HUMAnN Bray-Curtis PCoA1",
    y = "HUMAnN Bray-Curtis PCoA2",
    colour = "Elimination class"
  ) +
  theme_editorial() +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.55, "lines"),
    legend.margin = margin(0, 0, 0, 3)
  )

# -------------------------------------------------------------------------
# Functional association preparation
# -------------------------------------------------------------------------
humanize_pathway <- function(x) {
  y <- as.character(x)
  y <- sub("\\|.*$", "", y)
  y <- sub("^[A-Za-z0-9-]+:[[:space:]]*", "", y)
  y <- gsub("_", " ", y, fixed = TRUE)
  y <- gsub("&alpha;", "alpha", y, fixed = TRUE)
  y <- gsub("&beta;", "beta", y, fixed = TRUE)
  y <- gsub("&gamma;", "gamma", y, fixed = TRUE)
  y <- gsub("&delta;", "delta", y, fixed = TRUE)
  y <- gsub("[[:space:]]+", " ", y)
  y <- trimws(y)
  y <- sub("^superpathway of ", "", y, ignore.case = TRUE)
  y <- gsub("\\s*\\([^)]*\\)$", "", y)
  y <- gsub("[[:space:]]+", " ", y)
  trimws(y)
}

humanize_subsystem <- function(x) {
  y <- gsub("_", " ", as.character(x), fixed = TRUE)
  y <- gsub("[[:space:]]+", " ", y)
  trimws(y)
}

contrast_label <- function(term) {
  z <- tolower(as.character(term))
  out <- rep(NA_character_, length(z))
  out[grepl("low", z)] <- "Low T1 vs High T3"
  out[grepl("mid|middle", z)] <- "Middle T2 vs High T3"
  out
}

prepare_da <- function(path, family_pattern, label_fun, top_n, panel_code) {
  dt <- fread(path)
  required_cols <- c("feature", "predictor", "term", "estimate", "std_error", "p_value", "q_value")
  if (!all(required_cols %in% names(dt))) {
    stop("Unexpected ", panel_code, " association schema: ", paste(names(dt), collapse = ","))
  }

  dt <- dt[tolower(predictor) == "k_lpfos_tertile"]
  if ("feature_family" %in% names(dt)) {
    dt <- dt[grepl(family_pattern, feature_family, ignore.case = TRUE)]
  }
  dt[, contrast := contrast_label(term)]
  dt <- dt[!is.na(contrast)]
  dt <- dt[
    is.finite(as.numeric(estimate)) &
    is.finite(as.numeric(std_error)) &
    is.finite(as.numeric(p_value))
  ]
  if (!nrow(dt)) stop("No usable ", panel_code, " target rows")

  dt[, `:=`(
    estimate = as.numeric(estimate),
    std_error = as.numeric(std_error),
    p_value = as.numeric(p_value),
    q_value = as.numeric(q_value)
  )]
  dt[, feature_label := label_fun(feature)]
  dt <- dt[nzchar(feature_label)]
  if (panel_code == "HUMAnN") {
    dt <- dt[!grepl("\\|", feature)]
    dt <- dt[!grepl("UNMAPPED|UNINTEGRATED|UNCLASSIFIED", feature, ignore.case = TRUE)]
  }

  dt[, support := fifelse(
    is.finite(q_value) & q_value < 0.05, "FDR q<0.05",
    fifelse(p_value < 0.05, "Nominal p<0.05", "No nominal support")
  )]
  dt[, support_rank := fifelse(
    support == "FDR q<0.05", 1L,
    fifelse(support == "Nominal p<0.05", 2L, 3L)
  )]
  dt[, `:=`(
    ci_low = estimate - 1.96 * std_error,
    ci_high = estimate + 1.96 * std_error
  )]

  dt[, abs_estimate := abs(estimate)]
  setorder(dt, feature, contrast, support_rank, p_value, -abs_estimate)
  dt <- dt[, .SD[1], by = .(feature, contrast)]

  ranking <- dt[, .(
    support_rank = min(support_rank),
    min_q = safe_min(q_value, 1),
    min_p = safe_min(p_value, 1),
    max_abs_effect = max(abs(estimate), na.rm = TRUE)
  ), by = .(feature, feature_label)]
  setorder(ranking, support_rank, min_q, min_p, -max_abs_effect)
  ranking <- head(ranking, top_n)

  display <- dt[feature %in% ranking$feature]
  display <- merge(
    display,
    ranking[, .(feature, display_order = .I)],
    by = "feature",
    all.x = TRUE
  )
  label_map <- unique(display[, .(feature, feature_label, display_order)])
  setorder(label_map, display_order)
  label_map[, feature_label_unique := make.unique(feature_label)]
  display <- merge(
    display,
    label_map[, .(feature, feature_label_unique)],
    by = "feature",
    all.x = TRUE
  )
  display[, feature_label_unique := factor(
    feature_label_unique,
    levels = rev(label_map$feature_label_unique)
  )]
  display[, contrast := factor(
    contrast,
    levels = c("Low T1 vs High T3", "Middle T2 vs High T3")
  )]
  display[, support := factor(
    support,
    levels = c("FDR q<0.05", "Nominal p<0.05", "No nominal support")
  )]

  list(full = dt, display = display, ranking = ranking)
}

humann_da <- prepare_da(
  PATH_HUMANN_DA, "humann.*pathway|humann_pathway",
  humanize_pathway, 12L, "HUMAnN"
)
agora_da <- prepare_da(
  PATH_AGORA_DA, "agora2.*subsystem|agora2_subsystem",
  humanize_subsystem, 10L, "AGORA2"
)

support_levels_to_show <- {
  ordered <- c("FDR q<0.05", "Nominal p<0.05", "No nominal support")
  present <- unique(c(as.character(humann_da$display$support), as.character(agora_da$display$support)))
  keep <- ordered[ordered %in% present]
  if (!length(keep)) "No nominal support" else keep
}
SUPPORT_LEGEND_DYNAMIC=TRUE

forest_plot <- function(display, y_title, show_legend = TRUE, support_levels = support_levels_to_show) {
  max_ci <- max(abs(c(display$ci_low, display$ci_high)), na.rm = TRUE)
  max_ci <- max(5, ceiling(max_ci / 5) * 5)
  display <- copy(display)
  display[, support := factor(as.character(support), levels = support_levels)]
  shape_values <- support_shape[support_levels]
  fill_values <- support_fill[support_levels]
  ggplot(display, aes(y = feature_label_unique, x = estimate)) +
    geom_vline(xintercept = 0, colour = "#4B4650", linewidth = 0.62) +
    geom_segment(
      aes(x = ci_low, xend = ci_high, yend = feature_label_unique),
      colour = "#5D337E", linewidth = 0.72
    ) +
    geom_point(
      aes(shape = support, fill = support),
      size = 2.85, colour = "#3E2A4C", stroke = 0.72
    ) +
    facet_grid(. ~ contrast, scales = "fixed") +
    scale_shape_manual(values = shape_values, drop = FALSE) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_x_continuous(
      limits = c(-max_ci, max_ci),
      breaks = pretty_breaks(n = 5),
      expand = expansion(mult = c(0.035, 0.035))
    ) +
    labs(
      x = "CLR-LM coefficient (95% CI)",
      y = y_title,
      shape = "Statistical support",
      fill = "Statistical support"
    ) +
    theme_editorial(base_size = 10.1) +
    theme(
      legend.position = if (show_legend) "right" else "none",
      legend.key.height = unit(0.65, "lines"),
      axis.text.y = element_text(size = 8.8),
      strip.text = element_text(size = 9.3)
    )
}

p_c <- forest_plot(humann_da$display, "HUMAnN pathways", TRUE)
p_d <- forest_plot(agora_da$display, "AGORA2 subsystems", FALSE)

# -------------------------------------------------------------------------
# Source data and provenance
# -------------------------------------------------------------------------
write_tsv(alpha, file.path(SRC_DIR, "Figure_S03_panelA_HUMAnN_pathway_Shannon.tsv"))
write_tsv(alpha_stats, file.path(SRC_DIR, "Figure_S03_panelA_statistics.tsv"))
write_tsv(
  pcoa[, .(sample_id, PCoA1, PCoA2, exposure_group, elimination_class, source)],
  file.path(SRC_DIR, "Figure_S03_panelB_HUMAnN_PCoA.tsv")
)
write_tsv(humann_da$full, file.path(SRC_DIR, "Figure_S03_panelC_HUMAnN_pathways_full.tsv"))
write_tsv(humann_da$display, file.path(SRC_DIR, "Figure_S03_panelC_HUMAnN_pathways_display.tsv"))
write_tsv(humann_da$ranking, file.path(SRC_DIR, "Figure_S03_panelC_HUMAnN_pathways_ranking.tsv"))
write_tsv(agora_da$full, file.path(SRC_DIR, "Figure_S03_panelD_AGORA2_subsystems_full.tsv"))
write_tsv(agora_da$display, file.path(SRC_DIR, "Figure_S03_panelD_AGORA2_subsystems_display.tsv"))
write_tsv(agora_da$ranking, file.path(SRC_DIR, "Figure_S03_panelD_AGORA2_subsystems_ranking.tsv"))

panel_design <- data.table(
  panel = c("A", "B", "C", "D", "E"),
  content = c(
    "HUMAnN pathway Shannon diversity",
    "Locked HUMAnN Bray-Curtis PCoA",
    "Corrected HUMAnN pathway CLR-LM associations",
    "Corrected AGORA2 subsystem CLR-LM associations",
    "Functional concordance"
  ),
  status = c("INCLUDED", "INCLUDED", "INCLUDED", "INCLUDED", "OMITTED"),
  reason = c(
    "Derived from matched HUMAnN pathway matrix; taxonomic alpha table not used",
    "Locked row-level PCoA source, n=65",
    "Authoritative corrected k_LPFOS_tertile model",
    "Authoritative corrected k_LPFOS_tertile model",
    "No defensible direct HUMAnN-AGORA2 concordance source"
  )
)
write_tsv(panel_design, file.path(SRC_DIR, "Figure_S03_panel_design.tsv"))

provenance <- data.table(
  role = c(
    "panel_A_matrix", "panel_A_group_map", "panel_B",
    "panel_C", "panel_D", "panel_E"
  ),
  source_path = c(
    PATH_MATRIX, PATH_PCOA, PATH_PCOA,
    PATH_HUMANN_DA, PATH_AGORA_DA, ""
  ),
  source_type = c(
    "Matched HUMAnN pathway abundance matrix",
    "Locked HUMAnN PCoA cohort metadata",
    "Locked HUMAnN PCoA coordinates",
    "Corrected CLR-LM output",
    "Corrected CLR-LM output",
    "OMITTED"
  ),
  rows_or_samples = c(
    length(sample_cols), nrow(pcoa), nrow(pcoa),
    nrow(humann_da$full), nrow(agora_da$full), 0
  )
)
write_tsv(provenance, file.path(MAN_DIR, "provenance.tsv"))

# -------------------------------------------------------------------------
# Editorial layout
# -------------------------------------------------------------------------
top <- p_a + p_b + plot_layout(widths = c(1.0, 1.18))
bottom <- p_c / p_d + plot_layout(heights = c(1.15, 1.0))
figure <- top / bottom +
  plot_layout(heights = c(0.82, 2.25)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 19, colour = "#171419"),
      plot.tag.position = c(0, 1)
    )
  )

pdf_path <- file.path(FIG_DIR, "Figure_S03_excellence_candidate.pdf")
png_path <- file.path(FIG_DIR, "Figure_S03_excellence_candidate.png")
tif_path <- file.path(FIG_DIR, "Figure_S03_excellence_candidate.tiff")

ggsave(pdf_path, figure, width = 16.2, height = 14.6, units = "in", device = cairo_pdf)
ggsave(png_path, figure, width = 16.2, height = 14.6, units = "in", dpi = 300, bg = "white")
ggsave(
  tif_path, figure, width = 16.2, height = 14.6, units = "in",
  dpi = 300, compression = "lzw", bg = "white"
)

legend_text <- paste(
  "Figure S3. Functional microbiome structure and associations.",
  "Panel A shows HUMAnN pathway Shannon diversity derived from the matched",
  "pathway abundance matrix. Panel B shows the locked HUMAnN Bray-Curtis",
  "ordination. Panels C and D show corrected covariate-adjusted CLR-LM",
  "coefficients and 95% confidence intervals for k_LPFOS elimination tertiles",
  "in HUMAnN pathways and AGORA2 subsystems, respectively.",
  "HUMAnN pathway labels are shortened for readability by removing trailing",
  "context qualifiers. Only statistical-support categories present in the",
  "displayed results are shown in the legend."
)
writeLines(legend_text, file.path(LEG_DIR, "Figure_S03.txt"))

source_files <- list.files(SRC_DIR, full.names = TRUE)
manifest <- data.table(
  file = basename(source_files),
  bytes = file.info(source_files)$size,
  sha256 = vapply(source_files, sha256_file, character(1))
)
write_tsv(manifest, file.path(MAN_DIR, "source_data_manifest.tsv"))

visual_qc <- data.table(
  check = "MANUAL_VISUAL_QC",
  status = "PENDING",
  detail = paste(
    "Confirm editorial balance, A-D labels, functional rather than taxonomic",
    "panel A, readable pathway/subsystem labels, unclipped confidence intervals,",
    "legends, axes and purple palette."
  )
)
write_tsv(visual_qc, file.path(MAN_DIR, "visual_qc.tsv"))

cat("FUNCTIONAL_ALPHA_SOURCE=PASS source=humann_pathways_sample_matched.tsv n=", nrow(alpha), "\n", sep = "")
cat("HUMANN_PCOA_SOURCE=PASS n=", nrow(pcoa), "\n", sep = "")
cat("HUMANN_CLRLM_SOURCE=PASS target_rows=", nrow(humann_da$full), "\n", sep = "")
cat("AGORA2_CLRLM_SOURCE=PASS target_rows=", nrow(agora_da$full), "\n", sep = "")
cat("PANEL_E_DECISION=OMITTED_NO_DEFENSIBLE_DIRECT_SOURCE\n")
cat("FIGURE_S03_PANEL_DESIGN=A_FUNCTIONAL_ALPHA_B_PCOA_C_HUMANN_D_AGORA2\n")
cat("HUMANN_LABEL_QUALIFIER_STRIP=PASS\n")
cat("SUPPORT_LEGEND_DYNAMIC=PASS levels=", paste(support_levels_to_show, collapse = ";"), "\n", sep = "")
cat("GENERATOR_EXECUTION_PASS\n")
