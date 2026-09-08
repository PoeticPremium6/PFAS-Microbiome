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
  stop("Usage: 41A_build_Figure_S05_supplementary.R ROOT STAGE")
}

ROOT <- normalizePath(args[[1]], mustWork = TRUE)
STAGE <- normalizePath(args[[2]], mustWork = FALSE)

BUILD <- file.path(STAGE, "figure_builds", "Figure_S05")
FIG_DIR <- file.path(BUILD, "figures")
SRC_DIR <- file.path(STAGE, "source_data", "Figure_S05")
MAN_DIR <- file.path(STAGE, "manifests", "Figure_S05")
LEG_DIR <- file.path(STAGE, "legends")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MAN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LEG_DIR, recursive = TRUE, showWarnings = FALSE)

LOCKED <- file.path(
  ROOT,
  "11_figures_LOCKED/figure_set_locked_20260714_090127/source_data"
)

PATH_PERFORMANCE <- file.path(
  LOCKED,
  "Figure_5A_RF50_performance_by_feature_set_v3.tsv"
)
PATH_OOF <- file.path(
  ROOT,
  "12_machine_learning_candidate_prioritization/phase12_oof_predictions.tsv"
)
PATH_STABILITY <- file.path(
  LOCKED,
  "Figure_5C_RF50_stable_features_v8.tsv"
)
PATH_EVIDENCE <- file.path(
  LOCKED,
  "Figure_5E_candidate_evidence_map_v5.tsv"
)
PATH_NAMED <- file.path(
  LOCKED,
  "Figure_5D_top_candidates_per_family_v8.tsv"
)

required <- c(
  PATH_PERFORMANCE,
  PATH_OOF,
  PATH_STABILITY,
  PATH_EVIDENCE,
  PATH_NAMED
)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop(
    "Missing authoritative Figure S5 input(s): ",
    paste(missing, collapse = "; ")
  )
}

purple_dark <- "#4B1F6F"
purple_mid <- "#70409B"
purple_light <- "#A985CC"
purple_pale <- "#DDD3EC"
ink <- "#302B33"

support_fill <- c(
  "FDR q<0.05" = purple_dark,
  "Nominal p<0.05" = purple_mid,
  "No nominal support" = "white"
)
support_shape <- c(
  "FDR q<0.05" = 21,
  "Nominal p<0.05" = 22,
  "No nominal support" = 23
)

theme_editorial <- function(base_size = 12.0) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_line(
        colour = "#E5DDED",
        linewidth = 0.44
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        colour = "#4A4650",
        fill = NA,
        linewidth = 0.64
      ),
      axis.title = element_text(
        face = "bold",
        colour = "#242126",
        size = rel(1.08),
        margin = margin(4, 4, 4, 4)
      ),
      axis.text = element_text(
        colour = "#3D3941",
        size = rel(0.95)
      ),
      strip.background = element_rect(
        fill = "#EEE7F4",
        colour = "#B9A7C9",
        linewidth = 0.56
      ),
      strip.text = element_text(
        face = "bold",
        colour = "#332C38",
        size = rel(1.02),
        margin = margin(5, 6, 5, 6)
      ),
      legend.title = element_text(
        face = "bold",
        size = rel(1.00)
      ),
      legend.text = element_text(
        colour = "#3D3941",
        size = rel(0.93)
      ),
      legend.key = element_blank(),
      plot.margin = margin(9, 10, 9, 10)
    )
}

first_existing <- function(candidates, available) {
  hit <- candidates[candidates %in% available]
  if (!length(hit)) return(NA_character_)
  hit[[1]]
}

canonical_feature_set <- function(x) {
  raw <- as.character(x)
  y <- tolower(raw)
  out <- rep(NA_character_, length(y))

  out[grepl("agora2.*exchange", y)] <- "AGORA2 exchanges"
  out[
    grepl("agora2.*model.*feature|agora2_model_features", y) &
    is.na(out)
  ] <- "AGORA2 model features"
  out[
    grepl("agora2.*reaction", y) &
    is.na(out)
  ] <- "AGORA2 reactions"
  out[
    grepl("agora2.*subsystem", y) &
    is.na(out)
  ] <- "AGORA2 subsystems"
  out[
    grepl("combined.*capped", y) &
    is.na(out)
  ] <- "Combined capped"
  out[
    grepl("humann.*pathway", y) &
    is.na(out)
  ] <- "HUMAnN pathways"
  out[
    grepl("mag", y) &
    is.na(out)
  ] <- "MAGs"
  out[
    grepl("taxonomy.*genus|metaphlan.*genus", y) &
    is.na(out)
  ] <- "MetaPhlAn genus"
  out[
    grepl("taxonomy.*species|metaphlan.*species", y) &
    is.na(out)
  ] <- "MetaPhlAn species"

  fallback <- gsub("_", " ", raw, fixed = TRUE)
  fallback <- gsub("[[:space:]]+", " ", fallback)
  fallback <- trimws(fallback)
  out[is.na(out)] <- fallback[is.na(out)]
  out
}

humanize_label <- function(x) {
  y <- as.character(x)
  y <- sub("\\|.*$", "", y)
  y <- gsub("_", " ", y, fixed = TRUE)
  y <- gsub("&alpha;", "alpha", y, fixed = TRUE)
  y <- gsub("&beta;", "beta", y, fixed = TRUE)
  y <- gsub("&gamma;", "gamma", y, fixed = TRUE)
  y <- gsub("[[:space:]]+", " ", y)
  trimws(y)
}

wrap_label <- function(x, width = 34L) {
  vapply(
    as.character(x),
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1)
  )
}

format_metric <- function(x) {
  y <- tolower(as.character(x))
  out <- rep(NA_character_, length(y))
  out[grepl("spearman|rho", y)] <- "Spearman"
  out[grepl("r2|r_squared|rsquared", y)] <- "R²"
  out
}

write_tsv <- function(dt, path) {
  export <- copy(as.data.table(dt))

  for (column in names(export)) {
    values <- export[[column]]

    if (is.factor(values)) {
      values <- as.character(values)
    }

    if (is.character(values)) {
      values <- gsub(
        "[\\r\\n\\t]+",
        " ",
        values,
        perl = TRUE
      )
      values <- gsub(
        "[[:space:]]+",
        " ",
        values
      )
      values <- trimws(values)
      export[[column]] <- values
    }
  }

  fwrite(
    export,
    path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
}

sha256_file <- function(path) {
  output <- suppressWarnings(
    system2(
      "sha256sum",
      shQuote(path),
      stdout = TRUE,
      stderr = TRUE
    )
  )
  if (!length(output)) return(NA_character_)
  strsplit(output[[1]], "[[:space:]]+")[[1]][1]
}

# -------------------------------------------------------------------------
# Panel A: cross-validated performance
# -------------------------------------------------------------------------
performance_raw <- fread(PATH_PERFORMANCE)
required_performance <- c("feature_set", "metric", "value")
if (!all(required_performance %in% names(performance_raw))) {
  stop(
    "Unexpected performance schema: ",
    paste(names(performance_raw), collapse = ",")
  )
}

performance <- performance_raw[
  ,
  .(
    feature_set_raw = as.character(feature_set),
    feature_set = canonical_feature_set(feature_set),
    metric = format_metric(metric),
    value = as.numeric(value)
  )
]
performance <- performance[
  !is.na(metric) &
  nzchar(feature_set) &
  is.finite(value)
]

setorder(performance, feature_set, metric)
performance <- performance[, .SD[1L], by = .(feature_set, metric)]
if (nrow(performance) != 18L) {
  stop(
    "Expected 18 RF50 performance rows; observed ",
    nrow(performance)
  )
}
if (uniqueN(performance$feature_set) != 9L) {
  stop(
    "Expected 9 RF50 feature sets; observed ",
    uniqueN(performance$feature_set)
  )
}

performance_wide <- dcast(
  performance,
  feature_set ~ metric,
  value.var = "value"
)
if (!all(c("R²", "Spearman") %in% names(performance_wide))) {
  stop("Performance table lacks R² or Spearman")
}

setorder(performance_wide, -Spearman, -`R²`, feature_set)
feature_order <- performance_wide$feature_set
performance[
  ,
  feature_set := factor(
    feature_set,
    levels = rev(feature_order)
  )
]

performance[
  ,
  text_hjust := fifelse(value >= 0, -0.20, 1.20)
]

performance_min <- min(performance$value, na.rm = TRUE)
performance_max <- max(performance$value, na.rm = TRUE)
performance_pad <- max(
  0.035,
  0.10 * (performance_max - performance_min)
)

p_a <- ggplot(
  performance,
  aes(
    x = value,
    y = feature_set
  )
) +
  geom_vline(
    xintercept = 0,
    colour = "#827987",
    linewidth = 0.74,
    linetype = "dashed"
  ) +
  geom_segment(
    data = performance_wide[
      ,
      feature_set := factor(
        feature_set,
        levels = rev(feature_order)
      )
    ],
    aes(
      x = `R²`,
      xend = Spearman,
      y = feature_set,
      yend = feature_set
    ),
    inherit.aes = FALSE,
    colour = "#CDBFD8",
    linewidth = 1.18
  ) +
  geom_point(
    aes(
      shape = metric,
      fill = metric
    ),
    size = 3.72,
    colour = "#3E3445",
    stroke = 0.78
  ) +
  geom_label(
    data = performance[metric == "R²"],
    aes(
      label = sprintf("%.2f", value)
    ),
    nudge_x = -0.012,
    nudge_y = 0.18,
    hjust = 0.5,
    size = 3.18,
    colour = "#3D3740",
    fill = alpha("white", 0.94),
    label.padding = unit(0.11, "lines"),
    label.r = unit(0.08, "lines"),
    linewidth = 0.28,
    show.legend = FALSE
  ) +
  geom_label(
    data = performance[metric == "Spearman"],
    aes(
      label = sprintf("%.2f", value)
    ),
    nudge_x = 0.012,
    nudge_y = -0.18,
    hjust = 0.5,
    size = 3.18,
    colour = "#3D3740",
    fill = alpha("white", 0.94),
    label.padding = unit(0.11, "lines"),
    label.r = unit(0.08, "lines"),
    linewidth = 0.28,
    show.legend = FALSE
  ) +
  scale_shape_manual(
    values = c("R²" = 21, "Spearman" = 24)
  ) +
  scale_fill_manual(
    values = c("R²" = purple_mid, "Spearman" = purple_light)
  ) +
  scale_y_discrete(
    expand = expansion(add = 0.62)
  ) +
  scale_x_continuous(
    limits = c(
      performance_min - performance_pad * 1.32,
      performance_max + performance_pad * 1.32
    ),
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0.015, 0.015))
  ) +
  labs(
    x = "Cross-validated performance",
    y = "Feature set",
    shape = "Metric",
    fill = "Metric"
  ) +
  theme_editorial(base_size = 12.1) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.key.width = unit(1.05, "lines"),
    legend.key.height = unit(0.76, "lines"),
    legend.spacing.x = unit(5, "pt"),
    axis.text.y = element_text(size = 10.1),
    axis.text.x = element_text(size = 9.8)
  )

# -------------------------------------------------------------------------
# Panel B: participant-level aggregated OOF calibration
# -------------------------------------------------------------------------
oof_raw <- fread(PATH_OOF)
required_oof <- c(
  "feature_set",
  "prediction_repeat",
  "fold",
  "sample",
  "observed",
  "predicted"
)
if (!all(required_oof %in% names(oof_raw))) {
  stop(
    "Unexpected OOF schema: ",
    paste(names(oof_raw), collapse = ",")
  )
}

if ("target" %in% names(oof_raw)) {
  target_counts <- oof_raw[, .N, by = target][order(-N)]
  target_keep <- target_counts$target[[1]]
  oof_raw <- oof_raw[target == target_keep]
} else {
  target_keep <- "unspecified"
}

oof <- oof_raw[
  ,
  .(
    feature_set_raw = as.character(feature_set),
    feature_set = canonical_feature_set(feature_set),
    sample_id = as.character(sample),
    observed = as.numeric(observed),
    predicted = as.numeric(predicted),
    prediction_repeat = as.character(prediction_repeat),
    fold = as.character(fold)
  )
]
oof <- oof[
  nzchar(feature_set) &
  nzchar(sample_id) &
  is.finite(observed) &
  is.finite(predicted)
]

observed_consistency <- oof[
  ,
  .(
    unique_observed = uniqueN(round(observed, 12))
  ),
  by = .(feature_set, sample_id)
]
if (observed_consistency[unique_observed != 1L, .N] > 0L) {
  stop(
    "Observed outcome differs within ",
    observed_consistency[unique_observed != 1L, .N],
    " feature-set/sample combinations"
  )
}

oof_aggregated <- oof[
  ,
  .(
    observed = observed[[1]],
    predicted_mean = mean(predicted),
    predicted_median = median(predicted),
    predicted_sd = sd(predicted),
    predicted_q025 = as.numeric(
      quantile(predicted, 0.025, names = FALSE)
    ),
    predicted_q975 = as.numeric(
      quantile(predicted, 0.975, names = FALSE)
    ),
    n_predictions = .N,
    n_repeats = uniqueN(prediction_repeat),
    n_folds = uniqueN(fold)
  ),
  by = .(feature_set, sample_id)
]

if (uniqueN(oof_aggregated$sample_id) != 47L) {
  stop(
    "Expected 47 unique elimination participants; observed ",
    uniqueN(oof_aggregated$sample_id)
  )
}
if (uniqueN(oof_aggregated$feature_set) != 9L) {
  stop(
    "Expected 9 OOF feature sets; observed ",
    uniqueN(oof_aggregated$feature_set)
  )
}
if (nrow(oof_aggregated) != 423L) {
  stop(
    "Expected 423 participant-feature-set aggregates; observed ",
    nrow(oof_aggregated)
  )
}
if (
  min(oof_aggregated$n_predictions) != 50L ||
  max(oof_aggregated$n_predictions) != 50L
) {
  stop(
    "Expected exactly 50 predictions per participant-feature-set; range ",
    min(oof_aggregated$n_predictions),
    "-",
    max(oof_aggregated$n_predictions)
  )
}

oof_summary <- oof_aggregated[
  ,
  {
    sse <- sum((observed - predicted_mean)^2)
    sst <- sum((observed - mean(observed))^2)
    .(
      n_participants = .N,
      spearman_rho = suppressWarnings(
        cor(
          observed,
          predicted_mean,
          method = "spearman",
          use = "complete.obs"
        )
      ),
      r_squared = if (sst > 0) 1 - sse / sst else NA_real_,
      rmse = sqrt(mean((observed - predicted_mean)^2)),
      mean_repeated_predictions = mean(n_predictions)
    )
  },
  by = feature_set
]

selection <- merge(
  performance_wide,
  oof_summary,
  by = "feature_set",
  all = FALSE
)
setorder(selection, -Spearman, -`R²`, feature_set)
selection <- head(selection, 4L)

if (nrow(selection) != 4L) {
  stop(
    "Could not select four matched performance/OOF feature sets"
  )
}

selected_feature_sets <- selection$feature_set
oof_selected <- oof_aggregated[
  feature_set %in% selected_feature_sets
]
oof_selected[
  ,
  feature_set := factor(
    feature_set,
    levels = selected_feature_sets
  )
]

annotation_b <- merge(
  selection[
    ,
    .(
      feature_set,
      performance_spearman = Spearman,
      performance_r2 = `R²`
    )
  ],
  oof_summary,
  by = "feature_set",
  all.x = TRUE
)
annotation_b[
  ,
  label := sprintf(
    "n=47; rho=%.2f\nR²=%.2f; RMSE=%.2f",
    spearman_rho,
    r_squared,
    rmse
  )
]
annotation_b[
  ,
  feature_set := factor(
    feature_set,
    levels = selected_feature_sets
  )
]

calibration_min <- min(
  c(
    oof_selected$observed,
    oof_selected$predicted_q025
  ),
  na.rm = TRUE
)
calibration_max <- max(
  c(
    oof_selected$observed,
    oof_selected$predicted_q975
  ),
  na.rm = TRUE
)
calibration_pad <- 0.045 * (
  calibration_max - calibration_min
)
calibration_limits <- c(
  calibration_min - calibration_pad,
  calibration_max + calibration_pad
)

p_b <- ggplot(
  oof_selected,
  aes(
    x = observed,
    y = predicted_mean
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "#746B79",
    linewidth = 0.76,
    linetype = "dashed"
  ) +
  geom_linerange(
    aes(
      ymin = predicted_q025,
      ymax = predicted_q975
    ),
    colour = purple_light,
    linewidth = 0.48,
    alpha = 0.38
  ) +
  geom_point(
    shape = 21,
    size = 2.72,
    stroke = 0.58,
    fill = purple_mid,
    colour = "#3E2A4C",
    alpha = 0.84
  ) +
  geom_label(
    data = annotation_b,
    aes(
      x = -Inf,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = -0.05,
    vjust = 1.10,
    size = 3.22,
    lineheight = 0.98,
    fontface = "italic",
    colour = "#39323D",
    fill = alpha("white", 0.86),
    label.padding = unit(0.10, "lines"),
    label.r = unit(0.05, "lines"),
    linewidth = 0,
    label.size = 0
  ) +
  facet_wrap(
    ~ feature_set,
    ncol = 2
  ) +
  scale_x_continuous(
    limits = calibration_limits,
    breaks = pretty_breaks(n = 4)
  ) +
  scale_y_continuous(
    limits = calibration_limits,
    breaks = pretty_breaks(n = 4)
  ) +
  coord_equal() +
  labs(
    x = "Observed elimination score",
    y = "Mean out-of-fold prediction"
  ) +
  theme_editorial(base_size = 11.5) +
  theme(
    strip.text = element_text(size = 10.6),
    axis.text = element_text(size = 9.45),
    panel.spacing = unit(8, "pt"),
    plot.margin = margin(9, 9, 9, 9)
  )

# -------------------------------------------------------------------------
# Panel C: stable selected features
# -------------------------------------------------------------------------
stability_raw <- fread(PATH_STABILITY)
stability_col <- first_existing(
  c(
    "stability_plot",
    "stability",
    "mean_normalized_importance"
  ),
  names(stability_raw)
)
stability_label_col <- first_existing(
  c(
    "label_clean",
    "label",
    "label_raw",
    "feature"
  ),
  names(stability_raw)
)
stability_family_col <- first_existing(
  c(
    "feature_family_fixed",
    "feature_family",
    "family"
  ),
  names(stability_raw)
)

if (
  is.na(stability_col) ||
  is.na(stability_label_col) ||
  is.na(stability_family_col)
) {
  stop(
    "Unexpected stability schema: ",
    paste(names(stability_raw), collapse = ",")
  )
}

stability <- stability_raw[
  ,
  .(
    feature = as.character(feature),
    feature_label = humanize_label(
      get(stability_label_col)
    ),
    feature_family = canonical_feature_set(
      get(stability_family_col)
    ),
    stability = as.numeric(
      get(stability_col)
    )
  )
]
stability <- stability[
  nzchar(feature_label) &
  nzchar(feature_family) &
  is.finite(stability)
]
setorder(stability, -stability, feature_label)
stability <- stability[, .SD[1L], by = feature]

if (nrow(stability) < 15L) {
  stop(
    "Stable-feature source unexpectedly small: ",
    nrow(stability)
  )
}
stability <- head(stability, 18L)

stability[
  ,
  feature_label_plot := wrap_label(
    feature_label,
    width = 35L
  )
]
stability[
  ,
  feature_label_plot := factor(
    feature_label_plot,
    levels = rev(feature_label_plot)
  )
]

stability_families <- unique(stability$feature_family)
family_palette <- setNames(
  colorRampPalette(
    c(
      "#4B1F6F",
      "#70409B",
      "#9369B6",
      "#B596CD",
      "#D6C8E5"
    )
  )(length(stability_families)),
  stability_families
)

stability_max <- max(stability$stability, na.rm = TRUE)

p_c <- ggplot(
  stability,
  aes(
    x = stability,
    y = feature_label_plot
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = stability,
      yend = feature_label_plot
    ),
    colour = "#CFC3DA",
    linewidth = 1.12
  ) +
  geom_point(
    aes(
      fill = feature_family
    ),
    shape = 21,
    size = 3.88,
    colour = "#3E2A4C",
    stroke = 0.72
  ) +
  geom_label(
    aes(
      label = sprintf("%.3f", stability)
    ),
    nudge_x = stability_max * 0.024,
    hjust = 0,
    size = 3.08,
    colour = "#403943",
    fill = alpha("white", 0.94),
    label.padding = unit(0.10, "lines"),
    label.r = unit(0.07, "lines"),
    linewidth = 0.26,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = family_palette
  ) +
  scale_x_continuous(
    limits = c(0, stability_max * 1.28),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "RF50 stability / normalized importance",
    y = "Stable selected features",
    fill = "Feature family"
  ) +
  theme_editorial(base_size = 11.35) +
  theme(
    axis.text.y = element_text(
      size = 9.05,
      lineheight = 0.96
    ),
    axis.text.x = element_text(size = 9.5),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.key.width = unit(0.78, "lines"),
    legend.key.height = unit(0.72, "lines"),
    legend.spacing.x = unit(4, "pt"),
    plot.margin = margin(9, 10, 9, 10)
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(size = 3.8)
    )
  )

# -------------------------------------------------------------------------
# Panel D: candidate-prioritization evidence
# -------------------------------------------------------------------------
evidence_raw <- fread(PATH_EVIDENCE)
evidence_label_col <- first_existing(
  c(
    "label_plot",
    "label_raw",
    "feature"
  ),
  names(evidence_raw)
)
evidence_family_col <- first_existing(
  c(
    "feature_family",
    "feature_family_fixed",
    "family"
  ),
  names(evidence_raw)
)
evidence_score_col <- first_existing(
  c(
    "score_plot",
    "candidate_score",
    "score"
  ),
  names(evidence_raw)
)
evidence_stability_col <- first_existing(
  c(
    "stability_plot",
    "stability"
  ),
  names(evidence_raw)
)

if (
  is.na(evidence_label_col) ||
  is.na(evidence_family_col) ||
  is.na(evidence_score_col) ||
  is.na(evidence_stability_col)
) {
  stop(
    "Unexpected candidate-evidence schema: ",
    paste(names(evidence_raw), collapse = ",")
  )
}

evidence <- evidence_raw[
  ,
  .(
    feature = as.character(feature),
    feature_label = humanize_label(
      get(evidence_label_col)
    ),
    feature_family = canonical_feature_set(
      get(evidence_family_col)
    ),
    candidate_score = as.numeric(
      get(evidence_score_col)
    ),
    stability = as.numeric(
      get(evidence_stability_col)
    ),
    q_value = if ("q" %in% names(evidence_raw)) {
      as.numeric(q)
    } else {
      NA_real_
    },
    p_value = if ("p" %in% names(evidence_raw)) {
      as.numeric(p)
    } else {
      NA_real_
    },
    positive_flag = if ("positive_flag" %in% names(evidence_raw)) {
      as.character(positive_flag)
    } else if ("positive" %in% names(evidence_raw)) {
      as.character(positive)
    } else {
      NA_character_
    }
  )
]
evidence <- evidence[
  nzchar(feature_label) &
  nzchar(feature_family) &
  is.finite(candidate_score) &
  is.finite(stability)
]

evidence[
  ,
  support := fifelse(
    is.finite(q_value) & q_value < 0.05,
    "FDR q<0.05",
    fifelse(
      is.finite(p_value) & p_value < 0.05,
      "Nominal p<0.05",
      "No nominal support"
    )
  )
]
evidence[
  ,
  support_rank := fifelse(
    support == "FDR q<0.05",
    1L,
    fifelse(
      support == "Nominal p<0.05",
      2L,
      3L
    )
  )
]
setorder(
  evidence,
  -candidate_score,
  support_rank,
  -stability,
  feature_label
)
evidence <- evidence[, .SD[1L], by = feature]

if (nrow(evidence) < 10L) {
  stop(
    "Candidate-evidence source unexpectedly small: ",
    nrow(evidence)
  )
}

evidence_display <- head(evidence, 12L)
evidence_display[
  ,
  feature_label_plot := wrap_label(
    feature_label,
    width = 33L
  )
]
evidence_display[
  ,
  feature_label_plot := factor(
    feature_label_plot,
    levels = rev(feature_label_plot)
  )
]
support_order <- c(
  "FDR q<0.05",
  "Nominal p<0.05",
  "No nominal support"
)
support_levels_present <- support_order[
  support_order %chin% unique(evidence_display$support)
]
if (!length(support_levels_present)) {
  stop("No valid association-support level in displayed candidates")
}

evidence_display[
  ,
  support := factor(
    as.character(support),
    levels = support_levels_present
  )
]

support_fill_present <- support_fill[support_levels_present]
support_shape_present <- support_shape[support_levels_present]

support_level_manifest <- evidence_display[
  ,
  .N,
  by = support
]
support_level_manifest[
  ,
  support := as.character(support)
]
support_level_manifest[
  ,
  display_order := match(
    support,
    support_levels_present
  )
]
setorder(
  support_level_manifest,
  display_order
)

candidate_max <- max(
  evidence_display$candidate_score,
  na.rm = TRUE
)
size_legend_needed <- (
  uniqueN(
    round(
      evidence_display$stability,
      10
    )
  ) > 1L
)
candidate_size_range <- if (
  size_legend_needed
) {
  c(3.25, 5.25)
} else {
  c(4.15, 4.15)
}

p_d <- ggplot(
  evidence_display,
  aes(
    x = candidate_score,
    y = feature_label_plot
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = candidate_score,
      yend = feature_label_plot
    ),
    colour = "#CFC3DA",
    linewidth = 1.12
  ) +
  geom_point(
    aes(
      shape = support,
      fill = support,
      size = stability
    ),
    colour = "#3E2A4C",
    stroke = 0.76
  ) +
  geom_label(
    aes(
      label = sprintf("%.2f", candidate_score)
    ),
    nudge_x = candidate_max * 0.022,
    hjust = 0,
    size = 3.12,
    colour = "#403943",
    fill = alpha("white", 0.94),
    label.padding = unit(0.10, "lines"),
    label.r = unit(0.07, "lines"),
    linewidth = 0.26,
    show.legend = FALSE
  ) +
  scale_shape_manual(
    values = support_shape_present,
    drop = TRUE
  ) +
  scale_fill_manual(
    values = support_fill_present,
    drop = TRUE
  ) +
  scale_size_continuous(
    range = candidate_size_range
  ) +
  scale_x_continuous(
    limits = c(0, candidate_max * 1.20),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Candidate-prioritization score",
    y = "Prioritized candidate features",
    shape = "Association support",
    fill = "Association support",
    size = "RF stability"
  ) +
  theme_editorial(base_size = 11.35) +
  theme(
    axis.text.y = element_text(
      size = 9.05,
      lineheight = 0.96
    ),
    axis.text.x = element_text(size = 9.5),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.key.width = unit(0.80, "lines"),
    legend.key.height = unit(0.72, "lines"),
    legend.spacing.x = unit(4, "pt"),
    plot.margin = margin(9, 10, 9, 10)
  ) +
  guides(
    shape = guide_legend(
      order = 1,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(size = 4.1)
    ),
    fill = guide_legend(
      order = 1,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(size = 4.1)
    ),
    size = if (
      size_legend_needed
    ) {
      guide_legend(
        order = 2,
        nrow = 1,
        byrow = TRUE
      )
    } else {
      "none"
    }
  )

# -------------------------------------------------------------------------
# Source data, statistics and provenance
# -------------------------------------------------------------------------
performance_summary <- merge(
  performance_wide,
  oof_summary,
  by = "feature_set",
  all.x = TRUE
)
setorder(performance_summary, -Spearman, -`R²`)

write_tsv(
  performance,
  file.path(
    SRC_DIR,
    "Figure_S05_panelA_RF50_performance_long.tsv"
  )
)
write_tsv(
  performance_summary,
  file.path(
    SRC_DIR,
    "Figure_S05_panelA_RF50_performance_summary.tsv"
  )
)
write_tsv(
  oof_aggregated,
  file.path(
    SRC_DIR,
    "Figure_S05_panelB_OOF_participant_aggregated_all.tsv"
  )
)
write_tsv(
  oof_selected,
  file.path(
    SRC_DIR,
    "Figure_S05_panelB_OOF_participant_aggregated_display.tsv"
  )
)
write_tsv(
  oof_summary,
  file.path(
    SRC_DIR,
    "Figure_S05_panelB_OOF_model_statistics.tsv"
  )
)
write_tsv(
  selection,
  file.path(
    SRC_DIR,
    "Figure_S05_panelB_selected_feature_sets.tsv"
  )
)
write_tsv(
  stability,
  file.path(
    SRC_DIR,
    "Figure_S05_panelC_stable_features.tsv"
  )
)
write_tsv(
  evidence,
  file.path(
    SRC_DIR,
    "Figure_S05_panelD_candidate_evidence_full.tsv"
  )
)
write_tsv(
  evidence_display,
  file.path(
    SRC_DIR,
    "Figure_S05_panelD_candidate_evidence_display.tsv"
  )
)
write_tsv(
  support_level_manifest,
  file.path(
    SRC_DIR,
    "Figure_S05_panelD_support_levels.tsv"
  )
)
write_tsv(
  fread(PATH_NAMED),
  file.path(
    SRC_DIR,
    "Figure_S05_named_candidates_reference.tsv"
  )
)

aggregation_audit <- data.table(
  metric = c(
    "raw_prediction_rows",
    "unique_participants",
    "feature_sets",
    "prediction_repeats",
    "aggregated_rows",
    "selected_feature_sets",
    "selected_display_rows",
    "minimum_predictions_per_participant_model",
    "maximum_predictions_per_participant_model",
    "target"
  ),
  value = c(
    nrow(oof),
    uniqueN(oof_aggregated$sample_id),
    uniqueN(oof_aggregated$feature_set),
    uniqueN(oof$prediction_repeat),
    nrow(oof_aggregated),
    length(selected_feature_sets),
    nrow(oof_selected),
    min(oof_aggregated$n_predictions),
    max(oof_aggregated$n_predictions),
    target_keep
  )
)
write_tsv(
  aggregation_audit,
  file.path(
    SRC_DIR,
    "Figure_S05_panelB_aggregation_audit.tsv"
  )
)

panel_design <- data.table(
  panel = c("A", "B", "C", "D"),
  content = c(
    "Cross-validated RF50 R-squared and Spearman performance",
    "Participant-level aggregated out-of-fold calibration",
    "Stable selected RF50 features",
    "Candidate-prioritization evidence ranking"
  ),
  status = c(
    "INCLUDED",
    "INCLUDED",
    "INCLUDED",
    "INCLUDED"
  ),
  reason = c(
    "Locked feature-set performance source; negative and near-zero values retained",
    "Twenty-one thousand one hundred fifty repeated prediction rows aggregated to 423 participant-model rows; displayed facets each use 47 unique participants",
    "Locked stable-feature source with readable labels and feature-family provenance",
    "Locked candidate evidence source; ranking is hypothesis prioritization rather than causal validation"
  )
)
write_tsv(
  panel_design,
  file.path(
    SRC_DIR,
    "Figure_S05_panel_design.tsv"
  )
)

provenance <- data.table(
  role = c(
    "panel_A",
    "panel_B",
    "panel_C",
    "panel_D",
    "named_candidate_reference"
  ),
  source_path = c(
    PATH_PERFORMANCE,
    PATH_OOF,
    PATH_STABILITY,
    PATH_EVIDENCE,
    PATH_NAMED
  ),
  source_type = c(
    "Locked RF50 performance",
    "Repeated out-of-fold prediction records aggregated within participant and feature set",
    "Locked stable RF50 features",
    "Locked candidate evidence map",
    "Locked optional named-candidate reference"
  ),
  rows_or_samples = c(
    nrow(performance),
    uniqueN(oof_aggregated$sample_id),
    nrow(stability),
    nrow(evidence),
    nrow(fread(PATH_NAMED))
  )
)
write_tsv(
  provenance,
  file.path(MAN_DIR, "provenance.tsv")
)

# -------------------------------------------------------------------------
# Editorial layout
# -------------------------------------------------------------------------
top <- p_a + p_b + plot_layout(
  widths = c(0.96, 1.16)
)
bottom <- p_c + p_d + plot_layout(
  widths = c(1.02, 1.00)
)
figure <- top / bottom +
  plot_layout(
    heights = c(0.98, 1.08)
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        face = "bold",
        size = 22,
        colour = "#171419"
      ),
      plot.tag.position = c(0, 1)
    )
  )

pdf_path <- file.path(
  FIG_DIR,
  "Figure_S05_excellence_candidate.pdf"
)
png_path <- file.path(
  FIG_DIR,
  "Figure_S05_excellence_candidate.png"
)
tiff_path <- file.path(
  FIG_DIR,
  "Figure_S05_excellence_candidate.tiff"
)

ggsave(
  pdf_path,
  figure,
  width = 17.2,
  height = 13.8,
  units = "in",
  device = cairo_pdf
)
ggsave(
  png_path,
  figure,
  width = 17.2,
  height = 13.8,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  tiff_path,
  figure,
  width = 17.2,
  height = 13.8,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

legend_text <- paste(
  "Figure S5. Cross-validated random-forest performance and candidate",
  "prioritization diagnostics. Panel A shows cross-validated R-squared and",
  "Spearman performance across nine feature sets, retaining negative and",
  "near-zero values. Panel B shows participant-level out-of-fold calibration",
  "for the four feature sets with the highest locked Spearman performance.",
  "The 50 repeated predictions were aggregated within each of 47 participants",
  "and feature sets; points show mean predictions and vertical intervals show",
  "the 2.5th-97.5th percentiles of repeated predictions. Panel C shows stable",
  "features selected across repeated RF50 fits. Panel D shows ranked candidate",
  "evidence, with point size indicating RF stability and symbol indicating",
  "association support. Only association-support categories represented among",
  "the displayed candidates are included in the panel D legend. Negative or",
  "near-zero cross-validated performance indicates limited predictive utility.",
  "Candidate rankings are hypothesis priorities and do not establish causality",
  "or PFAS degradation."
)
writeLines(
  legend_text,
  file.path(LEG_DIR, "Figure_S05.txt")
)

source_files <- list.files(
  SRC_DIR,
  full.names = TRUE
)
source_manifest <- data.table(
  file = basename(source_files),
  bytes = file.info(source_files)$size,
  sha256 = vapply(
    source_files,
    sha256_file,
    character(1)
  )
)
write_tsv(
  source_manifest,
  file.path(
    MAN_DIR,
    "source_data_manifest.tsv"
  )
)

visual_qc <- data.table(
  check = "MANUAL_VISUAL_QC",
  status = "PENDING",
  detail = paste(
    "Confirm A-D panel labels; enlarged axis, facet, legend and annotation",
    "fonts; non-overlapping high-contrast numeric callouts in A, C and D;",
    "visible zero-reference line and negative performance in A; n=47 unique",
    "participants in each B facet; restrained repeated-prediction intervals;",
    "stable-feature family legend in C; panel D legend restricted to support",
    "categories actually represented among displayed candidates; balanced",
    "spacing and publication-grade purple editorial styling."
  )
)
write_tsv(
  visual_qc,
  file.path(MAN_DIR, "visual_qc.tsv")
)

cat(
  "RF50_PERFORMANCE_SOURCE=PASS rows=",
  nrow(performance),
  " feature_sets=",
  uniqueN(performance$feature_set),
  "\n",
  sep = ""
)
cat(
  "OOF_SOURCE=PASS raw_rows=",
  nrow(oof),
  " unique_participants=",
  uniqueN(oof_aggregated$sample_id),
  " feature_sets=",
  uniqueN(oof_aggregated$feature_set),
  "\n",
  sep = ""
)
cat(
  "OOF_AGGREGATION=PASS participant_model_rows=",
  nrow(oof_aggregated),
  " predictions_per_participant_model=",
  min(oof_aggregated$n_predictions),
  "\n",
  sep = ""
)
cat(
  "OOF_DISPLAY=PASS selected_models=",
  length(selected_feature_sets),
  " display_rows=",
  nrow(oof_selected),
  " unique_n_per_model=47\n",
  sep = ""
)
cat(
  "RF50_STABLE_FEATURE_SOURCE=PASS rows=",
  nrow(stability),
  "\n",
  sep = ""
)
cat(
  "CANDIDATE_EVIDENCE_SOURCE=PASS full_rows=",
  nrow(evidence),
  " display_rows=",
  nrow(evidence_display),
  "\n",
  sep = ""
)
cat(
  "PSEUDOREPLICATION_GATE=PASS_AGGREGATED_BEFORE_PLOTTING\n"
)
cat(
  "PANEL_D_SUPPORT_LEVELS=",
  paste(
    support_levels_present,
    collapse = "|"
  ),
  "\n",
  sep = ""
)
cat(
  "EDITORIAL_READABILITY_REFINEMENT=PASS numeric_callouts_nonoverlap=TRUE dynamic_support_legend=TRUE\n"
)
cat(
  "FIGURE_S05_PANEL_DESIGN=A_PERFORMANCE_B_AGGREGATED_OOF_C_STABILITY_D_EVIDENCE\n"
)
cat("TSV_SERIALIZATION_SANITIZER=PASS factors_newlines_tabs=TRUE\n")
cat("GENERATOR_EXECUTION_PASS\n")
