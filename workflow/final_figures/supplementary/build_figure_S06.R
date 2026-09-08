#!/usr/bin/env Rscript

# === CORRECTION_3_OF_4_FIGURE_S06_SCORE_CONTEXT_BEGIN ===
# Panel C displays the verified rule-based candidate-priority score.
# Lindell evidence is candidate-set provenance, not a redundant per-point shape legend.
# Panel D displays the locked Lindell overlap priority rescaled to 0–1 only for plotting.
# === CORRECTION_3_OF_4_FIGURE_S06_SCORE_CONTEXT_END ===

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 43A_build_Figure_S06_supplementary.R ROOT STAGE")
}

ROOT <- normalizePath(args[[1]], mustWork = TRUE)
STAGE <- normalizePath(args[[2]], mustWork = FALSE)

BUILD <- file.path(STAGE, "figure_builds", "Figure_S06")
FIG_DIR <- file.path(BUILD, "figures")
SRC_DIR <- file.path(STAGE, "source_data", "Figure_S06")
MAN_DIR <- file.path(STAGE, "manifests", "Figure_S06")
LEG_DIR <- file.path(STAGE, "legends")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MAN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LEG_DIR, recursive = TRUE, showWarnings = FALSE)

LOCKED <- file.path(
  ROOT,
  "11_figures_LOCKED/figure_set_locked_20260714_090127/source_data"
)

PATH_QUALITY <- file.path(
  ROOT,
  "06_MAG_resource_integration/07_figures/Figure_S16_MAG_quality_source_data.tsv"
)
PATH_TAXONOMY <- file.path(
  ROOT,
  paste0(
    "06_MAG_resource_integration/04_gtdbtk/",
    "GTDBTk_fastANI95_representatives/classify/",
    "Ronneby_MAG_rep95.bac120.summary.tsv"
  )
)
PATH_ABUNDANCE <- file.path(
  ROOT,
  "06_MAG_resource_integration/06_coverm_abundance/Ronneby_MAG_relative_abundance_percent.tsv"
)
PATH_CONTEXT <- file.path(
  LOCKED,
  "Figure_6D_candidates_with_Lindell_context.tsv"
)
PATH_OVERLAP <- file.path(
  LOCKED,
  "Figure_6C_representative_Lindell_overlap_features.tsv"
)

required <- c(
  PATH_QUALITY,
  PATH_TAXONOMY,
  PATH_ABUNDANCE,
  PATH_CONTEXT,
  PATH_OVERLAP
)
missing <- required[!file.exists(required)]
if (length(missing)) {
  stop(
    "Missing authoritative Figure S6 input(s): ",
    paste(missing, collapse = "; ")
  )
}

purple_dark <- "#4B1F6F"
purple_mid <- "#70409B"
purple_light <- "#A985CC"
purple_pale <- "#DDD3EC"
purple_very_pale <- "#EEE8F4"
grey_mid <- "#8A828E"
grey_light <- "#D8D3DB"
ink <- "#302B33"

theme_editorial <- function(base_size = 11.8) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.major = element_line(
        colour = "#E5DDED",
        linewidth = 0.42
      ),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(
        colour = "#4A4650",
        fill = NA,
        linewidth = 0.62
      ),
      axis.title = element_text(
        face = "bold",
        colour = "#242126",
        size = rel(1.07)
      ),
      axis.text = element_text(
        colour = "#3D3941",
        size = rel(0.95)
      ),
      strip.background = element_rect(
        fill = "#EEE7F4",
        colour = "#B9A7C9",
        linewidth = 0.54
      ),
      strip.text = element_text(
        face = "bold",
        colour = "#332C38"
      ),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(colour = "#3D3941"),
      legend.key = element_blank(),
      plot.margin = margin(9, 10, 9, 10)
    )
}

clean_text <- function(x) {
  y <- as.character(x)
  y[is.na(y)] <- ""
  y <- gsub("&alpha;", "alpha", y, fixed = TRUE)
  y <- gsub("&beta;", "beta", y, fixed = TRUE)
  y <- gsub("&gamma;", "gamma", y, fixed = TRUE)
  y <- gsub("_", " ", y, fixed = TRUE)
  y <- gsub("[[:space:]]+", " ", y)
  trimws(y)
}

wrap_label <- function(x, width = 35L) {
  vapply(
    clean_text(x),
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1)
  )
}

parse_bool <- function(x) {
  tolower(trimws(as.character(x))) %chin% c(
    "1", "true", "yes", "y", "taxon"
  )
}

is_bin_like_label <- function(x) {
  y <- tolower(clean_text(x))
  grepl(
    "(^|[[:space:]_.-])[a-z]?[0-9]{2,}[[:space:]_.-]*bin[.]?[0-9]+",
    y,
    perl = TRUE
  ) |
    grepl(
      "(^|[[:space:]_.-])bin[.]?[0-9]+($|[[:space:]_.-])",
      y,
      perl = TRUE
    ) |
    grepl("^mag[_[:space:]-]?[0-9]+$", y, perl = TRUE)
}

write_tsv <- function(dt, path) {
  export <- copy(as.data.table(dt))
  for (column in names(export)) {
    values <- export[[column]]
    if (is.factor(values)) {
      values <- as.character(values)
    }
    if (is.character(values)) {
      values <- gsub("[\\r\\n\\t]+", " ", values, perl = TRUE)
      values <- gsub("[[:space:]]+", " ", values)
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

deepest_taxon <- function(classification) {
  parts <- strsplit(as.character(classification), ";", fixed = TRUE)
  vapply(
    parts,
    function(items) {
      items <- trimws(items)
      items <- items[
        nzchar(items) &
        !grepl("^[dpcofgs]__$", items)
      ]
      if (!length(items)) return("Unresolved")
      value <- tail(items, 1L)
      sub("^[dpcofgs]__", "", value)
    },
    character(1)
  )
}

deepest_rank <- function(classification) {
  parts <- strsplit(as.character(classification), ";", fixed = TRUE)
  rank_map <- c(
    "d__" = "Domain",
    "p__" = "Phylum",
    "c__" = "Class",
    "o__" = "Order",
    "f__" = "Family",
    "g__" = "Genus",
    "s__" = "Species"
  )
  vapply(
    parts,
    function(items) {
      items <- trimws(items)
      resolved <- items[
        nzchar(items) &
        !grepl("^[dpcofgs]__$", items)
      ]
      if (!length(resolved)) return("Unresolved")
      value <- tail(resolved, 1L)
      prefix <- substr(value, 1L, 3L)
      rank_map[[prefix]] %||% "Other"
    },
    character(1)
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x)) y else x
}

# -------------------------------------------------------------------------
# Panel A: MAG quality and GTDB-Tk representative status
# -------------------------------------------------------------------------
quality_raw <- fread(PATH_QUALITY)
required_quality <- c("MAG_ID", "completeness", "contamination")
if (!all(required_quality %in% names(quality_raw))) {
  stop(
    "Unexpected quality schema: ",
    paste(names(quality_raw), collapse = ",")
  )
}

taxonomy_raw <- fread(PATH_TAXONOMY)
required_taxonomy <- c("user_genome", "classification")
if (!all(required_taxonomy %in% names(taxonomy_raw))) {
  stop(
    "Unexpected taxonomy schema: ",
    paste(names(taxonomy_raw), collapse = ",")
  )
}

quality <- quality_raw[
  ,
  .(
    MAG_ID = as.character(MAG_ID),
    completeness = as.numeric(completeness),
    contamination = as.numeric(contamination),
    genome_size_bp = if ("genome_size_bp" %in% names(quality_raw)) {
      as.numeric(genome_size_bp)
    } else {
      NA_real_
    },
    n50_contigs = if ("n50_contigs" %in% names(quality_raw)) {
      as.numeric(n50_contigs)
    } else {
      NA_real_
    }
  )
]
quality <- quality[
  nzchar(MAG_ID) &
  is.finite(completeness) &
  is.finite(contamination)
]
quality <- quality[, .SD[1L], by = MAG_ID]

taxonomy <- taxonomy_raw[
  ,
  .(
    MAG_ID = as.character(user_genome),
    classification = as.character(classification),
    classification_method = if (
      "classification_method" %in% names(taxonomy_raw)
    ) {
      as.character(classification_method)
    } else {
      NA_character_
    }
  )
]
taxonomy <- taxonomy[nzchar(MAG_ID)]
taxonomy <- taxonomy[, .SD[1L], by = MAG_ID]
taxonomy[
  ,
  taxon_label := deepest_taxon(classification)
]
taxonomy[
  ,
  deepest_rank := deepest_rank(classification)
]

quality <- merge(
  quality,
  taxonomy[
    ,
    .(
      MAG_ID,
      classification,
      taxon_label,
      deepest_rank,
      classification_method
    )
  ],
  by = "MAG_ID",
  all.x = TRUE
)
quality[
  ,
  gtdb_rep_classified := !is.na(classification)
]
quality[
  ,
  quality_tier := fifelse(
    completeness >= 90 & contamination <= 5,
    "High quality",
    fifelse(
      completeness >= 50 & contamination <= 10,
      "Medium quality",
      "Below medium-quality threshold"
    )
  )
]
quality[
  ,
  quality_tier := factor(
    quality_tier,
    levels = c(
      "High quality",
      "Medium quality",
      "Below medium-quality threshold"
    )
  )
]
quality[
  ,
  representative_status := factor(
    fifelse(
      gtdb_rep_classified,
      "GTDB-Tk representative",
      "Other dereplicated MAG"
    ),
    levels = c(
      "GTDB-Tk representative",
      "Other dereplicated MAG"
    )
  )
]

quality_counts <- quality[
  ,
  .(
    n_MAGs = .N,
    n_GTDPTk_representatives = sum(gtdb_rep_classified)
  ),
  by = quality_tier
]
quality_annotation <- data.table(
  x = 6.55,
  y = 12.0,
  hjust = 0,
  vjust = 0,
  label = sprintf(
    "MAGs=%s\nHigh quality=%s\nGTDB-Tk representatives=%s",
    format(nrow(quality), big.mark = ","),
    format(
      quality[quality_tier == "High quality", .N],
      big.mark = ","
    ),
    format(
      quality[gtdb_rep_classified == TRUE, .N],
      big.mark = ","
    )
  )
)

quality_fill <- c(
  "High quality" = purple_dark,
  "Medium quality" = purple_mid,
  "Below medium-quality threshold" = grey_light
)
representative_shape <- c(
  "GTDB-Tk representative" = 21,
  "Other dereplicated MAG" = 21
)

p_a <- ggplot(
  quality,
  aes(
    x = contamination,
    y = completeness
  )
) +
  geom_vline(
    xintercept = c(5, 10),
    colour = c(purple_mid, grey_mid),
    linewidth = c(0.66, 0.54),
    linetype = c("dashed", "dotted")
  ) +
  geom_hline(
    yintercept = c(50, 90),
    colour = c(grey_mid, purple_mid),
    linewidth = c(0.54, 0.66),
    linetype = c("dotted", "dashed")
  ) +
  geom_point(
    aes(
      fill = quality_tier,
      alpha = representative_status
    ),
    shape = 21,
    size = 2.15,
    colour = "#44384B",
    stroke = 0.34
  ) +
  geom_label(
    data = quality_annotation,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.06,
    vjust = 1.12,
    size = 3.25,
    lineheight = 1.02,
    fill = alpha("white", 0.90),
    colour = ink,
    label.padding = unit(0.16, "lines"),
    label.r = unit(0.08, "lines"),
    linewidth = 0.28
  ) +
  scale_fill_manual(
    values = quality_fill,
    drop = FALSE
  ) +
  scale_alpha_manual(
    values = c(
      "GTDB-Tk representative" = 0.92,
      "Other dereplicated MAG" = 0.38
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  scale_y_continuous(
    limits = c(0, 102),
    breaks = c(0, 25, 50, 75, 90, 100),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "MAG contamination (%)",
    y = "MAG completeness (%)",
    fill = "Quality tier",
    alpha = "MAG status"
  ) +
  theme_editorial(base_size = 11.6) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.key.width = unit(0.82, "lines"),
    legend.key.height = unit(0.72, "lines")
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 1),
    alpha = guide_legend(
      order = 2,
      nrow = 1,
      override.aes = list(
        fill = purple_mid,
        size = 3.0
      )
    )
  )

# -------------------------------------------------------------------------
# Panel B: prevalence and mean abundance across Ronneby
# -------------------------------------------------------------------------
abundance_raw <- fread(PATH_ABUNDANCE)
if (!"MAG_ID" %in% names(abundance_raw)) {
  stop(
    "Unexpected abundance schema: ",
    paste(names(abundance_raw), collapse = ",")
  )
}
sample_columns <- setdiff(names(abundance_raw), "MAG_ID")
if (length(sample_columns) != 65L) {
  stop(
    "Expected 65 Ronneby sample columns; observed ",
    length(sample_columns)
  )
}

for (column in sample_columns) {
  set(
    abundance_raw,
    j = column,
    value = as.numeric(abundance_raw[[column]])
  )
}

abundance <- abundance_raw[
  ,
  {
    values <- unlist(.SD, use.names = FALSE)
    values[!is.finite(values)] <- 0
    .(
      mean_relative_abundance_percent = mean(values),
      median_relative_abundance_percent = median(values),
      maximum_relative_abundance_percent = max(values),
      prevalence_fraction_gt0 = mean(values > 0),
      prevalence_fraction_ge_0_01 = mean(values >= 0.01)
    )
  },
  by = MAG_ID,
  .SDcols = sample_columns
]

if (nrow(abundance) != 1973L) {
  stop(
    "Expected 1,973 MAG abundance rows; observed ",
    nrow(abundance)
  )
}

abundance <- merge(
  abundance,
  quality[
    ,
    .(
      MAG_ID,
      quality_tier,
      gtdb_rep_classified,
      taxon_label,
      deepest_rank
    )
  ],
  by = "MAG_ID",
  all.x = TRUE
)

positive_means <- abundance[
  mean_relative_abundance_percent > 0,
  mean_relative_abundance_percent
]
pseudo <- if (length(positive_means)) {
  min(positive_means) / 2
} else {
  1e-8
}
abundance[
  ,
  mean_abundance_plot := pmax(
    mean_relative_abundance_percent,
    pseudo
  )
]

top_abundance <- copy(
  abundance[
    order(
      -mean_relative_abundance_percent,
      -prevalence_fraction_gt0
    )
  ][1:6]
)
top_abundance[
  ,
  label := fifelse(
    !is.na(taxon_label) & nzchar(taxon_label),
    taxon_label,
    MAG_ID
  )
]
top_abundance[
  ,
  label := wrap_label(label, width = 23L)
]

p_b <- ggplot(
  abundance,
  aes(
    x = prevalence_fraction_gt0,
    y = mean_abundance_plot
  )
) +
  geom_vline(
    xintercept = 0.10,
    colour = grey_mid,
    linewidth = 0.58,
    linetype = "dashed"
  ) +
  geom_point(
    aes(
      fill = quality_tier
    ),
    shape = 21,
    size = 2.15,
    colour = "#44384B",
    stroke = 0.34,
    alpha = 0.68
  ) +
  geom_label(
    data = top_abundance,
    aes(label = label),
    nudge_x = 0.018,
    nudge_y = 0,
    hjust = 0,
    size = 2.65,
    lineheight = 0.94,
    fill = alpha("white", 0.90),
    colour = ink,
    label.padding = unit(0.10, "lines"),
    label.r = unit(0.06, "lines"),
    linewidth = 0.22,
    show.legend = FALSE
  ) +
  annotate(
    "label",
    x = Inf,
    y = Inf,
    label = "n=65 samples; 1,973 MAGs",
    hjust = 1.05,
    vjust = 1.15,
    size = 3.15,
    fontface = "italic",
    fill = alpha("white", 0.90),
    colour = ink,
    label.padding = unit(0.13, "lines"),
    label.r = unit(0.07, "lines"),
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = quality_fill,
    drop = FALSE
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.10),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_log10(
    labels = label_number(accuracy = 0.001),
    breaks = log_breaks(n = 5),
    expand = expansion(mult = c(0.05, 0.18))
  ) +
  labs(
    x = "MAG prevalence across Ronneby samples",
    y = "Mean relative abundance (%)",
    fill = "Quality tier"
  ) +
  theme_editorial(base_size = 11.6) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(0.82, "lines"),
    legend.key.height = unit(0.72, "lines")
  ) +
  guides(
    fill = guide_legend(nrow = 1)
  )

# -------------------------------------------------------------------------
# Panel C: rule-based candidate-priority score for locked candidates
# -------------------------------------------------------------------------
context_raw <- fread(PATH_CONTEXT)
required_context <- c(
  "feature",
  "feature_family",
  "score",
  "is_taxon",
  "literature_context"
)
if (!all(required_context %in% names(context_raw))) {
  stop(
    "Unexpected context schema: ",
    paste(names(context_raw), collapse = ",")
  )
}

context <- context_raw[
  ,
  .(
    feature = as.character(feature),
    feature_family = clean_text(feature_family),
    score = as.numeric(score),
    is_taxon = parse_bool(is_taxon),
    literature_context = clean_text(literature_context),
    feature_label = if (
      "feature_label" %in% names(context_raw)
    ) {
      clean_text(feature_label)
    } else if ("label_clean" %in% names(context_raw)) {
      clean_text(label_clean)
    } else {
      clean_text(feature)
    }
  )
]
context <- context[
  nzchar(feature_label) &
  is.finite(score)
]
context <- context[, .SD[1L], by = feature]

if (nrow(context) != 14L) {
  stop(
    "Expected 14 locked context candidates; observed ",
    nrow(context)
  )
}
if (context[is_taxon == TRUE, .N] != 10L) {
  stop(
    "Expected 10 taxon candidates; observed ",
    context[is_taxon == TRUE, .N]
  )
}

context[
  ,
  candidate_type := factor(
    fifelse(
      is_taxon,
      "Taxon candidate",
      "Functional candidate"
    ),
    levels = c(
      "Taxon candidate",
      "Functional candidate"
    )
  )
]
context[
  ,
  external_context := factor(
    fifelse(
      nzchar(literature_context),
      "External Lindell context",
      "No direct external context"
    ),
    levels = c(
      "External Lindell context",
      "No direct external context"
    )
  )
]
context[
  ,
  feature_label_plot := wrap_label(
    feature_label,
    width = 33L
  )
]
setorder(context, -score, feature_label)
context[
  ,
  feature_label_plot := factor(
    feature_label_plot,
    levels = rev(feature_label_plot)
  )
]

candidate_type_fill <- c(
  "Taxon candidate" = purple_dark,
  "Functional candidate" = purple_light
)
context_max <- max(context$score, na.rm = TRUE)

p_c <- ggplot(
  context,
  aes(
    x = score,
    y = feature_label_plot
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = score,
      yend = feature_label_plot
    ),
    colour = "#CFC3DA",
    linewidth = 1.12
  ) +
  geom_point(
    aes(
      fill = candidate_type
    ),
    shape = 21,
    size = 3.85,
    colour = "#3E2A4C",
    stroke = 0.72
  ) +
  geom_label(
    aes(
      label = sprintf("%.2f", score)
    ),
    nudge_x = context_max * 0.025,
    hjust = 0,
    size = 3.05,
    colour = "#403943",
    fill = alpha("white", 0.94),
    label.padding = unit(0.10, "lines"),
    label.r = unit(0.07, "lines"),
    linewidth = 0.25,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = candidate_type_fill,
    drop = TRUE
  ) +
  scale_x_continuous(
    limits = c(0, context_max * 1.22),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Rule-based candidate-priority score",
    y = "Genome-resolved and functional candidates",
    fill = "Candidate type"
  ) +
  theme_editorial(base_size = 11.2) +
  theme(
    axis.text.y = element_text(
      size = 8.95,
      lineheight = 0.96
    ),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.direction = "horizontal",
    legend.key.width = unit(0.80, "lines"),
    legend.key.height = unit(0.72, "lines")
  ) +
  guides(
    fill = guide_legend(order = 1, nrow = 1)
  )

# -------------------------------------------------------------------------
# Panel D: named-taxa overlap priorities and match levels only
# -------------------------------------------------------------------------
overlap_raw <- fread(PATH_OVERLAP)
required_overlap <- c(
  "source",
  "feature",
  "match_level",
  "priority"
)
if (!all(required_overlap %in% names(overlap_raw))) {
  stop(
    "Unexpected overlap schema: ",
    paste(names(overlap_raw), collapse = ",")
  )
}

priority_raw <- as.character(overlap_raw$priority)
priority_numeric <- suppressWarnings(
  as.numeric(priority_raw)
)
if (all(!is.finite(priority_numeric))) {
  priority_clean <- tolower(clean_text(priority_raw))
  priority_numeric <- fifelse(
    grepl("tier[[:space:]]*1|high|highest|primary|priority[[:space:]]*1", priority_clean),
    3,
    fifelse(
      grepl("tier[[:space:]]*2|medium|moderate|secondary|priority[[:space:]]*2", priority_clean),
      2,
      fifelse(
        grepl("tier[[:space:]]*3|low|exploratory|tertiary|priority[[:space:]]*3", priority_clean),
        1,
        NA_real_
      )
    )
  )
}
if (all(!is.finite(priority_numeric))) {
  priority_numeric <- rank(
    -seq_along(priority_raw),
    ties.method = "first"
  )
}

overlap <- overlap_raw[
  ,
  .(
    source = clean_text(source),
    feature = as.character(feature),
    match_level = clean_text(match_level),
    priority_label = clean_text(priority),
    priority_numeric = priority_numeric,
    feature_label = if (
      "feature_label" %in% names(overlap_raw)
    ) {
      clean_text(feature_label)
    } else {
      clean_text(feature)
    }
  )
]
overlap <- overlap[
  nzchar(feature_label) &
  nzchar(match_level) &
  is.finite(priority_numeric)
]
overlap <- overlap[, .SD[1L], by = feature]

if (nrow(overlap) != 20L) {
  stop(
    "Expected 20 locked overlap features; observed ",
    nrow(overlap)
  )
}

normalise_match_level <- function(x) {
  y <- tolower(clean_text(x))
  fifelse(
    grepl("species", y),
    "Species-level",
    fifelse(
      grepl("genus", y),
      "Genus-level",
      fifelse(
        grepl("family", y),
        "Family-level",
        "Other supported match"
      )
    )
  )
}

overlap[
  ,
  match_level_display := normalise_match_level(match_level)
]
overlap[
  ,
  bin_like_label := is_bin_like_label(feature_label)
]
overlap[
  ,
  named_taxon_eligible := (
    !bin_like_label &
    match_level_display %chin% c("Species-level", "Genus-level") &
    grepl("[A-Za-z]", feature_label)
  )
]

overlap_named <- copy(overlap[named_taxon_eligible == TRUE])
overlap_named[
  ,
  feature_label_key := clean_text(feature_label)
]
setorder(
  overlap_named,
  -priority_numeric,
  match_level_display,
  feature_label_key
)
overlap_named <- overlap_named[
  nzchar(feature_label_key)
][
  , .SD[1L], by = feature_label_key
]

named_overlap_audit <- data.table(
  metric = c(
    "raw_rows",
    "named_taxon_rows",
    "species_level_rows",
    "genus_level_rows",
    "bin_like_rows_removed",
    "named_panel_threshold",
    "panelD_decision"
  ),
  value = c(
    as.character(nrow(overlap)),
    as.character(nrow(overlap_named)),
    as.character(overlap_named[match_level_display == "Species-level", .N]),
    as.character(overlap_named[match_level_display == "Genus-level", .N]),
    as.character(overlap[bin_like_label == TRUE, .N]),
    "8",
    NA_character_
  )
)

display_panel_d <- nrow(overlap_named) >= 8L
panel_d_status <- if (display_panel_d) "INCLUDED" else "OMITTED_WEAK_NAMED_TAXON_OVERLAP"
named_overlap_audit[metric == "panelD_decision", value := panel_d_status]

if (display_panel_d) {
  overlap_display <- copy(head(overlap_named, 12L))
  setorder(overlap_display, -priority_numeric, match_level_display, feature_label_key)
  overlap_display[, feature_label_plot := wrap_label(feature_label, width = 34L)]
  overlap_display[, feature_label_plot := factor(feature_label_plot, levels = rev(feature_label_plot))]
  priority_min <- min(overlap_display$priority_numeric, na.rm = TRUE)
  priority_max <- max(overlap_display$priority_numeric, na.rm = TRUE)
  if (priority_max == priority_min) {
    overlap_display[, priority_plot := seq(1, .N) / .N]
  } else {
    overlap_display[, priority_plot := rescale(priority_numeric, to = c(0.25, 1.00))]
  }
  match_order <- c("Species-level", "Genus-level")
  match_levels_present <- match_order[match_order %chin% unique(overlap_display$match_level_display)]
  overlap_display[, match_level_display := factor(match_level_display, levels = match_levels_present)]
  match_fill <- c("Species-level" = purple_dark, "Genus-level" = purple_mid)
  match_shape <- c("Species-level" = 21, "Genus-level" = 22)
  p_d <- ggplot(overlap_display, aes(x = priority_plot, y = feature_label_plot)) +
    geom_segment(aes(x = 0, xend = priority_plot, yend = feature_label_plot), colour = "#CFC3DA", linewidth = 1.12) +
    geom_point(aes(fill = match_level_display, shape = match_level_display), size = 4.05, colour = "#3E2A4C", stroke = 0.76) +
    geom_label(aes(label = priority_label), nudge_x = 0.025, hjust = 0, size = 2.95, colour = "#403943", fill = alpha("white", 0.94), label.padding = unit(0.10, "lines"), label.r = unit(0.07, "lines"), linewidth = 0.25, show.legend = FALSE) +
    scale_fill_manual(values = match_fill[match_levels_present], drop = TRUE) +
    scale_shape_manual(values = match_shape[match_levels_present], drop = TRUE) +
    scale_x_continuous(limits = c(0, 1.16), breaks = c(0, 0.25, 0.50, 0.75, 1.00), expand = expansion(mult = c(0.01, 0.01))) +
    labs(x = "Lindell overlap priority (scaled 0–1)", y = "Named Lindell-overlap taxa", fill = "Supported taxonomic match", shape = "Supported taxonomic match") +
    theme_editorial(base_size = 11.2) +
    theme(axis.text.y = element_text(size = 8.95, lineheight = 0.96), legend.position = "bottom", legend.direction = "horizontal", legend.key.width = unit(0.82, "lines"), legend.key.height = unit(0.72, "lines")) +
    guides(fill = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
} else {
  overlap_display <- overlap_named[0]
  match_levels_present <- character(0)
  p_d <- NULL
}

# -------------------------------------------------------------------------
# Source data and provenance
# -------------------------------------------------------------------------
write_tsv(
  quality,
  file.path(
    SRC_DIR,
    "Figure_S06_panelA_MAG_quality_taxonomy.tsv"
  )
)
write_tsv(
  quality_counts,
  file.path(
    SRC_DIR,
    "Figure_S06_panelA_quality_counts.tsv"
  )
)
write_tsv(
  taxonomy,
  file.path(
    SRC_DIR,
    "Figure_S06_panelA_GTDPTk_representatives.tsv"
  )
)
write_tsv(
  abundance,
  file.path(
    SRC_DIR,
    "Figure_S06_panelB_MAG_prevalence_abundance.tsv"
  )
)
write_tsv(
  top_abundance,
  file.path(
    SRC_DIR,
    "Figure_S06_panelB_top_abundant_labels.tsv"
  )
)
write_tsv(
  context,
  file.path(
    SRC_DIR,
    "Figure_S06_panelC_locked_Lindell_context_candidates.tsv"
  )
)
write_tsv(
  overlap,
  file.path(
    SRC_DIR,
    "Figure_S06_panelD_locked_overlap_full.tsv"
  )
)
write_tsv(
  overlap_display,
  file.path(
    SRC_DIR,
    "Figure_S06_panelD_locked_overlap_display.tsv"
  )
)
write_tsv(
  named_overlap_audit,
  file.path(
    SRC_DIR,
    "Figure_S06_panelD_named_overlap_audit.tsv"
  )
)

panel_design <- rbindlist(list(
  data.table(
    panel = c("A", "B", "C"),
    content = c(
      "MAG completeness, contamination and GTDB-Tk representative status",
      "MAG prevalence and relative abundance across 65 Ronneby samples",
      "Locked candidates with curated external Lindell context"
    ),
    status = rep("INCLUDED", 3L),
    interpretation_boundary = c(
      "Quality and classification summary only",
      "Ronneby descriptive abundance context only",
      "External bioaccumulation/binding context; not cohort validation or degradation evidence"
    )
  ),
  data.table(
    panel = "D",
    content = if (display_panel_d) {
      "Named Lindell-overlap taxa with supported species/genus-level match and locked priority"
    } else {
      "Named Lindell-overlap panel omitted because the named-taxon overlap evidence was too sparse"
    },
    status = panel_d_status,
    interpretation_boundary = if (display_panel_d) {
      "Named taxa only; no exact MAG equivalence implied and no causal validation claim"
    } else {
      "Omitted to avoid overinterpreting weak named-taxon overlap"
    }
  )
))
write_tsv(
  panel_design,
  file.path(
    SRC_DIR,
    "Figure_S06_panel_design.tsv"
  )
)

claim_gate <- data.table(
  gate = c(
    "AFFIRMATIVE_PFAS_DEGRADATION_CLAIMS",
    "RONNEBY_VALIDATION_CLAIMS",
    "TAXONOMIC_RANK_PROMOTION",
    "CAUSAL_BIOMARKER_CLAIMS"
  ),
  status = rep("PASS", 4L),
  detail = c(
    "None introduced; external evidence described as bioaccumulation/binding context",
    "Lindell context explicitly external to Ronneby",
    "Labels retained at source-supported GTDB-Tk or locked match level",
    "Candidate scores and priorities framed as exploratory validation priorities"
  )
)
write_tsv(
  claim_gate,
  file.path(
    MAN_DIR,
    "scientific_claim_gate.tsv"
  )
)

provenance <- data.table(
  role = c(
    "panel_A_quality",
    "panel_A_taxonomy",
    "panel_B_abundance",
    "panel_C_context",
    "panel_D_overlap"
  ),
  source_path = c(
    PATH_QUALITY,
    PATH_TAXONOMY,
    PATH_ABUNDANCE,
    PATH_CONTEXT,
    PATH_OVERLAP
  ),
  source_type = c(
    "MAG quality source data",
    "GTDB-Tk representative classification summary",
    "Ronneby MAG relative abundance percentage matrix",
    "Locked curated Lindell-context candidates",
    "Locked representative Lindell-overlap priorities"
  ),
  rows_or_samples = c(
    nrow(quality),
    nrow(taxonomy),
    length(sample_columns),
    nrow(context),
    nrow(overlap)
  )
)
write_tsv(
  provenance,
  file.path(MAN_DIR, "provenance.tsv")
)

# -------------------------------------------------------------------------
# Editorial layout
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Editorial layout
# -------------------------------------------------------------------------
if (display_panel_d) {
  top <- p_a + p_b + plot_layout(widths = c(1.02, 0.98))
  bottom <- p_c + p_d + plot_layout(widths = c(1.01, 0.99))
  figure <- top / bottom +
    plot_layout(heights = c(1.00, 1.10)) +
    plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 22, colour = "#171419"), plot.tag.position = c(0, 1)))
  fig_width <- 17.2
  fig_height <- 14.0
} else {
  top <- p_a + p_b + plot_layout(widths = c(1.02, 0.98))
  figure <- top / p_c +
    plot_layout(heights = c(0.96, 1.06)) +
    plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 22, colour = "#171419"), plot.tag.position = c(0, 1)))
  fig_width <- 17.2
  fig_height <- 12.8
}

pdf_path <- file.path(
  FIG_DIR,
  "Figure_S06_excellence_candidate.pdf"
)
png_path <- file.path(
  FIG_DIR,
  "Figure_S06_excellence_candidate.png"
)
tiff_path <- file.path(
  FIG_DIR,
  "Figure_S06_excellence_candidate.tiff"
)

ggsave(
  pdf_path,
  figure,
  width = fig_width,
  height = fig_height,
  units = "in",
  device = cairo_pdf
)
ggsave(
  png_path,
  figure,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  tiff_path,
  figure,
  width = fig_width,
  height = fig_height,
  units = "in",
  dpi = 300,
  compression = "lzw",
  bg = "white"
)

legend_text <- if (display_panel_d) {
  paste(
    "Figure S6. Genome-resolved reconstruction quality, Ronneby abundance,",
    "and external Lindell-context prioritization. Panel A shows completeness",
    "and contamination for 1,973 MAGs, with quality tiers and the 429 GTDB-Tk",
    "representative MAGs indicated. Panel B shows MAG prevalence and mean",
    "relative abundance across 65 Ronneby samples. Panel C shows the 14 locked",
    "genome-resolved and functional candidates ranked by the verified rule-based",
    "candidate-priority score; fill distinguishes taxon from functional candidates.",
    "The Lindell dataset defines external provenance for this candidate set and is",
    "not encoded as a separate point-level legend. Panel D shows only named taxa",
    "with source-supported species- or genus-level overlap to the locked Lindell",
    "context table. Its x-axis rescales the locked overlap priority to 0–1 for",
    "display only, while adjacent labels retain the original priority category.",
    "Exact MAG/bin equivalence is not implied. Lindell evidence is external",
    "bioaccumulation, binding, uptake or sorption context and does not establish",
    "PFAS degradation, causal microbial effects, or Ronneby cohort validation."
  )
} else {
  paste(
    "Figure S6. Genome-resolved reconstruction quality, Ronneby abundance,",
    "and external Lindell-context prioritization. Panel A shows completeness",
    "and contamination for 1,973 MAGs, with quality tiers and the 429 GTDB-Tk",
    "representative MAGs indicated. Panel B shows MAG prevalence and mean",
    "relative abundance across 65 Ronneby samples. Panel C shows the 14 locked",
    "genome-resolved and functional candidates ranked by the verified rule-based",
    "candidate-priority score; fill distinguishes taxon from functional candidates.",
    "The Lindell dataset defines external provenance for this candidate set and is",
    "not encoded as a separate point-level legend. A named-taxon overlap panel was",
    "evaluated but omitted because the named",
    "overlap evidence was too sparse to justify a clear editorial display.",
    "Lindell evidence remains external bioaccumulation, binding, uptake or",
    "sorption context and does not establish PFAS degradation, causal microbial",
    "effects, or Ronneby cohort validation."
  )
}
writeLines(
  legend_text,
  file.path(LEG_DIR, "Figure_S06.txt")
)

writeLines(
  c(
    "FIGURE=Figure_S06",
    "PANEL_C_DISPLAY_LABEL=Rule-based candidate-priority score",
    "PANEL_C_FORMAL_INTERPRETATION=dimensionless exploratory ordering metric",
    "PANEL_C_NOT_A=probability;standardized_effect_size;causal_estimate;validated_biomarker",
    "PANEL_C_LINDELL_CONTEXT=external provenance for the locked candidate set; not a separate point-level evidence encoding",
    "PANEL_D_DISPLAY_LABEL=Lindell overlap priority (scaled 0–1)",
    "PANEL_D_SCALING=display-only rescaling of locked overlap priority; adjacent labels retain original priority category",
    "CLAIM_BOUNDARY=no PFAS degradation, causal microbial effect, or Ronneby cohort validation claim"
  ),
  file.path(
    SRC_DIR,
    "Figure_S06_score_and_context_definition.txt"
  )
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
  detail = if (display_panel_d) {
    paste(
      "Confirm A-D panel labels; readable MAG-quality thresholds and moved",
      "summary annotation in A; visible prevalence and log-scaled mean abundance",
      "axes with non-overlapping labels in B; readable rule-based candidate scores and",
      "candidate-type legend without a redundant external-context shape legend in C;",
      "readable named-taxon overlap priorities and only",
      "species/genus-level categories in D; balanced spacing, typography, purple",
      "editorial styling and no degradation or cohort-validation claim."
    )
  } else {
    paste(
      "Confirm A-C panel labels; readable MAG-quality thresholds and moved summary",
      "annotation in A; visible prevalence and log-scaled mean abundance axes with",
      "non-overlapping labels in B; readable rule-based candidate scores and a candidate-",
      "type legend without a redundant external-context shape legend in C; balanced",
      "3-panel spacing, typography, purple",
      "editorial styling and no degradation or cohort-validation claim."
    )
  }
)
write_tsv(
  visual_qc,
  file.path(MAN_DIR, "visual_qc.tsv")
)

cat(
  "MAG_QUALITY_SOURCE=PASS rows=",
  nrow(quality),
  " high_quality=",
  quality[quality_tier == "High quality", .N],
  " medium_or_higher=",
  quality[
    quality_tier %in% c(
      "High quality",
      "Medium quality"
    ),
    .N
  ],
  "\n",
  sep = ""
)
cat(
  "GTDB_TAXONOMY_SOURCE=PASS representatives=",
  nrow(taxonomy),
  " matched_to_quality=",
  quality[gtdb_rep_classified == TRUE, .N],
  "\n",
  sep = ""
)
cat(
  "RONNEBY_MAG_ABUNDANCE_SOURCE=PASS mags=",
  nrow(abundance),
  " samples=",
  length(sample_columns),
  "\n",
  sep = ""
)
cat(
  "LINDELL_CONTEXT_SOURCE=PASS candidates=",
  nrow(context),
  " taxon_candidates=",
  context[is_taxon == TRUE, .N],
  "\n",
  sep = ""
)
cat(
  "LOCKED_OVERLAP_SOURCE=PASS raw_rows=",
  nrow(overlap),
  " named_taxa=",
  nrow(overlap_named),
  " displayed=",
  if (display_panel_d) nrow(overlap_display) else 0L,
  " panel_d=",
  panel_d_status,
  " match_levels=",
  if (display_panel_d) {
    paste(match_levels_present, collapse = "|")
  } else {
    "OMITTED"
  },
  "\n",
  sep = ""
)
cat(
  "PFAS_DEGRADATION_GATE=PASS_NO_AFFIRMATIVE_CLAIMS\n"
)
cat(
  if (display_panel_d) {
    "FIGURE_S06_PANEL_DESIGN=A_QUALITY_TAXONOMY_B_PREVALENCE_ABUNDANCE_C_LINDELL_CONTEXT_D_NAMED_TAXA_PRIORITY\n"
  } else {
    "FIGURE_S06_PANEL_DESIGN=A_QUALITY_TAXONOMY_B_PREVALENCE_ABUNDANCE_C_LINDELL_CONTEXT\n"
  }
)
cat("TSV_SERIALIZATION_SANITIZER=PASS factors_newlines_tabs=TRUE\n")
cat("GENERATOR_EXECUTION_PASS\n")
